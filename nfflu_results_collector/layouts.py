"""Where nf-flu publishes its outputs.

Every path read out of, or written into, an nf-flu output directory is declared
here as a template with a `{sample}` field, so changing the layout is a change
to this file alone. auto-nfflu imports this module too.

Two layouts are described. STAGE_CENTRIC groups outputs by pipeline stage
(`<outdir>/mapping/<sample>/`) and is what nf-flu 3.10 publishes. SAMPLE_CENTRIC
groups them by sample (`<outdir>/<sample>/mapping/`) and matches the per-sample
publishDir config in auto-nfflu's assets/config/nf-flu_nextflow_sample.config.
Callers select one by name; DEFAULT_LAYOUT applies when they do not.

One template serves both a whole-run glob and a single-sample path, since
`sample` defaults to the `*` wildcard:

    output_path(outdir, 'mapping_dir')               -> <outdir>/mapping/*
    output_path(outdir, 'mapping_dir', sample='S1')  -> <outdir>/mapping/S1

Under STAGE_CENTRIC the two BLAST paths differ between sequencing platforms, so
those templates carry a `{blast}` field filled from the platform. The per-sample
config publishes both platforms to one path, so SAMPLE_CENTRIC spells them out.

`nextclade_dir` and `pipeline_status` are auto-nfflu's own outputs rather than
nf-flu's, published alongside the run's nf-flu results.
"""
import glob
import json
import logging
import os
import re


STAGE_CENTRIC = {
    'mapping_dir':         'mapping/{sample}',
    'sample_bams':         'mapping/{sample}/{sample}*.bam',
    'pileup_dir':          'pileups/{sample}',
    'read_counts':         'irma/{sample}/tables/READ_COUNTS.txt',
    'mixtures_csv':        'mixtures/{sample}/{sample}_mixtures.csv',
    'mixtures_txt':        'mixtures/{sample}/{sample}_mixtures.txt',
    'irma_consensus':      'consensus/irma/{sample}.irma.consensus.fasta',
    'blastn_ref':          '{blast}/irma/{sample}*blastn.txt',
    'reference_sequences': 'reference_sequences/{sample}',
    'variants':            'variants/{sample}',
    'bcftools_consensus':  'consensus/bcftools/{sample}.consensus.fasta',
    'blastn_subtype':      '{blast}/consensus/{sample}*.blastn.txt',
    'fastq_dir':           'fastq',
    'subtype_results':     'subtyping_report/subtype_results.csv',
    'idxstats':            'mapping/{sample}/{sample}*.idxstats',
    'cleavage':            'annotation/{sample}/{sample}.cleavage.tsv',
    'genoflu':             'genoflu/{sample}.genoflu.tsv',
    'nextclade_tsvs':      'nextclade/{sample}',
    'software_versions':   'pipeline_info/software_versions.yml',
    'nextclade_dir':       'nextclade',
    'pipeline_status':     'pipeline_status.csv',
}

SAMPLE_CENTRIC = {
    'mapping_dir':         '{sample}/mapping',
    'sample_bams':         '{sample}/mapping/{sample}*.bam',
    'pileup_dir':          '{sample}/pileups',
    'read_counts':         '{sample}/irma/tables/READ_COUNTS.txt',
    'mixtures_csv':        '{sample}/mixtures/{sample}_mixtures.csv',
    'mixtures_txt':        '{sample}/mixtures/{sample}_mixtures.txt',
    'irma_consensus':      '{sample}/consensus/irma/{sample}*.irma.consensus.fasta',
    'blastn_ref':          '{sample}/blast/irma/{sample}*blastn.txt',
    'reference_sequences': '{sample}/reference_sequences/*',
    'variants':            '{sample}/variants/*',
    'bcftools_consensus':  '{sample}/consensus/bcftools/{sample}*.consensus.fasta',
    'blastn_subtype':      '{sample}/blast/consensus/{sample}*.blastn.txt',
    'fastq_dir':           '{sample}/fastq',
    'subtype_results':     'aggregate/bcftools/subtyping_report/subtype_results.csv',
    'idxstats':            '{sample}/mapping/{sample}*.idxstats',
    'cleavage':            '{sample}/annotation/{sample}.cleavage.tsv',
    'genoflu':             '{sample}/genoflu/{sample}.genoflu.tsv',
    'nextclade_tsvs':      '{sample}/nextclade',
    'software_versions':   'pipeline_info/software_versions.yml',
    'nextclade_dir':       'aggregate/nextclade',
    'pipeline_status':     'pipeline_status.csv',
}

PATHS_BY_LAYOUT = {'stage': STAGE_CENTRIC, 'sample': SAMPLE_CENTRIC}

DEFAULT_LAYOUT = 'stage'

# nf-flu publishes BLAST output one directory deeper on Illumina than on
# nanopore; see conf/modules_illumina.config and conf/modules_nanopore.config.
BLAST_DIR_BY_PLATFORM = {'illumina': os.path.join('blast', 'blastn'), 'nanopore': 'blast'}


def detect_platform(outdir):
    """Return 'nanopore' or 'illumina' for an nf-flu output directory.

    nf-flu writes pipeline_info/samplesheet.fixed.csv under both platforms, but
    with different headers: 'sample,barcode' for nanopore and
    'sample,fastq1,fastq2,single_end' for Illumina. Falls back to 'illumina'
    when the file is missing or unrecognised.
    """
    samplesheet_path = os.path.join(outdir, 'pipeline_info', 'samplesheet.fixed.csv')

    try:
        with open(samplesheet_path, 'r') as f:
            header = f.readline().strip()
    except OSError:
        logging.warning(json.dumps({
            "event_type": "platform_samplesheet_not_found",
            "samplesheet_path": samplesheet_path,
        }))
        return 'illumina'

    platform = 'nanopore' if 'barcode' in [field.strip() for field in header.split(',')] else 'illumina'

    logging.info(json.dumps({"event_type": "platform_detected", "platform": platform, "header": header}))

    return platform


def detect_layout(outdir):
    """Which layout an nf-flu output directory actually uses, or None when
    neither matches. `mapping_dir` is 'mapping/*' under stage and '*/mapping'
    under sample, so at most one can match a given directory."""
    for name in PATHS_BY_LAYOUT:
        if glob.glob(output_path(outdir, 'mapping_dir', layout=name)):
            return name
    return None


def output_path(outdir, output_name, sample='*', platform='illumina', layout=None):
    """Absolute path (or glob, when `sample` is left as the wildcard) for one
    nf-flu output within `outdir`."""
    template = PATHS_BY_LAYOUT[layout or DEFAULT_LAYOUT][output_name]
    return os.path.join(outdir, template.format(sample=sample, blast=BLAST_DIR_BY_PLATFORM[platform]))


def find_by_sample(outdir, output_name, platform='illumina', layout=None):
    """Every match of `output_name` across all samples, as (sample_name, path) pairs.

    The sample name is recovered from the same template that produced the glob,
    because where it sits in the path is itself layout-dependent -- under
    STAGE_CENTRIC it names a directory inside mapping/, under SAMPLE_CENTRIC it
    names the directory containing it.

    Only defined for outputs whose first `{sample}` is a whole path segment.
    Where it sits alongside a wildcard, as in '{sample}*blastn.txt', the split
    between the two is ambiguous and the recovered name would be unreliable.
    """
    template = PATHS_BY_LAYOUT[layout or DEFAULT_LAYOUT][output_name]

    if '{sample}' not in template.split(os.sep):
        raise ValueError(f"{output_name!r} does not carry the sample name as a whole path segment; use output_path() instead")

    full_template = os.path.join(outdir, template).replace('{blast}', BLAST_DIR_BY_PLATFORM[platform])
    literal_parts = [re.escape(part).replace(r'\*', '[^/]*') for part in full_template.split('{sample}')]
    sample_matcher = re.compile('([^/]+)'.join(literal_parts) + r'/?$')

    found = []
    for found_path in glob.glob(output_path(outdir, output_name, platform=platform, layout=layout)):
        match = sample_matcher.match(found_path)
        if match:
            found.append((match.group(1), found_path))

    return sorted(found)
