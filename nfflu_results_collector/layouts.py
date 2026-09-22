"""Path templates for stage-centric and sample-centric nf-flu outputs."""
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
    'bcftools_consensus':  '{sample}/consensus/bcftools/{sample}.consensus.fasta',
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
    """Detect the sequencing platform, defaulting to Illumina."""
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
    """Detect the output layout, or return None if neither layout matches."""
    for name in PATHS_BY_LAYOUT:
        if glob.glob(output_path(outdir, 'mapping_dir', layout=name)):
            return name
    return None


def output_path(outdir, output_name, sample='*', platform='illumina', layout=None):
    """Build an output path, using a glob when `sample` is `*`."""
    template = PATHS_BY_LAYOUT[layout or DEFAULT_LAYOUT][output_name]
    return os.path.join(outdir, template.format(sample=sample, blast=BLAST_DIR_BY_PLATFORM[platform]))


def find_by_sample(outdir, output_name, platform='illumina', layout=None):
    """Return `(sample_name, path)` matches for whole-segment sample fields."""
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
