"""Builds a synthetic nf-flu analysis output directory for tests.

Pure-Python (no pytest import) so it can be reused both by the pytest
fixture in conftest.py and by tests/regenerate_golden.py.

Four samples exercise the sample-ID parsing branches:
  - SAMPLE_STANDARD: standard "CID-Plate-Index-Well" format
  - SAMPLE_CONTROL:  control format, and partial/missing data (only
    2 of 8 segments sequenced, no subtype call, no cleavage, no genoflu)
  - SAMPLE_SALVAGE:  fails both standard/control regexes (lowercase CID),
    falls all the way through to the whole-name CID fallback
  - SAMPLE_EXTERNAL: external partner batch, whose ID contains underscores
    (auto-nfflu's sample_id_before_first_period) -- these must survive
    sequence-ID parsing, which splits on underscores
"""
import os
import textwrap

RUN_ID = "250101_M00123_0001_000000000-ABCDE"
ANALYSIS_TYPE = "short"

SAMPLE_STANDARD = "C0000000001-2024-A1-A01"
SAMPLE_CONTROL = "NTC20240101-Control-2024-B2-B02"
SAMPLE_SALVAGE = "c0000000003-2024-D4-D04"
SAMPLE_EXTERNAL = "PARTNER_BATCH7_0004"

# Sorted, because auto-discovery (tools.collect_nfflu_fastq_names) returns
# sorted names and tests compare that path against this list.
SAMPLE_IDS = [SAMPLE_STANDARD, SAMPLE_CONTROL, SAMPLE_EXTERNAL, SAMPLE_SALVAGE]

SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
SEGMENT_NUMBER = {seg: i + 1 for i, seg in enumerate(SEGMENTS)}
REF_SEG_LENGTH = {
    "PB2": 2341, "PB1": 2341, "PA": 2233, "HA": 1701,
    "NP": 1565, "NA": 1413, "M": 1027, "NS": 890,
}

# (segment -> N-count out of a 100bp synthetic consensus) per sample.
# Segments omitted from a sample's dict are absent from its FASTA entirely.
CONSENSUS_N_COUNTS = {
    SAMPLE_STANDARD: {
        "PB2": 1, "PB1": 1, "PA": 2, "HA": 3, "NP": 1, "NA": 4, "M": 0, "NS": 1,
    },
    SAMPLE_CONTROL: {
        "HA": 40, "NA": 7,
    },
    SAMPLE_SALVAGE: {
        "PB2": 5, "PB1": 6, "PA": 7, "HA": 2, "NP": 9, "NA": 3, "M": 15, "NS": 8,
    },
    SAMPLE_EXTERNAL: {
        "PB2": 3, "PB1": 4, "PA": 6, "HA": 1, "NP": 2, "NA": 5, "M": 11, "NS": 4,
    },
}

READS_MAPPED = {
    SAMPLE_STANDARD: {
        "PB2": 52000, "PB1": 51000, "PA": 49000, "HA": 61000,
        "NP": 58000, "NA": 55000, "M": 72000, "NS": 68000,
    },
    SAMPLE_CONTROL: {
        "HA": 1200, "NA": 900,
    },
    SAMPLE_SALVAGE: {
        "PB2": 41000, "PB1": 40500, "PA": 38000, "HA": 47000,
        "NP": 44000, "NA": 42000, "M": 59000, "NS": 51000,
    },
    SAMPLE_EXTERNAL: {
        "PB2": 33000, "PB1": 32500, "PA": 31000, "HA": 39000,
        "NP": 36000, "NA": 34000, "M": 48000, "NS": 43000,
    },
}

NEXTCLADE_DATASET = {
    SAMPLE_STANDARD: ("nextstrain/flu/h5nx/ha/EPI2101974", "2025-11-10--12-00-00Z"),
    SAMPLE_CONTROL: ("nextstrain/flu/h5nx/ha/EPI2101974", "2025-11-10--12-00-00Z"),
    SAMPLE_SALVAGE: ("nextstrain/flu/h3n2/ha/EPI1857216", "2025-11-04--15-46-13Z"),
    SAMPLE_EXTERNAL: ("nextstrain/flu/h5nx/ha/EPI2101974", "2025-11-10--12-00-00Z"),
}


def _write(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(content)


def _sample_path(analysis_dir, sample, *parts):
    """nf-flu publishes per-sample outputs under <analysis_dir>/<sample>/."""
    return os.path.join(analysis_dir, sample, *parts)


def _sequence_id(sample, segment):
    """Consensus FASTA header and Nextclade seqName, as nf-flu writes them:
    "<sample>_<n>_<GENE>" (bcftools.nf sets sequenceID = "${sample}_${segment}",
    and segment labels are "4_HA", "7_M", ... -- see IAV_SEGMENT_NAMES in
    nf-flu's bin/subtyping_report.py)."""
    return f"{sample}_{SEGMENT_NUMBER[segment]}_{segment}"


def build_fixture(root):
    """Build the fixture tree under `root` (a Path). Returns the nf-flu
    output directory path (i.e. the `analysis_dir` to pass to the collector)."""
    root = str(root)

    run_analysis_dir = os.path.join(root, "analysis_output", RUN_ID, ANALYSIS_TYPE)
    analysis_dir = os.path.join(run_analysis_dir, "nf-flu-3.10-output")

    _build_samplesheets(run_analysis_dir)
    _build_fastq(analysis_dir)
    _build_subtyping_report(analysis_dir)
    _build_mapping(analysis_dir)
    _build_consensus(analysis_dir)
    _build_annotation(analysis_dir)
    _build_genoflu(analysis_dir)
    _build_nextclade(analysis_dir)
    _build_mixtures(analysis_dir)
    _build_software_versions(analysis_dir)

    return analysis_dir


def _build_samplesheets(run_analysis_dir):
    path = os.path.join(run_analysis_dir, "samplesheets", f"{RUN_ID}_start_samplesheet.csv")
    lines = ["ID,R1,R2"]
    for sample in SAMPLE_IDS:
        lines.append(f"{sample},/dev/null/{sample}_R1.fastq.gz,/dev/null/{sample}_R2.fastq.gz")
    _write(path, "\n".join(lines) + "\n")


def _build_fastq(analysis_dir):
    for sample in SAMPLE_IDS:
        for read in ("1", "2"):
            path = _sample_path(analysis_dir, sample, "fastq", f"{sample}_{read}.merged.fastq.gz")
            _write(path, "")


def _build_subtyping_report(analysis_dir):
    # S1 -> H5N1, S3 -> H3N2, S4 -> H5N1. S2 (control) intentionally has no
    # subtype call.
    path = os.path.join(analysis_dir, "aggregate", "bcftools", "subtyping_report", "subtype_results.csv")
    content = textwrap.dedent(f"""\
        ,sample,Genotype,H_type,N_type
        0,{SAMPLE_STANDARD},3,5,1
        1,{SAMPLE_SALVAGE},7,3,2
        2,{SAMPLE_EXTERNAL},3,5,1
    """)
    _write(path, content)


def _build_mapping(analysis_dir):
    for sample, seg_reads in READS_MAPPED.items():
        for seg, reads in seg_reads.items():
            seg_num = SEGMENT_NUMBER[seg]
            ref_name = f"Segment_{seg_num}_{seg}"
            filename = f"{sample}.Segment_{seg_num}_{seg}.idxstats"
            path = _sample_path(analysis_dir, sample, "mapping", filename)
            unmapped = max(10, reads // 200)
            content = f"{ref_name}\t{REF_SEG_LENGTH[seg]}\t{reads}\t{unmapped}\n*\t0\t0\t{unmapped}\n"
            _write(path, content)


def _build_consensus(analysis_dir):
    for sample, seg_ns in CONSENSUS_N_COUNTS.items():
        path = _sample_path(analysis_dir, sample, "consensus", "bcftools", f"{sample}.consensus.fasta")
        records = []
        for seg, n_count in seg_ns.items():
            seq = ("N" * n_count) + ("A" * (100 - n_count))
            records.append(f">{_sequence_id(sample, seg)}\n{seq}\n")
        _write(path, "".join(records))


def _build_annotation(analysis_dir):
    # Cleavage site annotation only produced for the H5N1 sample (S1).
    sample = SAMPLE_STANDARD
    path = _sample_path(analysis_dir, sample, "annotation", f"{sample}.cleavage.tsv")
    header = f"{sample}_segment4_HA|misc_feature|HA|1035..1061 cleavage site"
    content = "Cleavage Sequence\tCleavage Site Sequence Header\n"
    content += f"PRRARRVSLVQERG\t{header}\n"
    _write(path, content)


def _build_genoflu(analysis_dir):
    # No genoflu call for the control sample (S2).
    genotypes = {
        SAMPLE_STANDARD: (
            "B3.6",
            "PB2: A3, PB1: A3, PA: A3, HA: A3, NP: A3, NA: A3, MP: A3, NS: A3",
        ),
        SAMPLE_SALVAGE: (
            "A1.2",
            "PB2: B1, PB1: A2, PA: A2, HA: A2, NP: A2, NA: B1, MP: A2, NS: A2",
        ),
        SAMPLE_EXTERNAL: (
            "B3.13",
            "PB2: C1, PB1: C1, PA: C1, HA: C1, NP: C1, NA: C1, MP: C1, NS: C1",
        ),
    }
    for sample, (genotype, seg_list) in genotypes.items():
        path = _sample_path(analysis_dir, sample, "genoflu", f"{sample}.genoflu.tsv")
        content = "Genotype\tGenotype List Used, >=98.0%\n"
        content += f"{genotype}\t{seg_list}\n"
        _write(path, content)


def _build_nextclade(analysis_dir):
    columns = ["seqName", "clade", "subclade", "legacy-clade", "qc.overallScore", "qc.overallStatus", "alignmentScore"]
    clade_by_sample = {
        SAMPLE_STANDARD: ("2.3.4.4b", "2.3.4.4b.1", "2.3.4.4b", "98.5", "good", "1850.0"),
        SAMPLE_CONTROL: ("2.3.4.4b", "2.3.4.4b.1", "2.3.4.4b", "72.0", "mediocre", "1600.0"),
        SAMPLE_SALVAGE: ("3C.2a1b.2a.2", "3C.2a1b.2a.2a", "3c2.A1b", "99.1", "good", "1920.0"),
        SAMPLE_EXTERNAL: ("2.3.4.4b", "2.3.4.4b.5", "2.3.4.4b", "97.2", "good", "1810.0"),
    }
    for sample in SAMPLE_IDS:
        clade, subclade, legacy, qc_score, qc_status, ha_score = clade_by_sample[sample]
        rows = ["\t".join(columns)]
        for seg in SEGMENTS:
            seq_name = _sequence_id(sample, seg)
            if seg == "HA":
                score = ha_score
            else:
                score = str(float(ha_score) - 400.0)
            rows.append("\t".join([seq_name, clade, subclade, legacy, qc_score, qc_status, score]))
        path = _sample_path(analysis_dir, sample, "nextclade", f"{sample}.nextclade.tsv")
        _write(path, "\n".join(rows) + "\n")

    # "Already published" dataset-version metadata: run-level, not per-sample
    # (no header, 4 columns).
    path = os.path.join(analysis_dir, "aggregate", "nextclade", "nextclade-dataset-versions.csv")
    lines = []
    for sample in SAMPLE_IDS:
        dataset_name, dataset_version = NEXTCLADE_DATASET[sample]
        lines.append(f"{sample},{dataset_name},{dataset_version},{sample}.nextclade.tsv")
    _write(path, "\n".join(lines) + "\n")


def _build_mixtures(analysis_dir):
    # No mixture report for the control sample (S2).
    mixture_rows = {
        SAMPLE_STANDARD: (True, False, True, 0.98, 0.62),
        SAMPLE_SALVAGE: (False, False, False, 0.99, 0.99),
    }
    for sample, (mixture, ha_mix, na_mix, ha_ratio, na_ratio) in mixture_rows.items():
        path = _sample_path(analysis_dir, sample, "mixtures", f"{sample}_mixtures.csv")
        content = "sample_name,mixture_present,ha_mixture_present,na_mixture_present,ha_read_ratio,na_read_ratio\n"
        content += f"{sample},{mixture},{ha_mix},{na_mix},{ha_ratio},{na_ratio}\n"
        _write(path, content)


def _build_software_versions(analysis_dir):
    path = os.path.join(analysis_dir, "pipeline_info", "software_versions.yml")
    content = textwrap.dedent("""\
        GENOFLU:
          genoflu: '1.8.1'
        NEXTCLADE_RUN:
          nextclade: '3.11.0'
        Workflow:
          CFIA-NCFAD/nf-flu: '3.10.0'
    """)
    _write(path, content)
