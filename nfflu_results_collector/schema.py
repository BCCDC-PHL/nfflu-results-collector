"""The run_summary.csv output schema: column set and order.

This is a frozen contract consumed by downstream ingestion, so
CANONICAL_COLUMNS is a literal, structural transcription of the column
order the tool has always produced -- not something clever derived from
segment/metric loops, because the real order doesn't follow either the
canonical genome-segment order (PB2 first, used everywhere else in this
package for processing) or alphabetical order. It's its own thing, so it
gets its own named constant here instead of being buried in a hand-edited
JSON array with no enforcement.
"""
import json
import logging

# Segment display order for output columns only (HA_*, GenoFLU_*, ...).
# Deliberately NOT the same as the canonical/processing segment order
# (config's "segments" list, PB2 first) -- this is just the order the
# report has always used. Don't "fix" this into alphabetical or genome
# order; downstream ingestion depends on the exact column sequence.
DISPLAY_SEGMENT_ORDER = ["HA", "NA", "M", "NP", "NS", "PA", "PB1", "PB2"]

SAMPLE_IDENTITY_COLUMNS = ["FastQID", "CID", "Plate", "Index", "Well", "Run"]

SUBTYPE_COMPUTED_COLUMNS = ["subtype_HA_NA_status"]
SUBTYPE_COLUMNS = ["subtype", "HA_subtype", "NA_subtype"]

SEGMENT_METRIC_COLUMNS = [
    f"{seg}_{metric}"
    for seg in DISPLAY_SEGMENT_ORDER
    for metric in ("reads_mapped", "seq_length", "consensus_completeness", "tree_pass")
]

HPAI_COLUMNS = ["HPAI_cleave_start", "HPAI_cleave_end", "HPAI_cleavage_site_motif"]

NEXTCLADE_COLUMNS = [
    "Nextclade_clade", "Nextclade_subclade", "Nextclade_legacy_clade",
    "Nextclade_qc.overallScore", "Nextclade_qc.overallStatus",
    "Nextclade_dataset_name", "Nextclade_dataset_version",
]

GENOFLU_COLUMNS = ["GenoFLU_Genotype"] + [f"GenoFLU_{seg}" for seg in DISPLAY_SEGMENT_ORDER]

PROVENANCE_COLUMNS = ["nextclade_version", "genoflu_version", "nfflu_version"]

CANONICAL_COLUMNS = (
    SAMPLE_IDENTITY_COLUMNS
    + SUBTYPE_COMPUTED_COLUMNS
    + SUBTYPE_COLUMNS
    + SEGMENT_METRIC_COLUMNS
    + HPAI_COLUMNS
    + NEXTCLADE_COLUMNS
    + GENOFLU_COLUMNS
    + PROVENANCE_COLUMNS
)


# Mirrors pileup_tools.mixture_detector.MixtureDetector.COLUMN_NAMES, which
# writes the per-sample files the mixture report concatenates.
MIXTURE_COLUMNS = [
    'sample_name', 'subtype', 'mixture_present', 'ha_mixture_present', 'na_mixture_present',
    'primary_ha_subtype', 'primary_ha_reads', 'secondary_ha_subtype', 'secondary_ha_reads', 'ha_read_ratio',
    'primary_na_subtype', 'primary_na_reads', 'secondary_na_subtype', 'secondary_na_reads', 'na_read_ratio',
    'initial_reads', 'pass_qc_reads', 'fail_qc_reads', 'match_reads', 'nomatch_reads',
]

# Every pipeline that can contribute a status column. auto-nfflu's
# pipeline_status.csv carries only the pipelines in one run's chain, so nanopore
# and Illumina runs would otherwise disagree on the header.
STATUS_COLUMNS = [
    'status_cutadapt-nf', 'status_basic-sequence-qc', 'status_downsample-reads',
    'status_basic-nanopore-qc', 'status_nf-flu',
]


def pad_status_columns(df):
    """Add any missing status_* column as blank and move the whole set to the
    end, so every run publishes the same header."""
    for col in STATUS_COLUMNS:
        if col not in df.columns:
            df[col] = None
    return df[[c for c in df.columns if c not in STATUS_COLUMNS] + STATUS_COLUMNS]


def order_and_validate(df, schema=None):
    """Reindex `df` to `schema`'s column order (default CANONICAL_COLUMNS).

    Columns declared in the schema but absent from `df` become all-null.
    Columns present in `df` but NOT declared in the schema are never
    silently dropped: they're appended after the declared schema (stable
    order) and a warning is logged. auto-nfflu's status_* columns arrive
    that way.
    """
    if schema is None:
        schema = CANONICAL_COLUMNS

    missing = [col for col in schema if col not in df.columns]
    for col in missing:
        logging.warning(json.dumps({"event_type": "expected_column_missing", "column": col}))
        df[col] = None

    extra = [col for col in df.columns if col not in schema]
    if extra:
        logging.warning(json.dumps({"event_type": "unexpected_columns_found", "columns": extra}))

    return df[schema + extra]
