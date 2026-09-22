"""The run_summary.csv output schema: column set and order.

Downstream ingestion reads this column set and order, so changing it is a
reviewed change pinned by a golden test, not an incidental one.
CANONICAL_COLUMNS is a literal, structural transcription of the column
order the tool produces -- not something clever derived from
segment/metric loops, because the real order doesn't follow either the
canonical genome-segment order (PB2 first, used everywhere else in this
package for processing) or alphabetical order. It's its own thing, so it
gets its own named constant here instead of being buried in a hand-edited
JSON array with no enforcement.
"""
import json
import logging

import pandas as pd

# Segment display order for output columns only (HA_*, GenoFLU_*, ...).
# Deliberately NOT the same as the canonical/processing segment order
# (config's "segments" list, PB2 first) -- this is just the order the
# report has always used. Don't "fix" this into alphabetical or genome
# order; downstream ingestion depends on the exact column sequence.
DISPLAY_SEGMENT_ORDER = ["HA", "NA", "M", "NP", "NS", "PA", "PB1", "PB2"]

SAMPLE_IDENTITY_COLUMNS = ["FastQID", "CID", "Plate", "Index", "Well", "Run"]

# Named as nf-flu and auto-nfflu name them: illumina/nanopore, short/long.
RUN_CONTEXT_COLUMNS = ["platform", "analysis_type"]

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
    + RUN_CONTEXT_COLUMNS
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

def _num(s):
    """Non-null values as numbers, with anything unparseable as NaN, so a
    check never raises on a column that has gone to strings."""
    return pd.to_numeric(s.dropna(), errors="coerce")


# Value sanity checks, keyed by column suffix so one entry covers all eight
# segments. These warn and never modify the frame. An all-null column is
# vacuously fine -- order_and_validate already warns separately about those.
CHECKS = {
    "_consensus_completeness": lambda s: _num(s).between(0, 100).all(),
    "_reads_mapped": lambda s: _num(s).notna().all(),
    "_tree_pass": lambda s: _num(s).isin([0, 1]).all(),
}


def order_and_validate(df, schema=None):
    """Reindex `df` to `schema`'s column order (default CANONICAL_COLUMNS).

    Columns declared in the schema but absent from `df` become all-null.
    Columns present in `df` but NOT declared in the schema are never
    silently dropped: they're appended after the declared schema (stable
    order) and a warning is logged. auto-nfflu's status_* columns arrive
    that way. Values are then checked against CHECKS, which only warns.
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

    out = df[schema + extra]

    for suffix, ok in CHECKS.items():
        for col in (c for c in out.columns if c.endswith(suffix)):
            if not ok(out[col]):
                logging.warning(json.dumps({"event_type": "column_failed_check", "column": col}))

    return out
