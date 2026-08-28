import os
import json
import logging
import pandas as pd
from nfflu_results_collector.tools import glob_single


def collect_auto_nfflu_names(nfflu_results_dir):
    """Sample discovery for auto-nfflu mode: read IDs from the
    start_samplesheet.csv that auto-nfflu writes alongside the nf-flu
    output directory. Only used as a fallback when RunContext.sample_ids
    isn't provided directly by the caller."""
    analysis_output_dir = os.path.dirname(nfflu_results_dir.rstrip(os.sep))
    start_df_path = glob_single(os.path.join(analysis_output_dir, 'samplesheets', '*start_samplesheet.csv'))

    if not start_df_path:
        logging.warning(json.dumps({
            "event_type": "start_samplesheet_not_found",
            "samplesheet_dir": os.path.join(analysis_output_dir, 'samplesheets')
        }))
        return []

    return pd.read_csv(start_df_path)['ID'].tolist()


def parse_pipeline_status(analysis_dir, config):
    """Merge in a pipeline_status.csv if present.

    Per-pipeline-stage completion status used to be re-derived here by
    globbing nf-flu's internal intermediate directories (irma, blast,
    variants, ...) and encoding stage order as magic status codes. That
    knowledge belongs to whatever orchestrates the pipeline chain (e.g.
    auto-nfflu, which already tracks completion per stage) -- this
    collector now just merges in a status file if the orchestrator wrote
    one, keyed by sample ID. Returns an empty DataFrame if absent; the
    collector works standalone without it.
    """
    status_path = os.path.join(analysis_dir, config["paths"]["pipeline_status"])

    if not os.path.exists(status_path):
        return pd.DataFrame()

    try:
        df = pd.read_csv(status_path)
    except Exception as e:
        logging.warning(json.dumps({"event_type": "pipeline_status_read_error", "status_path": status_path, "error": str(e)}))
        return pd.DataFrame()

    if 'ID' not in df.columns:
        logging.warning(json.dumps({"event_type": "pipeline_status_missing_id_column", "status_path": status_path}))
        return pd.DataFrame()

    return df.rename(columns={'ID': 'sample'})
