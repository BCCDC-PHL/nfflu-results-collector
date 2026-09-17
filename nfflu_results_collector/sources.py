"""Registry of result sources: where each upstream nf-flu result file
lives, how to parse it, and which output columns it contributes.

Adding a new result type means adding one ResultSource here (and a parser
function beside it) -- nothing else needs to change to have it flow
through collect_run_summary.

Three granularities:
  - "per_run":    parse(analysis_dir, config) -> DataFrame keyed on 'sample'
                  (a left-merge candidate; empty DataFrame means "no data")
  - "per_sample": parse(analysis_dir, sample, config) -> dict of columns
                  (called once per sample by the orchestrator)
  - "global":     parse(analysis_dir, config) -> dict of columns, broadcast
                  identically onto every row (e.g. software versions)
"""
import glob
import json
import logging
import os
import re
from dataclasses import dataclass
from typing import Callable, List

import pandas as pd
import yaml

import nfflu_results_collector.auto as auto
import nfflu_results_collector.nextclade as nextclade
import nfflu_results_collector.schema as schema
import nfflu_results_collector.tools as tools

SEGMENT_PATTERN = re.compile(r'Segment_\d+_([A-Za-z0-9]+)\.')
CLEAVAGE_RANGE_PATTERN = re.compile(r'\|(\d+)\.\.(\d+)')


@dataclass
class ResultSource:
    name: str
    granularity: str  # "per_run" | "per_sample" | "global"
    columns: List[str]
    parse: Callable


def _parse_subtype(analysis_dir, config):
    subtype_file = os.path.join(analysis_dir, config["paths"]["subtype_results"])

    if not os.path.exists(subtype_file):
        return pd.DataFrame()

    try:
        df = pd.read_csv(subtype_file, index_col=0)
        df = df[['sample', 'H_type', 'N_type']].rename(columns={
            'H_type': 'HA_subtype',
            'N_type': 'NA_subtype',
        })
        df['HA_subtype'] = df['HA_subtype'].apply(lambda x: "H" + str(int(x)) if pd.notna(x) else x)
        df['NA_subtype'] = df['NA_subtype'].apply(lambda x: "N" + str(int(x)) if pd.notna(x) else x)
        # fillna only for the readable subtype field; HA_subtype/NA_subtype
        # keep their NaNs since presence/absence drives subtype_HA_NA_status.
        df['subtype'] = df['HA_subtype'].fillna("") + df['NA_subtype'].fillna("")
        return df
    except Exception as e:
        logging.warning(json.dumps({"event_type": "subtype_file_read_error", "subtype_file": subtype_file, "error": str(e)}))
        return pd.DataFrame()


def _parse_idxstats(analysis_dir, sample, config):
    pattern = os.path.join(analysis_dir, config["paths"]["idxstats"].format(sample=sample))
    segments = config["segments"]

    result = {}
    for seg in segments:
        result[f'{seg}_reads_mapped'] = 0
        result[f'{seg}_seq_length'] = 0

    for filepath in glob.glob(pattern):
        try:
            df = pd.read_csv(filepath, sep='\t', header=None,
                              names=['ref', 'seq_length', 'reads_mapped', 'unmapped'])

            if df.shape[0] == 2:
                seq_length = int(df.iloc[0]['seq_length'])
                reads_mapped = int(df.iloc[0]['reads_mapped'])

                filename = os.path.basename(filepath)
                match = SEGMENT_PATTERN.search(filename)
                seg_name = match.group(1) if match else None

                if seg_name is None:
                    logging.error(json.dumps({"event_type": "segment_name_extraction_failed", "filename": filename}))
                    continue

                result[f'{seg_name}_reads_mapped'] = reads_mapped
                result[f'{seg_name}_seq_length'] = seq_length
            else:
                logging.warning(json.dumps({"event_type": "unexpected_idxstats_format", "filepath": filepath}))
                continue
        except Exception as e:
            logging.error(json.dumps({"event_type": "idxstats_file_read_error", "filepath": filepath, "error": str(e)}))
            continue

    return result


def _parse_consensus_completeness(analysis_dir, sample, config):
    fasta_path = os.path.join(analysis_dir, config["paths"]["consensus_fasta"].format(sample=sample))
    segments = config["segments"]

    if not os.path.exists(fasta_path):
        return {f'{seg}_consensus_completeness': 0.0 for seg in segments}

    try:
        df = tools.calculate_completeness(fasta_path, segments)
        return df.iloc[0].to_dict()
    except Exception as e:
        logging.warning(json.dumps({"event_type": "completeness_calculation_error", "sample_name": sample, "error": str(e)}))
        return {f'{seg}_consensus_completeness': 0.0 for seg in segments}


def _parse_cleavage(analysis_dir, sample, config):
    cleavage_file = os.path.join(analysis_dir, config["paths"]["cleavage"].format(sample=sample))

    result = {
        'HPAI_cleave_start': None,
        'HPAI_cleave_end': None,
        'HPAI_cleavage_site_motif': None,
    }

    if not os.path.exists(cleavage_file):
        return result

    try:
        df = pd.read_csv(cleavage_file, sep='\t')
        if not df.empty:
            row = df.iloc[0]
            result['HPAI_cleavage_site_motif'] = row['Cleavage Sequence']

            header = row['Cleavage Site Sequence Header']
            if not pd.isna(header):
                match = CLEAVAGE_RANGE_PATTERN.search(header)
                if match:
                    result['HPAI_cleave_start'] = int(match.group(1))
                    result['HPAI_cleave_end'] = int(match.group(2))
    except pd.errors.EmptyDataError:
        logging.warning(json.dumps({"event_type": "cleavage_file_empty", "sample_name": sample, "cleavage_file": cleavage_file}))
    except Exception as e:
        logging.warning(json.dumps({"event_type": "cleavage_file_read_error", "sample_name": sample, "cleavage_file": cleavage_file, "error": str(e)}))

    return result


def _parse_genotype(analysis_dir, sample, config):
    genoflu_file = os.path.join(analysis_dir, config["paths"]["genoflu"].format(sample=sample))
    segments = config["segments"]

    result = {f"GenoFLU_{key}": None for key in ['Genotype'] + segments}

    if not os.path.exists(genoflu_file):
        return result

    try:
        df = pd.read_csv(genoflu_file, sep='\t')

        result['GenoFLU_Genotype'] = df.loc[0, 'Genotype']
        segment_genotypes = df.loc[0, 'Genotype List Used, >=98.0%']

        if pd.isna(segment_genotypes) or segment_genotypes == '':
            return result

        segment_genotypes = dict([x.strip().split(":") for x in segment_genotypes.split(",")])
        segment_genotypes['M'] = segment_genotypes.pop('MP', None)

        result.update({f"GenoFLU_{seg}": genotype for seg, genotype in segment_genotypes.items()})

    except Exception as e:
        logging.warning(json.dumps({"event_type": "genoflu_file_read_error", "sample_name": sample, "genoflu_file": genoflu_file, "error": str(e)}))

    return result


def _parse_nextclade(analysis_dir, config):
    nextclade_rc = nextclade.Nextclade_Results_Collector(config)

    try:
        df = nextclade_rc.collect_nextclade_results(analysis_dir)

        if df is None or df.empty:
            return pd.DataFrame()

        rename_map = {
            'clade': 'Nextclade_clade',
            'subclade': 'Nextclade_subclade',
            'legacy-clade': 'Nextclade_legacy_clade',
            'qc.overallScore': 'Nextclade_qc.overallScore',
            'qc.overallStatus': 'Nextclade_qc.overallStatus',
            'nextclade_dataset_name': 'Nextclade_dataset_name',
            'nextclade_dataset_version': 'Nextclade_dataset_version',
        }

        cols_to_keep = ['sample'] + [col for col in rename_map if col in df.columns]
        df = df[cols_to_keep].rename(columns=rename_map)
        return df
    except Exception as e:
        logging.warning(json.dumps({"event_type": "nextclade_file_read_error", "error": str(e)}))
        return pd.DataFrame()


def _parse_provenance(analysis_dir, config):
    software_versions_path = os.path.join(analysis_dir, config["paths"]["software_versions"])

    result = {
        'genoflu_version': None,
        'nextclade_version': None,
        'nfflu_version': None,
    }

    if not os.path.exists(software_versions_path):
        return result

    try:
        with open(software_versions_path, 'r') as f:
            data = yaml.safe_load(f)

        if 'GENOFLU' in data and 'genoflu' in data['GENOFLU']:
            result['genoflu_version'] = data['GENOFLU']['genoflu']

        if 'NEXTCLADE_RUN' in data and 'nextclade' in data['NEXTCLADE_RUN']:
            result['nextclade_version'] = data['NEXTCLADE_RUN']['nextclade']

        if 'Workflow' in data and 'CFIA-NCFAD/nf-flu' in data['Workflow']:
            result['nfflu_version'] = data['Workflow']['CFIA-NCFAD/nf-flu']

    except Exception as e:
        logging.warning(json.dumps({"event_type": "software_versions_file_read_error", "software_versions_path": software_versions_path, "error": str(e)}))

    return result


SUBTYPE_SOURCE = ResultSource(
    name="subtype", granularity="per_run",
    columns=schema.SUBTYPE_COLUMNS, parse=_parse_subtype,
)
NEXTCLADE_SOURCE = ResultSource(
    name="nextclade", granularity="per_run",
    columns=schema.NEXTCLADE_COLUMNS, parse=_parse_nextclade,
)
# status_* columns are dynamic (one per upstream pipeline) and only ever
# populated when auto-nfflu mode is enabled; see collector.py.
PIPELINE_STATUS_SOURCE = ResultSource(
    name="pipeline_status", granularity="per_run",
    columns=[], parse=auto.parse_pipeline_status,
)
IDXSTATS_SOURCE = ResultSource(
    name="idxstats", granularity="per_sample",
    columns=(
        [f"{seg}_reads_mapped" for seg in schema.DISPLAY_SEGMENT_ORDER]
        + [f"{seg}_seq_length" for seg in schema.DISPLAY_SEGMENT_ORDER]
    ),
    parse=_parse_idxstats,
)
COMPLETENESS_SOURCE = ResultSource(
    name="consensus_completeness", granularity="per_sample",
    columns=[f"{seg}_consensus_completeness" for seg in schema.DISPLAY_SEGMENT_ORDER],
    parse=_parse_consensus_completeness,
)
CLEAVAGE_SOURCE = ResultSource(
    name="cleavage", granularity="per_sample",
    columns=schema.HPAI_COLUMNS, parse=_parse_cleavage,
)
GENOTYPE_SOURCE = ResultSource(
    name="genotype", granularity="per_sample",
    columns=schema.GENOFLU_COLUMNS, parse=_parse_genotype,
)
PROVENANCE_SOURCE = ResultSource(
    name="provenance", granularity="global",
    columns=schema.PROVENANCE_COLUMNS, parse=_parse_provenance,
)

PER_RUN_SOURCES = [SUBTYPE_SOURCE, NEXTCLADE_SOURCE]
PER_SAMPLE_SOURCES = [IDXSTATS_SOURCE, COMPLETENESS_SOURCE, CLEAVAGE_SOURCE, GENOTYPE_SOURCE]
GLOBAL_SOURCES = [PROVENANCE_SOURCE]

# Declared columns, across every registered source, that are part of the
# frozen run_summary.csv schema (excludes pipeline_status, whose columns
# are dynamic and only appended in auto-nfflu mode).
ALL_SOURCE_COLUMNS = [
    col
    for src in (PER_RUN_SOURCES + PER_SAMPLE_SOURCES + GLOBAL_SOURCES)
    for col in src.columns
]


def compute_tree_pass(per_sample_df, config):
    """tree_pass[seg] = 1 if {seg}_consensus_completeness > threshold else 0."""
    threshold = config["tree_pass"]["threshold"]
    for seg in config["segments"]:
        completeness_col = f'{seg}_consensus_completeness'
        tree_pass_col = f'{seg}_tree_pass'
        if completeness_col in per_sample_df.columns:
            per_sample_df[tree_pass_col] = (per_sample_df[completeness_col] > threshold).astype(int)
        else:
            per_sample_df[tree_pass_col] = 0
    return per_sample_df


def compute_subtype_status(df):
    """subtype_HA_NA_status: a description of HA/NA subtype call presence,
    computed from HA_subtype/NA_subtype regardless of what the subtype
    source itself reported (a prior "Genotype" passthrough value here was
    always overwritten by this computation anyway)."""
    description_map = {
        0: "No subtype for HA or NA",
        1: "HA subtype successful",
        2: "NA subtype successful",
        3: "HA and NA subtype successful",
    }
    if 'HA_subtype' not in df.columns or 'NA_subtype' not in df.columns:
        return df
    status_code = df['HA_subtype'].notna().astype(int) + df['NA_subtype'].notna().astype(int) * 2
    df['subtype_HA_NA_status'] = status_code.map(description_map)
    return df
