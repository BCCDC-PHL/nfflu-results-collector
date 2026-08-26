import os
import json
import pandas as pd
from glob import glob
import logging

import nfflu_results_collector.config as config_module


class Nextclade_Results_Collector:
    def __init__(self, config=None):
        """
        Initializes a Nextclade_Results_Collector instance.

        Falls back to the packaged defaults when no config is given, since
        collect_nextclade_results resolves its input paths from config["paths"].
        """
        self.config = config
        if self.config is None:
            self.config = config_module.load_default_config()

    def _ha_only(self):
        return self.config.get("nextclade", {}).get("ha_only", True)

    def _legacy_clade(self):
        return self.config.get("nextclade", {}).get("legacy_clade", False)

    def _filter_nextclade_df(self, df):
        df.insert(0, 'sample', df['seqName'].str.split('_').str[0])
        df = df.dropna(subset=['clade'])

        if self._ha_only():
            df = df.loc[df['seqName'].str.split("_").str[-1] == 'HA']

        if df.empty:
            return df

        df = df.loc[[df['alignmentScore'].idxmax()], :]

        logging.debug(json.dumps({"event_type": "nextclade_alignment_filtered", "alignment_score": float(df['alignmentScore'].values[0])}))

        return df

    def _read_nextclade_datasets(self, tsv_file):
        logging.info(json.dumps({"event_type": "nextclade_datasets_reading", "tsv_file": tsv_file}))
        try:
            df = pd.read_csv(tsv_file, header=None, names=['sample_name', 'nextclade_dataset_name', 'nextclade_dataset_version', 'nextclade_filename'])
            logging.debug(json.dumps({"event_type": "nextclade_datasets_loaded", "dataset_count": len(df)}))
            return df
        except FileNotFoundError:
            logging.warning(json.dumps({"event_type": "nextclade_datasets_file_not_found", "tsv_file": tsv_file}))
            raise
        except Exception as e:
            logging.error(json.dumps({"event_type": "nextclade_datasets_read_error", "tsv_file": tsv_file, "error": str(e)}))
            raise

    def collect_nextclade_results(self, analysis_dir):
        """
        Collect and aggregate nextclade results from multiple samples.

        Per-sample TSVs live in each sample's own output directory
        (paths.nextclade_tsvs), while the dataset name/version metadata
        is a run-level file published in paths.nextclade_dir. This
        collector does not reach into Nextflow execution internals (logs,
        work directories) to find that metadata. Publishing it is the
        pipeline orchestrator's job (e.g. auto-nfflu's
        post_analysis_nf_flu, which has the run context needed to do so).
        If the file is absent, dataset name/version are reported as 'N/A'
        for every sample.

        Args:
            analysis_dir: The nf-flu output directory
        """
        logging.info(json.dumps({"event_type": "nextclade_results_collection_started", "analysis_dir": analysis_dir}))

        if self._ha_only():
            logging.info(json.dumps({"event_type": "nextclade_ha_only_filtering_enabled"}))

        paths = self.config["paths"]
        datasets_csv_path = os.path.join(analysis_dir, paths["nextclade_dir"], 'nextclade-dataset-versions.csv')

        try:
            datasets_df = self._read_nextclade_datasets(datasets_csv_path)
            datasets_dict = datasets_df.set_index('nextclade_filename').to_dict('index')
        except FileNotFoundError:
            datasets_dict = {}
            logging.warning(json.dumps({"event_type": "no_nextclade_datasets_data", "path": datasets_csv_path}))

        # One nextclade directory per sample: <analysis_dir>/<sample>/nextclade
        sample_dirs = [d for d in glob(os.path.join(analysis_dir, paths["nextclade_tsvs"].format(sample='*')))
                       if os.path.isdir(d)]

        logging.info(json.dumps({"event_type": "sample_directories_found", "sample_directory_count": len(sample_dirs)}))

        collect_dfs = []

        for sample_dir in sample_dirs:
            sample_name = os.path.basename(os.path.dirname(sample_dir))
            logging.info(json.dumps({"event_type": "sample_processing_started", "sample_name": sample_name}))

            # Collect all nextclade.tsv files for this sample
            tsv_files = glob(os.path.join(sample_dir, '*.nextclade.tsv'))

            if len(tsv_files) == 0:
                logging.warning(json.dumps({"event_type": "no_nextclade_tsv_files_found", "sample_name": sample_name}))
                continue

            logging.debug(json.dumps({"event_type": "nextclade_tsv_files_found", "sample_name": sample_name, "tsv_file_count": len(tsv_files)}))

            # Process each TSV file
            dataframes = []
            for tsv_file in tsv_files:
                logging.debug(json.dumps({"event_type": "tsv_file_processing_started", "tsv_file": tsv_file}))
                # Read the TSV file
                try:
                    df = pd.read_csv(tsv_file, sep='\t')
                except Exception as e:
                    logging.error(json.dumps({"event_type": "tsv_file_read_failed", "tsv_file": tsv_file, "error": str(e)}))
                    continue

                # Parse dataset name from filename
                filename = os.path.basename(tsv_file)

                dataset_name = datasets_dict[filename]['nextclade_dataset_name'] if datasets_dict and filename in datasets_dict else 'N/A'
                dataset_version = datasets_dict[filename]['nextclade_dataset_version'] if datasets_dict and filename in datasets_dict else 'N/A'

                # Add dataset column
                filtered_df = self._filter_nextclade_df(df)

                insert_position = filtered_df.columns.get_loc('clade')

                filtered_df.insert(insert_position, 'dataset', dataset_name)
                filtered_df.insert(insert_position, 'segment', filtered_df['seqName'].str.split('_').str[-1])

                filtered_df['nextclade_dataset_name'] = dataset_name
                filtered_df['nextclade_dataset_version'] = dataset_version

                dataframes.append(filtered_df)

            if not dataframes:
                logging.warning(json.dumps({"event_type": "no_valid_dataframes_collected", "sample_name": sample_name}))
                continue

            # Concatenate all dataframes for this sample
            sample_df = pd.concat(dataframes, ignore_index=True)

            sample_df = sample_df.drop(columns=['index'], errors='ignore')

            sample_df = sample_df.loc[sample_df.groupby('sample')['alignmentScore'].idxmax()]

            collect_dfs.append(sample_df)

        logging.info(json.dumps({"event_type": "nextclade_results_collection_completed"}))

        # Concatenate all collected dataframes
        if not collect_dfs:
            logging.warning(json.dumps({"event_type": "no_valid_nextclade_dataframes_collected"}))
            return None

        final_df = pd.concat(collect_dfs, ignore_index=True)
        logging.debug(json.dumps({"event_type": "nextclade_results_concatenated", "config": self.config}))

        if self._legacy_clade():
            final_df['clade'] = final_df['legacy-clade']
            logging.info(json.dumps({"event_type": "remapping_to_legacy_clade_column"}))

        return final_df
