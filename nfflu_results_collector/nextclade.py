import os 
import json
import pandas as pd 
import shutil
from glob import glob 
import logging 

import nfflu_results_collector.tools as tools

class Nextclade_Results_Collector:
    def __init__(self, config=None):
        """
        Initializes a Nextclade_Results_Collector instance.
        
        This class provides methods to collect, filter, and publish Nextclade results.
        No parameters are required for initialization.
        """
        self.config = config
        if self.config is None:
            self.config = {}

    def _filter_nextclade_df(self, df):
        df.insert(0, 'sample', df['seqName'].str.split('_').str[0])
        df = df.dropna(subset=['clade'])

        if df.empty:
            return df
        
        df = df.loc[[df['alignmentScore'].idxmax()],:]

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

    def _publish_nextclade_datasets(self, nextclade_results_dir, output_datasets_path):
        log_path = tools.glob_single(os.path.join(os.path.dirname(nextclade_results_dir), '*nextflow.log'))
        if not log_path or not os.path.exists(log_path):
            logging.error(json.dumps({"event_type": "nextflow_log_not_found", "log_path": str(log_path), "search_dir": os.path.dirname(nextclade_results_dir)}))
            raise FileNotFoundError(f"Nextflow log file not found or ambiguous: {log_path}")
        work_dirs = tools.extract_workdirs(log_path)

        agg_process_name = None
        for process in work_dirs:
            if "AGG_NEXTCLADE_TSV" in process:
                agg_process_name = process
                break
        else:
            logging.error(json.dumps({"event_type": "agg_nextclade_tsv_process_not_found", "available_processes": list(work_dirs.keys())}))
            raise ValueError(f'Could not find AGG_NEXTCLADE_TSV process')

        datasets_path = os.path.join(work_dirs[agg_process_name], 'nextclade-tsv-outputs.csv')

        if not os.path.exists(datasets_path):
            logging.warning(json.dumps({"event_type": "nc_datasets_not_found_in_work_dir", "datasets_path": datasets_path}))
            raise FileNotFoundError(f"Nextclade datasets file not found: {datasets_path}")
        
        logging.info(json.dumps({"event_type": "nextclade_datasets_publishing", "output_path": output_datasets_path}))
        shutil.copy2(datasets_path, output_datasets_path)
    

    def collect_nextclade_results(self, nextclade_results_dir):
        """
        Collect and aggregate nextclade results from multiple samples.
        
        Args:
            nextclade_results_dir: Directory containing sample subdirectories with nextclade TSV files
        """
        logging.info(json.dumps({"event_type": "nextclade_results_collection_started", "nextclade_results_dir": nextclade_results_dir}))

        datasets_csv_path = os.path.join(nextclade_results_dir, 'nextclade-dataset-versions.csv')

        try:
            self._publish_nextclade_datasets(nextclade_results_dir, datasets_csv_path)
        except FileNotFoundError as e:
            logging.warning(json.dumps({"event_type": "nextclade_datasets_extraction_failed", "error": str(e)}))
        
        try:
            datasets_df = self._read_nextclade_datasets(datasets_csv_path)
            datasets_dict = datasets_df.set_index('nextclade_filename').to_dict('index')
        except FileNotFoundError as e:
            datasets_dict = {}
            logging.error(json.dumps({"event_type": "no_nextclade_datasets_data", "path": datasets_csv_path}))
            

        # Get all sample subdirectories
        sample_dirs = [d for d in glob(os.path.join(nextclade_results_dir, '*')) 
                       if os.path.isdir(d)]
        
        logging.info(json.dumps({"event_type": "sample_directories_found", "sample_directory_count": len(sample_dirs)}))

        collect_dfs = []
        
        for sample_dir in sample_dirs:
            sample_name = os.path.basename(sample_dir)
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
            
            # Write output CSV file
            # output_file = os.path.join(nextclade_results_dir, f'{sample_name}_nextclade.tsv')
            # try:
            #     sample_df.to_csv(output_file, sep='\t', index=False)
            #     logging.info(json.dumps({"event_type": "nextclade_results_written", "sample_name": sample_name, "output_file": output_file}))
            # except Exception as e:
            #     logging.error(json.dumps({"event_type": "output_file_write_failed", "output_file": output_file, "error": str(e)}))
            #     raise

            collect_dfs.append(sample_df)
        
        logging.info(json.dumps({"event_type": "nextclade_results_collection_completed"}))

        # Concatenate all collected dataframes
        if not collect_dfs:
            logging.warning(json.dumps({"event_type": "no_valid_nextclade_dataframes_collected"}))
            return None

        final_df = pd.concat(collect_dfs, ignore_index=True)
        logging.debug(json.dumps({"event_type": "nextclade_results_concatenated", "config": self.config}))
        
        if self.config.get("legacy-clade", False):
            final_df['clade'] = final_df['legacy-clade']
            logging.info(json.dumps({"event_type": "remapping_to_legacy_clade_column"}))
        
        return final_df