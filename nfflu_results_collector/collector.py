import argparse
import glob
import json
import os
import re
import sys
import pandas as pd
import logging
import yaml

import nfflu_results_collector.tools as tools
import nfflu_results_collector.config as config
import nfflu_results_collector.auto as auto
import nfflu_results_collector.nextclade as nextclade

# Compile regex patterns once for performance
SEGMENT_PATTERN = re.compile(r'Segment_\d+_([A-Za-z0-9]+)\.')
CLEAVAGE_RANGE_PATTERN = re.compile(r'\|(\d+)\.\.(\d+)')
SEGMENTS = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']

class Nfflu_Results_Collector:
    def __init__(self, user_config=None):
        if user_config is None:
            user_config = {}
        self.config = config.load_default_config()
        self.config.update(user_config)

    def collect_run_summary(self, analysis_dir, output_summary_file):
        """Main method to collect all results and merge them into a single CSV."""
        
        output_summary_file = output_summary_file.rstrip(os.sep)
        
        # Step 1: Collect sample names and parse them
        if self.config.get('auto-nfflu', False):
            sample_names = auto.collect_auto_nfflu_names(analysis_dir)
            logging.info(json.dumps({"event_type": "samples_collected", "sample_count": len(sample_names), "source": "auto_nfflu_start_samplesheet"}))
        else:
            sample_names = tools.collect_nfflu_fastq_names(analysis_dir)
            logging.info(json.dumps({"event_type": "samples_collected", "sample_count": len(sample_names), "source": "fastq_directory"}))

        if len(sample_names) == 0:
            logging.warning(json.dumps({"event_type": "no_samples_found", "fastq_directory": os.path.join(analysis_dir, "fastq")}))
            return

        samples_df = pd.DataFrame(sample_names, columns=['sample'])

        # Parse sample name into components using pandas str.split
        # Format is typically: CID-Plate-Index-Well (where CID can contain dashes)

        def extract_fields(sample_name):

            standard_format_search = re.match(r"^[A-Z][0-9]{10}-[0-9]{4,}-[A-Z0-9]{1,}-[A-Z][0-9]{2}$", sample_name)

            if standard_format_search:
                fields = sample_name.split("-")
                return fields
            
            control_format_search = re.match(r"^([A-Za-z]{3,}[0-9]{8})-[A-Za-z]+-([0-9]{4,})-([A-Z0-9]{1,})-([A-Z][0-9]{2})$", sample_name)

            if control_format_search:
                return [control_format_search.group(1), control_format_search.group(2), control_format_search.group(3), control_format_search.group(4)]

            salvage_fields_search = re.search(r"-[0-9]{4,9}-[A-Z0-9]{1,}-[A-Z][0-9]{2}$", sample_name)
            
            if salvage_fields_search:
                plate, index, well = salvage_fields_search.group(0).lstrip("-").split("-")
            else:
                plate, index, well = pd.NA, pd.NA, pd.NA
            
            salvage_sample_cid_search = re.search(r"^([A-Z][0-9]{10}).*", sample_name)
            salvage_control_cid_search = re.search(r"^([A-Za-z]{3,}[0-9]{8}).*", sample_name)

            if salvage_sample_cid_search:
                cid = salvage_sample_cid_search.group(1)
            elif salvage_control_cid_search:
                cid = salvage_control_cid_search.group(1)
            else:
                cid = sample_name

            return [cid, plate, index, well]

        samples_df[['CID', 'Plate', 'Index', 'Well']] = samples_df['sample'].apply(extract_fields).tolist()
        samples_df['FastQID'] = samples_df['sample']

        # Add Run column
        run_parts = analysis_dir.rstrip(os.sep).split(os.sep)
        if len(run_parts) >= 3:
            samples_df['Run'] = run_parts[-3]
        else:
            logging.warning(json.dumps({"event_type": "run_name_extraction_failed", "analysis_dir": analysis_dir}))
            samples_df['Run'] = pd.NA
        
        # Step 2: Collect subtype results
        subtype_df = self._collect_subtype(analysis_dir)
        
        # Step 3: Collect nextclade results
        nextclade_df = self._collect_nextclade(analysis_dir)
        
        # Step 4: For each sample, collect per-sample data
        per_sample_data = []
        
        for _, row in samples_df.iterrows():
            sample_name = row['sample']
            sample_data = {'sample': sample_name}
            
            # Collect all data using unpacking for cleaner code
            sample_data.update(
                **self._collect_idxstats(analysis_dir, sample_name),
                **self._collect_consensus_completeness(analysis_dir, sample_name),
                **self._collect_cleavage(analysis_dir, sample_name),
                **self._collect_genotype(analysis_dir, sample_name),
            )
            
            per_sample_data.append(sample_data)
        
        per_sample_df = pd.DataFrame(per_sample_data)
        
        # Calculate tree_pass columns using vectorized operations
        for seg in SEGMENTS:
            completeness_col = f'{seg}_consensus_completeness'
            tree_pass_col = f'{seg}_tree_pass'
            if completeness_col in per_sample_df.columns:
                per_sample_df[tree_pass_col] = (per_sample_df[completeness_col] > 90).astype(int)
            else:
                per_sample_df[tree_pass_col] = 0
        
        # Step 5: Merge all dataframes
        output_df = samples_df.copy()
        
        description_map = {
            0: "No subtype for HA or NA",
            1: "HA subtype successful",
            2: "NA subtype successful",
            3: "HA and NA subtype successful"
        }
        
        # Merge subtype data
        if not subtype_df.empty:
            output_df = output_df.merge(subtype_df, on='sample', how='left')
            output_df['subtype_HA_NA_status'] = output_df['HA_subtype'].notna().astype(int) + output_df['NA_subtype'].notna().astype(int) * 2
            output_df['subtype_HA_NA_status'] = output_df['subtype_HA_NA_status'].map(description_map)
        
        # Merge nextclade data
        if not nextclade_df.empty:
            output_df = output_df.merge(nextclade_df, on='sample', how='left')
        
        # Merge per-sample data
        output_df = output_df.merge(per_sample_df, on='sample', how='left')
        
        # Step 6: Reorder columns to match expected output
        # FastQID, CID, Plate, Index, Well, Run, subtype_HA_NA_status, HA_subtype, ...

        output_df = output_df.assign(**self._collect_provenance(analysis_dir))

        if "sample" in output_df.columns:
            output_df = output_df.drop(columns=['sample'])

        expected_columns = self.config.get('expected_columns', [])

        if expected_columns:
            # Ensure all columns exist using pandas assignment
            missing_cols = [col for col in expected_columns if col not in output_df.columns]
            for col in missing_cols:
                logging.warning(json.dumps({"event_type": "expected_column_missing", "column": col}))
                output_df[col] = None

            unknown_cols = [col for col in output_df.columns if col not in expected_columns]
            for col in unknown_cols:
                logging.warning(json.dumps({"event_type": "unexpected_column_found", "column": col}))

            # Reorder columns
            output_df = output_df[expected_columns]

        else:
            logging.warning(json.dumps({"event_type": "no_expected_columns_defined"}))

        if self.config.get('auto-nfflu', False):
            pipeline_status_df = auto.compute_pipeline_status_columns(analysis_dir)

            if not pipeline_status_df.empty:
                output_df = output_df.merge(pipeline_status_df, left_on='FastQID', right_on='ID', how='left')
                output_df = output_df.drop(columns=['ID'])

        # Step 7: Write to CSV
        output_summary_dir = os.path.dirname(output_summary_file)
        
        if output_summary_dir != '' and not os.path.exists(output_summary_dir):
            os.makedirs(output_summary_dir, exist_ok=True)

        output_df.to_csv(output_summary_file, index=False)
        logging.info(json.dumps({"event_type": "results_written", "output_file": output_summary_file}))


    def _collect_subtype(self, analysis_dir):
        """Collect subtype data from subtype_results.csv."""
        subtype_file = os.path.join(analysis_dir, "subtyping_report", "subtype_results.csv")
        
        if not os.path.exists(subtype_file):
            return pd.DataFrame()

        
        try:
            # Read the CSV file with pandas, handling the first unnamed index column
            df = pd.read_csv(subtype_file, index_col=0)
            
            # Select and rename the columns we need
            df = df[['sample', 'Genotype', 'H_type', 'N_type']].rename(columns={
                'Genotype': 'subtype_HA_NA_status',
                'H_type': 'HA_subtype',
                'N_type': 'NA_subtype'
            })
            df['HA_subtype'] = df['HA_subtype'].apply(lambda x : "H" + str(int(x)) if pd.notna(x) else x)
            df['NA_subtype'] = df['NA_subtype'].apply(lambda x : "N" + str(int(x)) if pd.notna(x) else x)
            df['subtype'] = df['HA_subtype'].fillna("") + df['NA_subtype'].fillna("")     # choose to fillna only for readable subtype field, need nans later 
            return df
            
        except Exception as e:
            logging.warning(json.dumps({"event_type": "subtype_file_read_error", "subtype_file": subtype_file, "error": str(e)}))
            return pd.DataFrame()

    def _collect_idxstats(self, analysis_dir, sample_name):
        """Collect idxstats data for a single sample across all segments."""
        # Pattern: ${analysis_run_dir}/mapping/${sample_name}/${sample_name}*idxstats
        pattern = os.path.join(analysis_dir, "mapping", sample_name, f"{sample_name}*.idxstats")
        
        result = {}
        
        # Initialize all segments to 0
        for seg in SEGMENTS:
            result[f'{seg}_reads_mapped'] = 0
            result[f'{seg}_seq_length'] = 0
        
        for filepath in glob.glob(pattern):
            try:
                # Read idxstats file using pandas (tab-separated, no header)
                df = pd.read_csv(filepath, sep='\t', header=None, 
                                names=['ref', 'seq_length', 'reads_mapped', 'unmapped'])
                
                if df.shape[0] == 2:
                    # Get the first row data
                    seq_length = int(df.iloc[0]['seq_length'])
                    reads_mapped = int(df.iloc[0]['reads_mapped'])
                    
                    # Extract segment number from filename using regex
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

    def _collect_consensus_completeness(self, analysis_dir, sample_name):
        """Calculate consensus completeness for a single sample."""
        # Path: ${analysis_run_dir}/consensus/bcftools/${sample_name}.consensus.fasta
        fasta_path = os.path.join(analysis_dir, "consensus", "bcftools", f"{sample_name}.consensus.fasta")
        
        if not os.path.exists(fasta_path):
            # Return zeros for all segments
            result = {f'{seg}_consensus_completeness': 0.0 for seg in SEGMENTS}
            return result
        
        try:
            df = tools.calculate_completeness(fasta_path)
            # Convert to dict with first row values
            result = df.iloc[0].to_dict()
            return result
        except Exception as e:
            logging.warning(json.dumps({"event_type": "completeness_calculation_error", "sample_name": sample_name, "error": str(e)}))
            return {f'{seg}_consensus_completeness': 0.0 for seg in SEGMENTS}

    def _collect_cleavage(self, analysis_dir, sample_name):
        """Collect HPAI cleavage data for a single sample."""
        # Pattern: ${analysis_run_dir}/annotation/${sample_name}/${sample_name}.cleavage.tsv
        cleavage_file = os.path.join(analysis_dir, "annotation", sample_name, f"{sample_name}.cleavage.tsv")
        
        result = {
            'HPAI_cleave_start': None,
            'HPAI_cleave_end': None,
            'HPAI_cleavage_site_motif': None
        }
        
        if not os.path.exists(cleavage_file):
            return result
        
        try:
            df = pd.read_csv(cleavage_file, sep='\t')
            if not df.empty:
                row = df.iloc[0]
                
                # Extract cleavage site motif
                result['HPAI_cleavage_site_motif'] = row['Cleavage Sequence']
                
                # Parse start and end from "Cleavage Site Sequence Header"
                # Example: R3940252709-3321-12-D09_segment4_HA|misc_feature|HA|1035..1061 cleavage site
                header = row['Cleavage Site Sequence Header']
                if not pd.isna(header):
                    # Extract the range like "1035..1061"
                    match = CLEAVAGE_RANGE_PATTERN.search(header)
                    if match:
                        result['HPAI_cleave_start'] = int(match.group(1))
                        result['HPAI_cleave_end'] = int(match.group(2))

        except pd.errors.EmptyDataError as e:
            logging.warning(json.dumps({"event_type": "cleavage_file_empty", "sample_name": sample_name, "cleavage_file": cleavage_file}))
        except Exception as e:
            logging.warning(json.dumps({"event_type": "cleavage_file_read_error", "sample_name": sample_name, "cleavage_file": cleavage_file, "error": str(e)}))
        
        return result

    def _collect_nextclade(self, analysis_dir):
        """Collect nextclade results from nextclade.tsv."""
        nextclade_path = os.path.join(analysis_dir, "nextclade")
        nextclade_rc = nextclade.Nextclade_Results_Collector(self.config)
        logging.debug(json.dumps({"event_type": "collecting_nextclade_results", "config": self.config}))

        try:
            df = nextclade_rc.collect_nextclade_results(nextclade_path)

            if df is None or df.empty:
                return pd.DataFrame()
            
            # final run_summary is only interested in HA clade call 
            #df = df.loc[df['seqName'].str.split("_").str[-1] == 'HA']
            
            # Rename columns to match target
            rename_map = {
                'clade': 'Nextclade_clade',
                'subclade': 'Nextclade_subclade',
                'legacy-clade': 'Nextclade_legacy_clade',
                'qc.overallScore': 'Nextclade_qc.overallScore',
                'qc.overallStatus': 'Nextclade_qc.overallStatus',
                'nextclade_dataset_name': 'Nextclade_dataset_name',
                'nextclade_dataset_version': 'Nextclade_dataset_version'
            }
            
            # Select only columns that exist and rename them
            cols_to_keep = ['sample'] + [col for col in rename_map.keys() if col in df.columns]
            df = df[cols_to_keep].rename(columns=rename_map)
            
            return df
        except Exception as e:
            logging.warning(json.dumps({"event_type": "nextclade_file_read_error", "error": str(e)}))
            return pd.DataFrame()

    def _collect_genotype(self, analysis_dir, sample_name):
        """Collect genotype data for a single sample from genoflu file."""
        # Pattern: ${analysis_run_dir}/genoflu/${sample_name}.genoflu.tsv
        genoflu_file = os.path.join(analysis_dir, "genoflu", f"{sample_name}.genoflu.tsv")

        result = {f"GenoFLU_{key}" : None for key in ['Genotype'] + SEGMENTS}

        if not os.path.exists(genoflu_file):
            return result
        
        try:
            df = pd.read_csv(genoflu_file, sep='\t')
            
            # Get the overall genotype - column may have leading/trailing spaces
            result['GenoFLU_Genotype'] = df.loc[0,'Genotype']
            segment_genotypes = df.loc[0,'Genotype List Used, >=98.0%']

            if pd.isna(segment_genotypes) or segment_genotypes == '':
                return result
        
            segment_genotypes = dict([x.strip().split(":") for x in segment_genotypes.split(",")])
            segment_genotypes['M'] = segment_genotypes.pop('MP', None)

            result.update({f"GenoFLU_{seg}": genotype for seg, genotype in segment_genotypes.items()})

        except Exception as e:
            logging.warning(json.dumps({"event_type": "genoflu_file_read_error", "sample_name": sample_name, "genoflu_file": genoflu_file, "error": str(e)}))
        
        return result

    def _collect_provenance(self, analysis_dir):
        """Collect software version information from software_versions.yml."""
        software_versions_path = os.path.join(analysis_dir, "pipeline_info", "software_versions.yml")
        
        result = {
            'genoflu_version': None,
            'nextclade_version': None,
            'nfflu_version': None
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

    def collect_nextclade_results(self, analysis_dir, nextclade_output_path):
        nextclade_path = os.path.join(analysis_dir, "nextclade")
        nextclade_rc = nextclade.Nextclade_Results_Collector(self.config)
        
        logging.debug(json.dumps({"event_type": "collecting_nextclade_results", "config": self.config}))

        nextclade_df = nextclade_rc.collect_nextclade_results(nextclade_path)

        nextclade_df.to_csv(nextclade_output_path, sep='\t', index=False)
        logging.info(json.dumps({"event_type": "nextclade_results_written", "output_path": nextclade_output_path}))

    def collect_mixture_report(self, analysis_dir, output_mixture_file):
        """Collect mixture report data for all samples in the analysis directory and concatenate them into a single output file."""
        # Pattern: ${analysis_run_dir}/mixture_analysis/${sample_name}_mixture_report.tsv
        mixture_files = glob.glob(os.path.join(analysis_dir, "mixtures", "*", "*_mixtures.csv"))
        
        dfs = []

        for mixture_path in mixture_files:
            sample_name = os.path.basename(mixture_path).replace('_mixtures.csv', '')
            
            try:
                df = pd.read_csv(mixture_path)

                # columns = ['FastQID', 'mixture_present', 'ha_mixture_present', 'na_mixture_present', 'ha_read_ratio', 'na_read_ratio']

                if df.shape[0] != 1:
                    logging.warning(json.dumps({"event_type": "unexpected_mixture_report_rows", "sample_name": sample_name, "row_count": df.shape[0]}))
                    continue

                dfs.append(df)
                
            except Exception as e:
                logging.warning(json.dumps({"event_type": "mixture_report_read_error", "sample_name": sample_name, "mixture_path": mixture_path, "error": str(e)}))
                continue
        
        if not dfs:
            logging.warning(json.dumps({"event_type": "no_valid_mixture_reports_found"}))
            return 

        # Concatenate all dataframes and return the result
        final_df = pd.concat(dfs, ignore_index=True).rename(columns={'sample_name':'FastQID'})
        final_df.to_csv(output_mixture_file, index=False)
        logging.info(json.dumps({"event_type": "mixture_report_written", "output_file": output_mixture_file}))
        
        return 

    def symlink_consensus_fastas(self, analysis_dir, output_dir):
        """Symlink consensus FASTA files to the output directory."""
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Glob pattern to find all consensus FASTA files
        pattern = os.path.join(analysis_dir, "consensus", "bcftools", "*.consensus.fasta")
        
        for fasta_file in glob.glob(pattern):
            basename = os.path.basename(fasta_file)
            dest = os.path.join(output_dir, basename)
            
            # Remove existing symlink/file if it exists
            if os.path.exists(dest) or os.path.islink(dest):
                os.remove(dest)
            
            # Create symlink
            os.symlink(os.path.abspath(fasta_file), dest)
