import glob
import json
import os
import logging

import pandas as pd

import nfflu_results_collector.auto as auto
import nfflu_results_collector.config as config
import nfflu_results_collector.nextclade as nextclade
import nfflu_results_collector.sample_id as sample_id
import nfflu_results_collector.schema as schema
import nfflu_results_collector.sources as sources
import nfflu_results_collector.tools as tools


class Nfflu_Results_Collector:
    def __init__(self, user_config=None, config_path=None):
        self.config = config.load_config(user_config_path=config_path, overrides=user_config)

    def _resolve_sample_ids(self, analysis_dir):
        if self.config.get('auto-nfflu', False):
            sample_ids = auto.collect_auto_nfflu_names(analysis_dir)
            logging.info(json.dumps({"event_type": "samples_collected", "sample_count": len(sample_ids), "source": "auto_nfflu_start_samplesheet"}))
        else:
            sample_ids = tools.collect_nfflu_fastq_names(analysis_dir, self.config["paths"]["fastq_dir"])
            logging.info(json.dumps({"event_type": "samples_collected", "sample_count": len(sample_ids), "source": "fastq_directory"}))
        return sample_ids

    def _resolve_run_id(self, analysis_dir):
        run_parts = analysis_dir.rstrip(os.sep).split(os.sep)
        if len(run_parts) >= 3:
            return run_parts[-3]
        logging.warning(json.dumps({"event_type": "run_name_extraction_failed", "analysis_dir": analysis_dir}))
        return pd.NA

    def _build_sample_identity(self, sample_ids, run_id):
        records = []
        for sample in sample_ids:
            record = {"sample": sample, "FastQID": sample, "Run": run_id}
            record.update(sample_id.parse_sample_id(sample))
            records.append(record)
        return pd.DataFrame(records)

    def _collect_per_sample(self, analysis_dir, sample_ids):
        records = []
        for sample in sample_ids:
            record = {"sample": sample}
            for src in sources.PER_SAMPLE_SOURCES:
                record.update(src.parse(analysis_dir, sample, self.config))
            records.append(record)
        per_sample_df = pd.DataFrame(records)
        return sources.compute_tree_pass(per_sample_df, self.config)

    def _per_run_sources(self):
        active = list(sources.PER_RUN_SOURCES)
        if self.config.get('auto-nfflu', False):
            active.append(sources.PIPELINE_STATUS_SOURCE)
        return active

    def collect_run_summary(self, analysis_dir, output_summary_file, *, run_id=None, sample_ids=None):
        """Collect all results and merge them into a single summary CSV.

        `run_id` and `sample_ids` are normally derived from `analysis_dir`
        (its path and contents), matching the original behavior of this
        method. A caller that already knows both (e.g. auto-nfflu, which
        orchestrated the run) can pass them directly instead, skipping
        that re-derivation.
        """
        output_summary_file = output_summary_file.rstrip(os.sep)

        if sample_ids is None:
            sample_ids = self._resolve_sample_ids(analysis_dir)

        if len(sample_ids) == 0:
            logging.warning(json.dumps({"event_type": "no_samples_found", "analysis_dir": analysis_dir}))
            return

        if run_id is None:
            run_id = self._resolve_run_id(analysis_dir)

        output_df = self._build_sample_identity(sample_ids, run_id)

        per_sample_df = self._collect_per_sample(analysis_dir, sample_ids)
        output_df = output_df.merge(per_sample_df, on='sample', how='left')

        for src in self._per_run_sources():
            src_df = src.parse(analysis_dir, self.config)
            if src_df is not None and not src_df.empty:
                output_df = output_df.merge(src_df, on='sample', how='left')

        output_df = sources.compute_subtype_status(output_df)

        for src in sources.GLOBAL_SOURCES:
            output_df = output_df.assign(**src.parse(analysis_dir, self.config))

        if "sample" in output_df.columns:
            output_df = output_df.drop(columns=['sample'])

        expected_columns = self.config.get('expected_columns', schema.CANONICAL_COLUMNS)
        output_df = schema.order_and_validate(output_df, expected_columns)

        output_summary_dir = os.path.dirname(output_summary_file)
        if output_summary_dir != '' and not os.path.exists(output_summary_dir):
            os.makedirs(output_summary_dir, exist_ok=True)

        output_df.to_csv(output_summary_file, index=False)
        logging.info(json.dumps({"event_type": "results_written", "output_file": output_summary_file}))

    def collect_nextclade_results(self, analysis_dir, nextclade_output_path):
        nextclade_rc = nextclade.Nextclade_Results_Collector(self.config)

        logging.debug(json.dumps({"event_type": "collecting_nextclade_results", "config": self.config}))

        nextclade_df = nextclade_rc.collect_nextclade_results(analysis_dir)

        # collect_nextclade_results returns None when no sample produced a
        # usable nextclade call (e.g. nextclade skipped, or no HA segments).
        # Fall back to an empty frame so we always write a (headerless)
        # file rather than crashing on None.to_csv -- important because
        # auto-nfflu calls this before collect_run_summary, so a crash here
        # would prevent the run summary from ever being written.
        if nextclade_df is None:
            logging.warning(json.dumps({"event_type": "nextclade_results_empty", "analysis_dir": analysis_dir}))
            nextclade_df = pd.DataFrame()

        output_dir = os.path.dirname(nextclade_output_path)
        if output_dir != '' and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        nextclade_df.to_csv(nextclade_output_path, sep='\t', index=False)
        logging.info(json.dumps({"event_type": "nextclade_results_written", "output_path": nextclade_output_path}))

    def collect_mixture_report(self, analysis_dir, output_mixture_file):
        """Collect mixture report data for all samples in the analysis directory and concatenate them into a single output file."""
        mixture_glob = os.path.join(analysis_dir, self.config["paths"]["mixtures"].format(sample="*"))
        mixture_files = glob.glob(mixture_glob)

        dfs = []

        for mixture_path in mixture_files:
            sample_name = os.path.basename(mixture_path).replace('_mixtures.csv', '')

            try:
                df = pd.read_csv(mixture_path)

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

        final_df = pd.concat(dfs, ignore_index=True).rename(columns={'sample_name': 'FastQID'})
        final_df.to_csv(output_mixture_file, index=False)
        logging.info(json.dumps({"event_type": "mixture_report_written", "output_file": output_mixture_file}))

        return

    def symlink_consensus_fastas(self, analysis_dir, output_dir):
        """Symlink consensus FASTA files to the output directory."""
        os.makedirs(output_dir, exist_ok=True)

        pattern = os.path.join(analysis_dir, self.config["paths"]["consensus_fasta"].format(sample="*"))

        for fasta_file in glob.glob(pattern):
            basename = os.path.basename(fasta_file)
            dest = os.path.join(output_dir, basename)

            if os.path.exists(dest) or os.path.islink(dest):
                os.remove(dest)

            os.symlink(os.path.abspath(fasta_file), dest)
