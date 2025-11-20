import os 
import json
import logging
import pandas as pd 
from glob import glob
from nfflu_results_collector.tools import glob_single, collect_nfflu_fastq_names


def collect_auto_nfflu_names(nfflu_results_dir):
    analysis_output_dir = os.path.dirname(nfflu_results_dir.rstrip(os.sep))
    start_df_path = glob_single(os.path.join(analysis_output_dir, 'samplesheets', '*start_samplesheet.csv'))

    if not start_df_path:
        logging.warning(json.dumps({
            "event_type": "start_samplesheet_not_found",
            "samplesheet_dir": os.path.join(analysis_output_dir, 'samplesheets')
        }))
        return []
    
    return pd.read_csv(start_df_path)['ID'].tolist()


def _compute_nfflu_status(results_dir):
    
    required_files = [
        ('irma', os.path.join(results_dir, 'irma', "*.irma.consensus.fasta"), lambda x: x.split(".")[0]),                # 2 
        ('blastn_ref', os.path.join(results_dir, 'blast', 'blastn', 'irma', "*blastn.txt"), lambda x: x.split(".")[0]),  # 3 
        ('reference', os.path.join(results_dir, 'reference_sequences', '*'), lambda x : x),                              # 4 
        ('minimap2', os.path.join(results_dir, 'mapping', "*"), lambda x : x),
        ('freebayes', os.path.join(results_dir, 'variants', "*"), lambda x : x),
        ('bcftools', os.path.join(results_dir, 'consensus', 'bcftools', "*.consensus.fasta"), lambda x: x.split(".")[0]),
        ('blastn_subtype', os.path.join(results_dir, 'blast', 'blastn', 'consensus', "*.blastn.txt"), lambda x: x.split(".")[0]),
        ('mixture', os.path.join(results_dir, 'mixtures', "*", "*_mixtures.txt"), lambda x: x.split("_")[0]),
    ]

    status_codes = {x : n for n, x in enumerate([stage for stage, _, _ in required_files], 2)}

    base_df = pd.DataFrame(collect_nfflu_fastq_names(results_dir), columns=['ID'])

    for stage, glob_expr, transform in required_files:
        search_results = glob(glob_expr)

        if len(search_results) == 0:
            logging.error(json.dumps({
                "event_type": "no_files_found_for_stage",
                "stage": stage,
                "glob_expression": glob_expr
            }))
            continue 

        sample_names = [transform(os.path.basename(result)) for result in search_results]

        # mapping it so that True = pass and False = fail for this stage
        if stage != 'mixture':
            stage_df = pd.DataFrame(sample_names, columns=['ID'])
            stage_df[stage] = True

            base_df = base_df.merge(stage_df, on='ID', how='left')
        else:
            mixture_status = []
            for fp in search_results:
                with open(fp, 'r') as infile:
                    mixture_status.append(infile.read() != '1')
                    
            base_df = base_df.merge(pd.DataFrame({'ID': sample_names, stage: mixture_status}), on='ID', how='left')

    base_df = base_df.set_index('ID')
    base_df = base_df.fillna(value=False).astype(bool)

    # here is where I flip the logic 
    # 'error' will be True if there is a failure, False if no error 
    # 'index' will be the first stage that failed (since failed stages are True now)
    status_df = base_df.apply(lambda row : {'index': (~row).idxmax(), 'error': (~row).max()}, axis=1, result_type='expand')

    base_df['status_nf-flu'] = status_df.apply(lambda x : status_codes[x['index']] if x['error'] else 0, axis=1).tolist()
    base_df = base_df.reset_index()
    
    return base_df[['ID', 'status_nf-flu']]

def compute_pipeline_status_columns(nfflu_results_dir):
    """Combine all samplesheets in a directory into a single DataFrame."""

    analysis_output_dir = os.path.dirname(nfflu_results_dir.rstrip(os.sep))

    samplesheet_paths = glob(os.path.join(analysis_output_dir, 'samplesheets', '*_output_samplesheet.csv'))

    sample_names = collect_auto_nfflu_names(nfflu_results_dir)
    if not sample_names:
        logging.warning(json.dumps({
            "event_type": "no_sample_names_found",
            "nfflu_results_dir": nfflu_results_dir
        }))
        return pd.DataFrame(columns=['ID'])
    main_df = pd.DataFrame(sample_names, columns=['ID'])

    for samplesheet_path in samplesheet_paths:
        try:
            parts = os.path.basename(samplesheet_path).rsplit("_", 3)
            if len(parts) < 2:
                logging.warning(json.dumps({
                    "event_type": "samplesheet_filename_malformed",
                    "samplesheet_path": samplesheet_path,
                    "filename": os.path.basename(samplesheet_path)
                }))
                continue
            pipeline_name = parts[1]
            df = pd.read_csv(samplesheet_path)[['ID']]
            df["status_" + pipeline_name] = 0

            main_df = main_df.merge(df, on='ID', how='left')
        except Exception as e:
            logging.warning(json.dumps({
                "event_type": "samplesheet_read_error",
                "samplesheet_path": samplesheet_path,
                "error": str(e)
            }))

    main_df = main_df.fillna(1)

    nfflu_df = _compute_nfflu_status(nfflu_results_dir)
    main_df = main_df.merge(nfflu_df, on='ID', how='left')
    main_df['status_nf-flu'] = main_df['status_nf-flu'].fillna(1).astype(int)
    
    return main_df