
from glob import glob
import json
import logging
import os
from Bio import SeqIO
import pandas as pd 
import re


pd.set_option('future.no_silent_downcasting', True)

def calculate_completeness(fasta_path):
    """Calculate completeness of sequences in a FASTA file."""

    segments = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']

    completeness = dict(zip(segments, [0]*len(segments)))

    for record in SeqIO.parse(fasta_path, "fasta"):
        segment = record.id.split('_')[-1]

        if len(record.seq) > 0:
            n_count = str(record.seq).upper().count('N')
            if segment in completeness:
                completeness[segment] = round(100 - (n_count * 100.0 / len(record.seq)), 2)
            else:
                logging.warning(json.dumps({"event_type": "unknown_segment_found", "segment": segment, "record_id": record.id}))
    df = pd.DataFrame(completeness, index=[0])
    df.columns = [f'{col}_consensus_completeness' for col in df.columns]

    return df

def glob_single(glob_expr):
    """
    Return a single file matching a glob expression.

    Returns:
        str: The file path string if exactly one file matches.
        'none': If no files match.
        'multiple': If more than one file matches.
    """
    files = glob(glob_expr)
    if len(files) == 0:
        logging.warning(json.dumps({
            "event_type": "no_files_found",
            "glob_expression": glob_expr
        }))
        return None
    if len(files) > 1:
        logging.warning(json.dumps({
            "event_type": "multiple_files_found",
            "glob_expression": glob_expr,
            "file_count": len(files)
        }))
        return None
    return files[0].rstrip(os.sep)


def collect_nfflu_fastq_names(analysis_dir):
    """Extract sample names from fastq directory, handling both Nanopore and Illumina paired-end naming."""
    analysis_dir = os.path.abspath(analysis_dir)
    sample_names = set()  # Using a set to avoid duplicates from paired-end reads
    fastq_dir = os.path.join(analysis_dir, "fastq")
    
    if not os.path.exists(fastq_dir):
        return []
    
    pattern = os.path.join(fastq_dir, "*.merged.fastq.gz")
    for filepath in glob(pattern):
        filename = os.path.basename(filepath)
        
        # Handle Illumina paired-end reads with _*1 or _*2 pattern
        if '_1.merged.fastq.gz' in filename or '_2.merged.fastq.gz' in filename:
            # This is an Illumina paired-end read
            # Extract everything before the last underscore
            sample_name = filename.rsplit('_', 1)[0]
            sample_names.add(sample_name)
        else:
            # Assume this is a Nanopore read or another format
            sample_name = filename.replace(".merged.fastq.gz", "")
            sample_names.add(sample_name)

    return sorted(list(sample_names))


def detect_platform(analysis_dir):
    """Return 'nanopore' or 'illumina' for an nf-flu output directory.

    nf-flu writes pipeline_info/samplesheet.fixed.csv under both platforms, but
    with different headers: 'sample,barcode' for nanopore and
    'sample,fastq1,fastq2,single_end' for Illumina. Falls back to 'illumina'
    when the file is missing or unrecognised.
    """
    samplesheet_path = os.path.join(analysis_dir, "pipeline_info", "samplesheet.fixed.csv")

    try:
        with open(samplesheet_path, 'r') as f:
            header = f.readline().strip()
    except OSError:
        logging.warning(json.dumps({"event_type": "platform_samplesheet_not_found", "samplesheet_path": samplesheet_path}))
        return "illumina"

    fields = [field.strip() for field in header.split(',')]
    platform = "nanopore" if "barcode" in fields else "illumina"

    logging.info(json.dumps({"event_type": "platform_detected", "platform": platform, "header": header}))

    return platform


def extract_workdirs(log_file_path: str) -> dict:
    """
    Extract all workDirs from a Nextflow log file, returning a dictionary mapping process names to workDir paths.
    
    Args:
        log_file_path: Path to the Nextflow log file
    
    Returns:
        dict: Mapping of process names to their corresponding workDir paths.
    """
    pattern = r'Task completed.+name: ([^;]+).+workDir: ([^\s]+)'
    workdirs = {}
    
    with open(log_file_path, 'r') as f:
        for line in f:
            match = re.search(pattern, line)
            if match:
                name = match.group(1).strip()
                workdir = match.group(2).strip()
                
                workdirs[name] = workdir
    return workdirs
