# nfflu-results-collector

A Python module for collecting, aggregating, and summarizing results from [nf-flu](https://github.com/CFIA-NCFAD/nf-flu) influenza analysis pipeline runs.

## Overview

`nfflu-results-collector` extracts and merges data from various output files produced by nf-flu pipeline runs into a consolidated summary CSV file. The module handles multiple data types including:

- Sample metadata and identifiers
- Subtyping results (HA/NA)
- Nextclade phylogenetic classification
- Read mapping statistics (idxstats)
- Consensus sequence completeness
- HPAI cleavage site analysis
- GenoFLU genotyping results
- Software version provenance
- Mixture analysis reports

## Features

- **Automated Data Collection**: Extracts results from standardized nf-flu output directory structure
- **Multi-segment Analysis**: Processes all 8 influenza genome segments (PB2, PB1, PA, HA, NP, NA, M, NS)
- **Quality Metrics**: Calculates consensus completeness and tree-pass status per segment
- **Flexible Configuration**: Supports custom column ordering and auto-nfflu integration
- **Logging**: Comprehensive logging for debugging and tracking data collection progress
- **Mixture Reporting**: Optional mixture analysis report generation
- **Consensus Linking**: Creates symlinks to consensus FASTA files for downstream analysis

## Installation

### From Source

```bash
git clone https://github.com/BCCDC-PHL/nfflu-results-collector.git
cd nfflu-results-collector
pip install -e .
```

### Requirements

- Python >= 3.8
- pandas
- biopython
- pyyaml

## Usage

### Command Line Interface

Basic usage:

```bash
nfflu-results-collector \
  -d /path/to/nfflu/analysis/output \
  -o /path/to/output/summary.csv
```

With all options:

```bash
nfflu-results-collector \
  -d /path/to/nfflu/analysis/output \
  -o /path/to/output/summary.csv \
  -O /path/to/output/mixture_report.csv \
  -s /path/to/output/consensus_symlinks \
  --auto-nfflu \
  --log-level DEBUG
```

### Arguments

- `-d, --analysis-dir`: **(Required)** Path to nf-flu analysis output directory
- `-o, --output-summary`: **(Required)** Path for output summary CSV file
- `-O, --output-mixture`: Path for optional mixture report CSV file
- `-s, --output-symlinks`: Directory path for creating consensus FASTA symlinks
- `-a, --auto-nfflu`: Enable auto-nfflu mode (for automated pipeline integration)
- `--log-level`: Logging verbosity (default: INFO; options: DEBUG, INFO, WARNING, ERROR)

### Python API

```python
from nfflu_results_collector.collector import Nfflu_Results_Collector

# Initialize collector
collector = Nfflu_Results_Collector({'auto-nfflu': False})

# Collect summary results
collector.collect_run_summary(
    analysis_dir='/path/to/nfflu/output',
    output_summary_file='/path/to/summary.csv'
)

# Optional: Collect mixture report
collector.collect_mixture_report(
    analysis_dir='/path/to/nfflu/output',
    output_mixture_file='/path/to/mixture.csv'
)

# Optional: Create symlinks to consensus sequences
collector.symlink_consensus_fastas(
    analysis_dir='/path/to/nfflu/output',
    output_dir='/path/to/symlinks'
)
```

## Input Directory Structure

The collector expects the nf-flu output directory to have the following structure:

```
analysis_output/
├── fastq/                          # Input FASTQ files (for sample name extraction)
├── subtyping_report/
│   └── subtype_results.csv         # Subtyping results
├── nextclade/                      # Nextclade results by sample
│   └── {sample}/
│       └── *.nextclade.tsv
├── mapping/                        # Read mapping statistics
│   └── {sample}/
│       └── {sample}*.idxstats
├── consensus/
│   └── bcftools/
│       └── {sample}.consensus.fasta
├── annotation/
│   └── {sample}/
│       └── {sample}.cleavage.tsv   # HPAI cleavage site
├── genoflu/
│   └── {sample}.genoflu.tsv        # GenoFLU genotyping
├── mixtures/
│   └── {sample}/
│       └── {sample}_mixtures.csv
└── pipeline_info/
    └── software_versions.yml       # Software provenance
```

## Output Format

### Summary CSV

The primary output is a CSV file with one row per sample containing:

**Sample Identifiers:**
- `FastQID`: Full sample identifier from FASTQ filename
- `CID`: Client ID
- `Plate`: Plate identifier
- `Index`: Index number
- `Well`: Well position
- `Run`: Sequencing run identifier

**Subtyping:**
- `subtype_HA_NA_status`: Description of subtyping success
- `HA_subtype`: Hemagglutinin subtype (e.g., H1, H3)
- `NA_subtype`: Neuraminidase subtype (e.g., N1, N2)
- `subtype`: Combined subtype (e.g., H1N1)

**Nextclade:**
- `Nextclade_clade`: Phylogenetic clade
- `Nextclade_subclade`: Phylogenetic subclade
- `Nextclade_qc.overallScore`: Quality score
- `Nextclade_qc.overallStatus`: Quality status
- `nextclade_dataset_name`: Dataset used for classification
- `nextclade_dataset_version`: Dataset version

**Per-Segment Metrics (for each of PB2, PB1, PA, HA, NP, NA, M, NS):**
- `{segment}_reads_mapped`: Number of reads mapped to segment
- `{segment}_seq_length`: Reference sequence length
- `{segment}_consensus_completeness`: Percentage completeness (non-N bases)
- `{segment}_tree_pass`: Binary flag (1 if completeness > 90%, else 0)

**HPAI Cleavage Site:**
- `HPAI_cleave_start`: Start position of cleavage site
- `HPAI_cleave_end`: End position of cleavage site
- `HPAI_cleavage_site_motif`: Amino acid sequence of cleavage site

**GenoFLU Genotyping:**
- `GenoFLU_Genotype`: Overall genotype classification
- `GenoFLU_{segment}`: Genotype for each segment

**Software Versions:**
- `genoflu_version`: GenoFLU version used
- `nextclade_version`: Nextclade version used
- `nfflu_version`: nf-flu pipeline version

### Mixture Report CSV

Optional output containing mixture analysis results:
- `FastQID`: Sample identifier
- `mixture_present`: Boolean indicating mixture detection
- `ha_mixture_present`: HA segment mixture detection
- `na_mixture_present`: NA segment mixture detection
- `ha_read_ratio`: Read ratio for HA segment
- `na_read_ratio`: Read ratio for NA segment

## Configuration

Default configuration is loaded from `nfflu_results_collector/config/defaults.json`. Users can override settings programmatically:

```python
collector = Nfflu_Results_Collector({
    'auto-nfflu': True,
    'expected_columns': ['FastQID', 'CID', 'Plate', ...]
})
```

## Logging

The module uses Python's built-in logging framework. Set log level via command line:

```bash
nfflu-results-collector -d /path/to/data -o output.csv --log-level DEBUG
```

Log messages include:
- Sample collection progress
- File parsing status
- Missing data warnings
- Error details for failed operations

## Development

### Project Structure

```
nfflu-results-collector/
├── nfflu_results_collector/
│   ├── __init__.py
│   ├── __main__.py          # CLI entry point
│   ├── collector.py         # Main collection logic
│   ├── nextclade.py         # Nextclade-specific collection
│   ├── tools.py             # Utility functions
│   ├── config.py            # Configuration management
│   ├── auto.py              # auto-nfflu integration
│   └── config/
│       └── defaults.json    # Default configuration
├── pyproject.toml           # Package metadata
└── README.md
```

### Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes with appropriate logging
4. Test with real nf-flu output data
5. Submit a pull request

## License

[Specify license]

## Authors

- John Palmer (john.palmer@bccdc.ca)

## Links

- [nf-flu Pipeline](https://github.com/CFIA-NCFAD/nf-flu)
- [GitHub Repository](https://github.com/BCCDC-PHL/nfflu-results-collector)