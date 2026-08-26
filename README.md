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
- **Pluggable result sources**: each result type (subtype, nextclade, idxstats, ...) is a small, independently testable parser registered in `sources.py` -- adding a new one doesn't require touching the orchestration logic
- **Config-driven**: file paths, segment list, tree-pass threshold, and Nextclade behaviour all live in config (packaged defaults, optionally overridden by a YAML/JSON file and/or a caller-supplied dict), not hardcoded across the codebase
- **Frozen output schema**: the summary CSV's column set and order (`nfflu_results_collector.schema.CANONICAL_COLUMNS`) is validated on every run; unexpected columns are never silently dropped
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
- pandas >= 2.0
- biopython >= 1.80
- pyyaml >= 6.0

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
  -n /path/to/output/nextclade.tsv \
  -c /path/to/config.yaml \
  --auto-nfflu \
  --log-level DEBUG
```

### Arguments

- `-d, --analysis-dir`: **(Required)** Path to nf-flu analysis output directory
- `-o, --output-summary`: **(Required)** Path for output summary CSV file
- `-O, --output-mixture`: Path for optional mixture report CSV file
- `-s, --output-symlinks`: Directory path for creating consensus FASTA symlinks
- `-n, --output-nextclade`: Path for optional aggregated Nextclade TSV export
- `-a, --auto-nfflu`: Enable auto-nfflu mode (sample discovery from `start_samplesheet.csv`; merges in `pipeline_status.csv` if present)
- `-c, --config`: Path to a user config file (YAML or JSON) overriding packaged defaults -- see [Configuration](#configuration)
- `--log-level`: Logging verbosity (default: INFO; options: DEBUG, INFO, WARNING, ERROR)

### Python API

```python
from nfflu_results_collector.collector import Nfflu_Results_Collector

# Initialize collector (optionally with a config file and/or dict overrides)
collector = Nfflu_Results_Collector({'auto-nfflu': False}, config_path='config.yaml')

# Collect summary results. run_id/sample_ids are normally derived from
# analysis_dir, but a caller that already knows both (e.g. an orchestrator
# that just ran the pipeline) can pass them directly:
collector.collect_run_summary(
    analysis_dir='/path/to/nfflu/output',
    output_summary_file='/path/to/summary.csv',
    run_id='250101_M00123_0001_000000000-ABCDE',   # optional
    sample_ids=['SAMPLE-001', 'SAMPLE-002'],        # optional
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

The collector expects the nf-flu output directory to have the following structure (each path is configurable -- see [Configuration](#configuration)):

Results are published per sample, one directory per sample, with run-level
aggregates alongside them:

```
analysis_output/
├── {sample}/
│   ├── fastq/
│   │   └── {sample}*.merged.fastq.gz   # used for sample name extraction
│   ├── mapping/
│   │   └── {sample}*.idxstats          # Read mapping statistics
│   ├── consensus/
│   │   └── bcftools/
│   │       └── {sample}.consensus.fasta
│   ├── annotation/
│   │   └── {sample}.cleavage.tsv       # HPAI cleavage site
│   ├── genoflu/
│   │   └── {sample}.genoflu.tsv        # GenoFLU genotyping
│   ├── nextclade/
│   │   └── *.nextclade.tsv
│   └── mixtures/
│       └── {sample}_mixtures.csv
├── aggregate/
│   ├── bcftools/
│   │   └── subtyping_report/
│   │       └── subtype_results.csv     # Subtyping results
│   └── nextclade/
│       └── nextclade-dataset-versions.csv  # dataset name/version per sample (optional)
├── pipeline_status.csv                 # per-pipeline-stage status, written by the
│                                        # orchestrator (e.g. auto-nfflu); optional,
│                                        # only merged in when auto-nfflu mode is on
└── pipeline_info/
    └── software_versions.yml           # Software provenance
```

Note the two `nextclade` locations: the per-sample TSVs are published under each
sample, while `nextclade-dataset-versions.csv` is a single run-level file written
by the orchestrator under `aggregate/`.

Note: this collector no longer reaches into Nextflow execution internals (logs,
work directories) for anything. `nextclade-dataset-versions.csv` is expected to
already be published at the path above by whatever orchestrated the pipeline run;
if it's absent, dataset name/version are reported as `'N/A'`.

## Output Format

### Summary CSV

The primary output is a CSV file with one row per sample. Its column set and
order (`nfflu_results_collector.schema.CANONICAL_COLUMNS`) is a frozen contract
consumed by downstream ingestion -- it's pinned by a golden-output test
(`tests/test_golden_output.py`) so refactors can't accidentally change it. If
`--auto-nfflu` is used and a `pipeline_status.csv` is found, dynamically-named
`status_*` columns (one per upstream pipeline) are appended after the frozen
columns.

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

**Per-Segment Metrics** (for each of HA, NA, M, NP, NS, PA, PB1, PB2 -- this
display order, not genome order, is what the report has always used):
- `{segment}_reads_mapped`: Number of reads mapped to segment
- `{segment}_seq_length`: Reference sequence length
- `{segment}_consensus_completeness`: Percentage completeness (non-N bases)
- `{segment}_tree_pass`: Binary flag (1 if completeness > threshold, default 90%, else 0)

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

Configuration is layered: packaged defaults (`nfflu_results_collector/config/defaults.json`)
are deep-merged with an optional user file (`-c/--config`, YAML or JSON, selected
by extension), then deep-merged with any caller-supplied overrides (e.g. the dict
passed to `Nfflu_Results_Collector(...)`, or `--auto-nfflu` on the CLI).
`config-template.json` is a copy of the packaged defaults, kept in sync by a
test; copy it and override only the keys you need.

```yaml
# example config.yaml
tree_pass:
  threshold: 85           # default: 90
segments: [PB2, PB1, PA, HA, NP, NA, M, NS]
nextclade:
  ha_only: true           # default: true
  legacy_clade: false     # default: false
paths:
  idxstats: "{sample}/mapping/{sample}*.idxstats"
  # ... see defaults.json for the full set of path templates
```

For backward compatibility, the pre-nested flat keys (`auto-nfflu`,
`legacy-clade`, `nextclade-ha-only`) are still accepted anywhere a config layer
is provided; they're normalized onto their nested locations with a deprecation
warning logged. Prefer the nested keys in new configs.

```python
collector = Nfflu_Results_Collector({
    'auto-nfflu': True,
    'nextclade': {'ha_only': False},
})
```

## Logging

The module uses Python's built-in logging framework, emitting structured
JSON messages. Set log level via command line:

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
│   ├── collector.py         # Thin orchestrator: resolves sample_ids/run_id,
│   │                         # runs registered sources, applies computed
│   │                         # columns, validates against the schema, writes CSV
│   ├── sources.py           # ResultSource registry: one parser per result
│   │                         # type (subtype, idxstats, completeness, ...)
│   ├── sample_id.py         # FastQID -> {CID, Plate, Index, Well} parsing
│   ├── schema.py            # CANONICAL_COLUMNS (the frozen output contract)
│   │                         # and order_and_validate()
│   ├── nextclade.py         # Nextclade-specific collection/filtering
│   ├── auto.py               # auto-nfflu sample discovery + pipeline_status.csv merge
│   ├── tools.py              # Shared utility functions
│   ├── config.py             # load_config(): layered defaults/file/overrides
│   └── config/
│       └── defaults.json     # Packaged default configuration
├── tests/
│   ├── fixture_builder.py    # Builds a synthetic nf-flu output directory
│   ├── regenerate_golden.py  # Maintenance script for the golden-output test
│   └── fixtures/golden_run_summary.csv
├── pyproject.toml            # Package metadata + dependencies
└── README.md
```

### Running tests

```bash
pip install -e ".[test]"
pytest
```

`tests/test_golden_output.py` is the most important test in the suite: it
byte-compares `collect_run_summary`'s output against a checked-in golden CSV,
guarding the frozen output contract. If you intentionally change what columns
are produced, regenerate it deliberately and review the diff:

```bash
python -m tests.regenerate_golden
git diff tests/fixtures/golden_run_summary.csv
```

### Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes with appropriate logging and tests
4. Test with real nf-flu output data
5. Submit a pull request

## License

[Specify license]

## Authors

- John Palmer (john.palmer@bccdc.ca)

## Links

- [nf-flu Pipeline](https://github.com/CFIA-NCFAD/nf-flu)
- [GitHub Repository](https://github.com/BCCDC-PHL/nfflu-results-collector)
