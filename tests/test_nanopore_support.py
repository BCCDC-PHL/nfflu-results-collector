import csv
import os

import pytest

from nfflu_results_collector.collector import MIXTURE_COLUMNS, Nfflu_Results_Collector
from nfflu_results_collector.tools import detect_platform


ILLUMINA_HEADER = "sample,fastq1,fastq2,single_end"
NANOPORE_HEADER = "sample,barcode"


def _write_samplesheet(analysis_dir, header):
    path = os.path.join(analysis_dir, "pipeline_info", "samplesheet.fixed.csv")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(header + "\nS1,x\n")
    return path


@pytest.mark.parametrize("header,expected", [
    (NANOPORE_HEADER, "nanopore"),
    (ILLUMINA_HEADER, "illumina"),
])
def test_detect_platform_reads_samplesheet_header(tmp_path, header, expected):
    _write_samplesheet(str(tmp_path), header)
    assert detect_platform(str(tmp_path)) == expected


def test_detect_platform_defaults_to_illumina_when_missing(tmp_path):
    assert detect_platform(str(tmp_path)) == "illumina"


def test_mixture_report_written_with_na_rows_when_no_reports(tmp_path):
    output = str(tmp_path / "mixture_report.csv")

    Nfflu_Results_Collector().collect_mixture_report(
        str(tmp_path), output, sample_ids=["S1", "S2"]
    )

    with open(output) as f:
        rows = list(csv.DictReader(f))

    assert [row["FastQID"] for row in rows] == ["S1", "S2"]
    assert list(rows[0]) == ["FastQID"] + [c for c in MIXTURE_COLUMNS if c != "sample_name"]
    assert all(row["mixture_present"] == "" for row in rows)


def test_mixture_report_keeps_samples_without_a_report(tmp_path):
    analysis_dir = str(tmp_path)
    mixtures_path = os.path.join(analysis_dir, "mixtures", "S1", "S1_mixtures.csv")
    os.makedirs(os.path.dirname(mixtures_path))
    with open(mixtures_path, "w") as f:
        f.write(",".join(MIXTURE_COLUMNS) + "\n")
        f.write("S1,H5N1,0," + ",".join(["0"] * (len(MIXTURE_COLUMNS) - 3)) + "\n")

    output = str(tmp_path / "mixture_report.csv")
    Nfflu_Results_Collector().collect_mixture_report(
        analysis_dir, output, sample_ids=["S1", "S2"]
    )

    with open(output) as f:
        rows = list(csv.DictReader(f))

    assert [row["FastQID"] for row in rows] == ["S1", "S2"]
    assert rows[0]["subtype"] == "H5N1"
    assert rows[1]["subtype"] == ""


def test_collect_run_summary_accepts_run_id_and_sample_ids(tmp_path):
    """The call auto-nfflu makes: neither value is derivable from tmp_path."""
    output = str(tmp_path / "run_summary.csv")

    Nfflu_Results_Collector().collect_run_summary(
        str(tmp_path), output, run_id="RUN123", sample_ids=["S1", "S2"]
    )

    with open(output) as f:
        rows = list(csv.DictReader(f))

    assert [row["FastQID"] for row in rows] == ["S1", "S2"]
    assert {row["Run"] for row in rows} == {"RUN123"}
