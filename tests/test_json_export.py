import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from nfflu_results_collector.collector import Nfflu_Results_Collector, _write_summary_json
from nfflu_results_collector.schema import CANONICAL_COLUMNS
from tests.fixture_builder import SAMPLE_IDS, SAMPLE_STANDARD


def test_export_matches_csv_and_golden(analysis_dir, tmp_path):
    csv_path = tmp_path / "summary.csv"
    json_path = tmp_path / "nested" / "summary.json"
    Nfflu_Results_Collector().collect_run_summary(str(analysis_dir), str(csv_path), output_summary_json_file=json_path)
    assert csv_path.read_bytes() == (Path(__file__).parent / "fixtures/golden_run_summary.csv").read_bytes()
    rows = json.loads(json_path.read_text())
    assert {row["FastQID"] for row in rows} == set(SAMPLE_IDS)
    assert list(rows[0]) == CANONICAL_COLUMNS
    expected = pd.read_csv(csv_path, dtype={field: str for field in ["FastQID", "CID", "Plate", "Index", "Well", "Run"]})
    for actual, record in zip(rows, expected.to_dict("records")):
        for key, value in record.items():
            if pd.isna(value):
                assert actual[key] is None
            elif isinstance(value, (int, float)):
                assert actual[key] == pytest.approx(value)
            else:
                assert actual[key] == value


def test_dynamic_status_and_extra_columns(analysis_dir, tmp_path, monkeypatch):
    (Path(analysis_dir) / "pipeline_status.csv").write_text(f"ID,status_demo\n{SAMPLE_STANDARD},0\n")
    collector = Nfflu_Results_Collector({"auto-nfflu": True})
    original = collector._collect_per_sample
    def with_extra(*args):
        return original(*args).assign(extra_note="001-extra")
    monkeypatch.setattr(collector, "_collect_per_sample", with_extra)
    destination = tmp_path / "summary.json"
    collector.collect_run_summary(str(analysis_dir), str(tmp_path / "summary.csv"), sample_ids=list(SAMPLE_IDS), output_summary_json_file=destination)
    rows = json.loads(destination.read_text())
    assert rows[0]["status_demo"] == 0
    assert rows[1]["status_demo"] is None
    assert all(row["extra_note"] == "001-extra" for row in rows)
    assert list(rows[0])[:len(CANONICAL_COLUMNS)] == CANONICAL_COLUMNS


def test_string_identifiers_nulls_nonfinite_and_zero(tmp_path):
    destination = tmp_path / "summary.json"
    frame = pd.DataFrame({"FastQID": ["001.a-b", "0002"], "CID": ["001", "002"],
                          "zero": [0, 0.0], "missing": [pd.NA, None],
                          "not_finite": [float("inf"), float("nan")],
                          "Nextclade_qc.overallStatus": ["good", None]})
    _write_summary_json(frame, destination)
    rows = json.loads(destination.read_text(), parse_constant=lambda value: pytest.fail(value))
    assert rows[0]["FastQID"] == "001.a-b"
    assert rows[0]["CID"] == "001"
    assert rows[0]["zero"] == 0
    assert rows[0]["Nextclade_qc.overallStatus"] == "good"
    assert all(row["missing"] is None and row["not_finite"] is None for row in rows)


def test_no_samples_keeps_csv_and_writes_empty_json(tmp_path):
    csv_path = tmp_path / "summary.csv"
    csv_path.write_text("existing csv")
    json_path = tmp_path / "summary.json"
    Nfflu_Results_Collector().collect_run_summary(str(tmp_path), str(csv_path), sample_ids=[], output_summary_json_file=json_path)
    assert json.loads(json_path.read_text()) == []
    assert csv_path.read_text() == "existing csv"


@pytest.mark.parametrize("failure", ["serialize", "replace"])
def test_failure_preserves_destination_and_cleans_temporary(analysis_dir, tmp_path, monkeypatch, failure):
    destination = tmp_path / "json" / "summary.json"
    destination.parent.mkdir()
    destination.write_text('[{"previous":true}]')
    def fail(*args, **kwargs):
        if failure == "serialize":
            args[1].write("partial")
        raise OSError("simulated write failure")
    if failure == "serialize":
        monkeypatch.setattr(pd.DataFrame, "to_json", fail)
    else:
        monkeypatch.setattr("nfflu_results_collector.collector.os.replace", fail)
    with pytest.raises(OSError, match="simulated"):
        Nfflu_Results_Collector().collect_run_summary(str(analysis_dir), str(tmp_path / "summary.csv"), output_summary_json_file=destination)
    assert destination.read_text() == '[{"previous":true}]'
    assert list(destination.parent.iterdir()) == [destination]


def test_cli_exports_json(analysis_dir, tmp_path):
    destination = tmp_path / "summary.json"
    subprocess.run([sys.executable, "-m", "nfflu_results_collector", "-d", str(analysis_dir),
                    "-o", str(tmp_path / "summary.csv"), "--output-summary-json", str(destination)], check=True, capture_output=True)
    assert len(json.loads(destination.read_text())) == len(SAMPLE_IDS)
