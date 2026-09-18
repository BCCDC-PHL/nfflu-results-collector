import csv
import os

import pandas as pd
import pytest

from nfflu_results_collector.collector import Nfflu_Results_Collector
from tests.fixture_builder import build_fixture, RUN_ID, SAMPLE_IDS, SAMPLE_STANDARD, SAMPLE_CONTROL, SAMPLE_SALVAGE

GOLDEN_PATH = os.path.join(os.path.dirname(__file__), "fixtures", "golden_run_summary.csv")


def test_explicit_run_context_matches_auto_derived_output(analysis_dir, tmp_path):
    """Passing run_id/sample_ids explicitly (the auto-nfflu integration
    path) must produce the same output as letting the collector derive
    them itself (the standalone-CLI path), when both resolve to the same
    values."""
    output_path = tmp_path / "run_summary.csv"
    collector = Nfflu_Results_Collector({"auto-nfflu": False})
    collector.collect_run_summary(
        str(analysis_dir), str(output_path), run_id=RUN_ID, sample_ids=list(SAMPLE_IDS)
    )

    with open(output_path) as f:
        actual = f.read()
    with open(GOLDEN_PATH) as f:
        expected = f.read()
    assert actual == expected


def test_explicit_sample_ids_subset_only_includes_those_samples(analysis_dir, tmp_path):
    output_path = tmp_path / "run_summary.csv"
    collector = Nfflu_Results_Collector({"auto-nfflu": False})
    collector.collect_run_summary(
        str(analysis_dir), str(output_path), run_id=RUN_ID, sample_ids=[SAMPLE_STANDARD]
    )
    df = pd.read_csv(output_path)
    assert df["FastQID"].tolist() == [SAMPLE_STANDARD]


def test_auto_nfflu_mode_merges_pipeline_status_columns(analysis_dir, tmp_path):
    status_path = os.path.join(str(analysis_dir), "pipeline_status.csv")
    with open(status_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "status_cutadapt-nf", "status_basic-sequence-qc", "status_nf-flu"])
        writer.writerow([SAMPLE_STANDARD, 0, 0, 0])
        writer.writerow([SAMPLE_CONTROL, 0, 1, 1])
        writer.writerow([SAMPLE_SALVAGE, 0, 0, 0])

    output_path = tmp_path / "run_summary.csv"
    collector = Nfflu_Results_Collector({"auto-nfflu": True})
    collector.collect_run_summary(
        str(analysis_dir), str(output_path), run_id=RUN_ID, sample_ids=list(SAMPLE_IDS)
    )

    df = pd.read_csv(output_path)
    assert "status_nf-flu" in df.columns
    row = df.set_index("FastQID").loc[SAMPLE_CONTROL]
    assert row["status_basic-sequence-qc"] == 1
    assert row["status_nf-flu"] == 1

    # the 64 declared columns still appear first, in order, ahead of
    # the dynamically-named status_* columns
    from nfflu_results_collector.schema import CANONICAL_COLUMNS
    assert list(df.columns)[: len(CANONICAL_COLUMNS)] == CANONICAL_COLUMNS


def test_auto_nfflu_mode_without_pipeline_status_file_still_succeeds(analysis_dir, tmp_path):
    """pipeline_status.csv is optional -- the collector must not crash if the
    orchestrator didn't write one (e.g. running the CLI by hand). Without it
    there are no status_* columns to report, and the orchestrator owns which
    ones exist, so none are invented here."""
    output_path = tmp_path / "run_summary.csv"
    collector = Nfflu_Results_Collector({"auto-nfflu": True})
    collector.collect_run_summary(
        str(analysis_dir), str(output_path), run_id=RUN_ID, sample_ids=list(SAMPLE_IDS)
    )
    df = pd.read_csv(output_path)
    assert len(df) == len(SAMPLE_IDS)
    assert not [c for c in df.columns if c.startswith("status_")]


def test_collect_mixture_report(analysis_dir, tmp_path):
    """Every sample gets a row, including the two with no mixture report,
    because downstream ingestion reads the file unconditionally."""
    output_path = tmp_path / "mixture_report.csv"
    collector = Nfflu_Results_Collector()
    collector.collect_mixture_report(str(analysis_dir), str(output_path))

    df = pd.read_csv(output_path)
    assert set(df["FastQID"]) == set(SAMPLE_IDS)
    indexed = df.set_index("FastQID")
    assert bool(indexed.loc[SAMPLE_STANDARD, "mixture_present"]) is True
    assert pd.isna(indexed.loc[SAMPLE_CONTROL, "mixture_present"])


def test_symlink_consensus_fastas(analysis_dir, tmp_path):
    output_dir = tmp_path / "symlinks"
    collector = Nfflu_Results_Collector()
    collector.symlink_consensus_fastas(str(analysis_dir), str(output_dir))

    linked = sorted(os.listdir(output_dir))
    assert linked == sorted(f"{s}.consensus.fasta" for s in SAMPLE_IDS)
    for name in linked:
        assert os.path.islink(output_dir / name)


def test_collect_nextclade_results_writes_tsv(analysis_dir, tmp_path):
    output_path = tmp_path / "nextclade.tsv"
    collector = Nfflu_Results_Collector()
    collector.collect_nextclade_results(str(analysis_dir), str(output_path))

    df = pd.read_csv(output_path, sep="\t")
    assert set(df["sample"]) == set(SAMPLE_IDS)


def test_collect_nextclade_results_handles_empty_nextclade_without_crashing(tmp_path):
    """When nextclade produced no usable calls (empty/absent nextclade
    dir), collect_nextclade_results must not raise -- auto-nfflu calls it
    before collect_run_summary, so a crash here would block the summary."""
    analysis_dir = tmp_path / "empty_analysis"
    (analysis_dir / "nextclade").mkdir(parents=True)
    output_path = tmp_path / "nested" / "nextclade.tsv"  # parent dir does not exist yet

    collector = Nfflu_Results_Collector()
    collector.collect_nextclade_results(str(analysis_dir), str(output_path))

    assert output_path.exists()


def test_collect_run_summary_creates_missing_output_directory(analysis_dir, tmp_path):
    output_path = tmp_path / "does" / "not" / "exist" / "run_summary.csv"
    collector = Nfflu_Results_Collector()
    collector.collect_run_summary(str(analysis_dir), str(output_path))
    assert output_path.exists()


def test_both_layouts_produce_the_same_run_summary(tmp_path):
    """The per-sample layout auto-nfflu migrates to next is exercised here
    against the same fixture data, so switching config["layout"] is a
    tested change rather than an untried one."""
    outputs = {}
    for layout in ("stage", "sample"):
        analysis_dir = build_fixture(tmp_path / layout, layout=layout)
        output_path = tmp_path / f"{layout}_run_summary.csv"
        collector = Nfflu_Results_Collector({"layout": layout})
        collector.collect_run_summary(str(analysis_dir), str(output_path), run_id=RUN_ID)
        outputs[layout] = output_path.read_text()

    assert outputs["stage"] == outputs["sample"]
