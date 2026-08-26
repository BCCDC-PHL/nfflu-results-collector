import os

from nfflu_results_collector.collector import Nfflu_Results_Collector

GOLDEN_PATH = os.path.join(os.path.dirname(__file__), "fixtures", "golden_run_summary.csv")


def test_run_summary_matches_golden_output(analysis_dir, tmp_path):
    """The run_summary.csv column set, order, and values are a frozen
    contract consumed by downstream ingestion. This test pins that
    contract byte-for-byte against a checked-in golden file, so a change
    to it has to be deliberate. See regenerate_golden.py."""
    output_path = tmp_path / "run_summary.csv"

    collector = Nfflu_Results_Collector({"auto-nfflu": False})
    collector.collect_run_summary(str(analysis_dir), str(output_path))

    with open(output_path) as f:
        actual = f.read()
    with open(GOLDEN_PATH) as f:
        expected = f.read()

    assert actual == expected
