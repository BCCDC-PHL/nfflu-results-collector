"""Maintenance script: regenerates tests/fixtures/golden_run_summary.csv by
running the collector against the synthetic fixture.

The golden file is a snapshot of what the collector currently produces for
the synthetic fixture. Only regenerate it deliberately, when a change to
expected_columns or collection logic is an *intentional*, reviewed change to
the output contract -- not as part of routine development, and review the
resulting diff. The golden-output test (test_golden_output.py) exists
specifically to catch *unintentional* drift.

Usage (from repo root):
    python -m tests.regenerate_golden
"""
import os
import tempfile

from nfflu_results_collector.collector import Nfflu_Results_Collector
from tests.fixture_builder import build_fixture

GOLDEN_PATH = os.path.join(os.path.dirname(__file__), "fixtures", "golden_run_summary.csv")


def main():
    with tempfile.TemporaryDirectory() as tmp_dir:
        analysis_dir = build_fixture(tmp_dir)
        collector = Nfflu_Results_Collector({"auto-nfflu": False})
        collector.collect_run_summary(analysis_dir, GOLDEN_PATH)
    print(f"Wrote golden output to {GOLDEN_PATH}")


if __name__ == "__main__":
    main()
