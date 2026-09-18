import json

import pandas as pd
import pytest

import nfflu_results_collector.config as config
import nfflu_results_collector.schema as schema


def test_order_and_validate_reorders_and_fills_missing_columns():
    df = pd.DataFrame([{"b": 2, "a": 1}])
    result = schema.order_and_validate(df, schema=["a", "b", "c"])
    assert list(result.columns) == ["a", "b", "c"]
    assert result["c"].isna().all()


def test_order_and_validate_appends_undeclared_columns_instead_of_dropping(caplog):
    import logging
    df = pd.DataFrame([{"a": 1, "surprise": 42}])
    with caplog.at_level(logging.WARNING):
        result = schema.order_and_validate(df, schema=["a"])
    assert list(result.columns) == ["a", "surprise"]
    assert result["surprise"].tolist() == [42]
    events = [json.loads(r.message)["event_type"] for r in caplog.records]
    assert "unexpected_columns_found" in events


def test_order_and_validate_warns_on_out_of_range_value_without_changing_it(caplog):
    import logging
    df = pd.DataFrame([{"HA_consensus_completeness": 150.0}])
    with caplog.at_level(logging.WARNING):
        result = schema.order_and_validate(df, schema=["HA_consensus_completeness"])
    assert result["HA_consensus_completeness"].tolist() == [150.0]
    events = [json.loads(r.message)["event_type"] for r in caplog.records]
    assert "column_failed_check" in events


def test_order_and_validate_warns_rather_than_raising_on_a_non_numeric_value(caplog):
    """A check must never take down a run: a column gone to strings gets a
    warning, not a TypeError out of the comparison."""
    import logging
    df = pd.DataFrame([{"M_consensus_completeness": "100%"}])
    with caplog.at_level(logging.WARNING):
        result = schema.order_and_validate(df, schema=["M_consensus_completeness"])
    assert result["M_consensus_completeness"].tolist() == ["100%"]
    events = [json.loads(r.message)["event_type"] for r in caplog.records]
    assert "column_failed_check" in events


def test_order_and_validate_does_not_double_warn_on_a_missing_column(caplog):
    import logging
    df = pd.DataFrame([{"a": 1}])
    with caplog.at_level(logging.WARNING):
        schema.order_and_validate(df, schema=["a", "HA_reads_mapped"])
    events = [json.loads(r.message)["event_type"] for r in caplog.records]
    assert "expected_column_missing" in events
    assert "column_failed_check" not in events


def test_order_and_validate_passes_valid_values(caplog):
    import logging
    df = pd.DataFrame([{"HA_consensus_completeness": 97.0, "HA_reads_mapped": 61000, "HA_tree_pass": 1}])
    with caplog.at_level(logging.WARNING):
        schema.order_and_validate(df, schema=list(df.columns))
    events = [json.loads(r.message)["event_type"] for r in caplog.records]
    assert "column_failed_check" not in events


def test_canonical_columns_have_no_duplicates():
    assert len(schema.CANONICAL_COLUMNS) == len(set(schema.CANONICAL_COLUMNS))


def test_display_segment_order_is_a_permutation_of_processing_segments():
    """DISPLAY_SEGMENT_ORDER and config's "segments" are deliberately different
    orders of the same eight segments; adding one to either alone is a bug."""
    assert set(schema.DISPLAY_SEGMENT_ORDER) == set(config.load_default_config()["segments"])


def test_mixture_columns_match_pileup_tools():
    """MIXTURE_COLUMNS is a copy of the column list pileup-tools writes, kept
    here so the collector doesn't have to depend on it. Skipped where
    pileup-tools isn't installed."""
    mixture_detector = pytest.importorskip("pileup_tools.mixture_detector")
    assert schema.MIXTURE_COLUMNS == mixture_detector.MixtureDetector.COLUMN_NAMES
