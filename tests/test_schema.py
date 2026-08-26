import json

import pandas as pd
import pytest

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


def test_order_and_validate_strict_raises_on_undeclared_columns():
    df = pd.DataFrame([{"a": 1, "surprise": 42}])
    with pytest.raises(ValueError, match="surprise"):
        schema.order_and_validate(df, schema=["a"], strict=True)


def test_canonical_columns_have_no_duplicates():
    assert len(schema.CANONICAL_COLUMNS) == len(set(schema.CANONICAL_COLUMNS))


def test_display_segment_order_is_a_permutation_of_processing_segments():
    processing_segments = {"PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"}
    assert set(schema.DISPLAY_SEGMENT_ORDER) == processing_segments
