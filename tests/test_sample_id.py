import pandas as pd
import pytest

from nfflu_results_collector.sample_id import parse_sample_id


@pytest.mark.parametrize(
    "sample_name,expected",
    [
        # standard: CID-Plate-Index-Well
        (
            "C0000000001-2024-A1-A01",
            {"CID": "C0000000001", "Plate": "2024", "Index": "A1", "Well": "A01"},
        ),
        # control: labelled CID-label-Plate-Index-Well
        (
            "NTC20240101-Control-2024-B2-B02",
            {"CID": "NTC20240101", "Plate": "2024", "Index": "B2", "Well": "B02"},
        ),
        # salvage, full fallback: neither standard nor control CID pattern
        # matches (lowercase leading letter), so CID falls back to the
        # whole sample name, but the Plate/Index/Well suffix still parses.
        (
            "c0000000003-2024-D4-D04",
            {
                "CID": "c0000000003-2024-D4-D04",
                "Plate": "2024",
                "Index": "D4",
                "Well": "D04",
            },
        ),
        # salvage, control-style CID recovered (elif branch): fails the
        # strict control format (no "-{label}-" segment) but the
        # letters+8digits CID prefix and the trailing suffix both salvage.
        (
            "NTC201401039999-9999-Z9-Z99",
            {"CID": "NTC20140103", "Plate": "9999", "Index": "Z9", "Well": "Z99"},
        ),
        # salvage, standard-style CID recovered, extra trailing digits
        # truncated to the first 10 (a common real-world data-entry slip).
        (
            "C000000000123-2024-A1-A01",
            {"CID": "C0000000001", "Plate": "2024", "Index": "A1", "Well": "A01"},
        ),
    ],
)
def test_parse_sample_id_known_formats(sample_name, expected):
    assert parse_sample_id(sample_name) == expected


def test_parse_sample_id_full_fallback_has_no_recoverable_fields():
    result = parse_sample_id("unrecognizable_sample_42")
    assert result["CID"] == "unrecognizable_sample_42"
    assert pd.isna(result["Plate"])
    assert pd.isna(result["Index"])
    assert pd.isna(result["Well"])
