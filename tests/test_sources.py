import pandas as pd
import pytest

import nfflu_results_collector.config as config
import nfflu_results_collector.schema as schema
import nfflu_results_collector.sources as sources
from tests.fixture_builder import SAMPLE_STANDARD, SAMPLE_CONTROL, SAMPLE_SALVAGE, SAMPLE_EXTERNAL


@pytest.fixture
def cfg():
    return config.load_default_config()


def test_registered_source_columns_plus_computed_columns_cover_full_schema():
    """Every column in the frozen schema is accounted for either by a
    registered source or by a named computed-column step. Catches drift
    if a parser gains/loses a field without schema.py being updated, or
    vice versa."""
    computed = (
        set(schema.SAMPLE_IDENTITY_COLUMNS)
        | set(schema.SUBTYPE_COMPUTED_COLUMNS)
        | {f"{seg}_tree_pass" for seg in schema.DISPLAY_SEGMENT_ORDER}
    )
    covered = set(sources.ALL_SOURCE_COLUMNS) | computed
    assert covered == set(schema.CANONICAL_COLUMNS)


def test_idxstats_source_zero_fills_missing_segments(analysis_dir, cfg):
    # S2 (control) only has HA + NA idxstats files.
    result = sources._parse_idxstats(str(analysis_dir), SAMPLE_CONTROL, cfg)
    assert result["HA_reads_mapped"] == 1200
    assert result["NA_reads_mapped"] == 900
    assert result["PB2_reads_mapped"] == 0
    assert result["PB2_seq_length"] == 0


def test_consensus_completeness_matches_designed_n_counts(analysis_dir, cfg):
    result = sources._parse_consensus_completeness(str(analysis_dir), SAMPLE_STANDARD, cfg)
    assert result["HA_consensus_completeness"] == 97.0
    assert result["M_consensus_completeness"] == 100.0


def test_consensus_completeness_defaults_to_zero_when_fasta_missing(analysis_dir, cfg, tmp_path):
    result = sources._parse_consensus_completeness(str(tmp_path), "nonexistent-sample", cfg)
    assert all(v == 0.0 for v in result.values())
    assert set(result.keys()) == {f"{seg}_consensus_completeness" for seg in cfg["segments"]}


def test_tree_pass_threshold_is_strict_greater_than(cfg):
    df = pd.DataFrame([
        {"sample": "s1", "HA_consensus_completeness": 90.0},
        {"sample": "s2", "HA_consensus_completeness": 90.01},
        {"sample": "s3", "HA_consensus_completeness": 89.99},
    ])
    result = sources.compute_tree_pass(df, cfg)
    assert result.set_index("sample")["HA_tree_pass"].to_dict() == {"s1": 0, "s2": 1, "s3": 0}


def test_genotype_parser_remaps_mp_to_m_and_preserves_value_quirk(analysis_dir, cfg):
    # GenoFLU values have a leading space quirk from the original parser
    # ("PB2: A3, ...".split(":") -> " A3"); preserved deliberately since
    # it's part of the frozen output contract.
    result = sources._parse_genotype(str(analysis_dir), SAMPLE_STANDARD, cfg)
    assert result["GenoFLU_M"] == " A3"
    assert result["GenoFLU_Genotype"] == "B3.6"


def test_genotype_parser_defaults_to_none_when_file_missing(analysis_dir, cfg):
    result = sources._parse_genotype(str(analysis_dir), SAMPLE_CONTROL, cfg)
    assert all(v is None for v in result.values())


def test_cleavage_parser_extracts_motif_and_position(analysis_dir, cfg):
    result = sources._parse_cleavage(str(analysis_dir), SAMPLE_STANDARD, cfg)
    assert result["HPAI_cleavage_site_motif"] == "PRRARRVSLVQERG"
    assert result["HPAI_cleave_start"] == 1035
    assert result["HPAI_cleave_end"] == 1061


def test_subtype_parser_derives_prefixed_subtype_columns(analysis_dir, cfg):
    df = sources._parse_subtype(str(analysis_dir), cfg)
    row = df.set_index("sample").loc[SAMPLE_STANDARD]
    assert row["HA_subtype"] == "H5"
    assert row["NA_subtype"] == "N1"
    assert row["subtype"] == "H5N1"


def test_subtype_parser_omits_samples_without_a_call(analysis_dir, cfg):
    df = sources._parse_subtype(str(analysis_dir), cfg)
    assert SAMPLE_CONTROL not in set(df["sample"])


def test_compute_subtype_status_maps_presence_to_description():
    df = pd.DataFrame([
        {"HA_subtype": "H5", "NA_subtype": "N1"},
        {"HA_subtype": "H5", "NA_subtype": pd.NA},
        {"HA_subtype": pd.NA, "NA_subtype": "N1"},
        {"HA_subtype": pd.NA, "NA_subtype": pd.NA},
    ])
    result = sources.compute_subtype_status(df)
    assert result["subtype_HA_NA_status"].tolist() == [
        "HA and NA subtype successful",
        "HA subtype successful",
        "NA subtype successful",
        "No subtype for HA or NA",
    ]


def test_provenance_reads_software_versions(analysis_dir, cfg):
    result = sources._parse_provenance(str(analysis_dir), cfg)
    assert result == {
        "genoflu_version": "1.8.1",
        "nextclade_version": "3.11.0",
        "nfflu_version": "3.10.0",
    }


def test_nextclade_source_picks_ha_only_and_attaches_dataset_metadata(analysis_dir, cfg):
    df = sources._parse_nextclade(str(analysis_dir), cfg)
    row = df.set_index("sample").loc[SAMPLE_SALVAGE]
    assert row["Nextclade_clade"] == "3C.2a1b.2a.2"
    assert row["Nextclade_dataset_name"] == "nextstrain/flu/h3n2/ha/EPI1857216"
    # exactly one row per sample (HA-only filter + alignment-score max)
    assert df["sample"].is_unique


def test_nextclade_source_keeps_sample_ids_containing_underscores(analysis_dir, cfg):
    """Sequence IDs are "<sample>_<n>_<GENE>", so deriving the sample by
    splitting on the first underscore truncated external partner IDs,
    which contain underscores themselves, and every Nextclade column for
    those samples silently came out null."""
    df = sources._parse_nextclade(str(analysis_dir), cfg)
    assert SAMPLE_EXTERNAL in set(df["sample"])
    row = df.set_index("sample").loc[SAMPLE_EXTERNAL]
    assert row["Nextclade_clade"] == "2.3.4.4b"
    assert row["Nextclade_subclade"] == "2.3.4.4b.5"
