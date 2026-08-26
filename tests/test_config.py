import json

import nfflu_results_collector.config as config
import nfflu_results_collector.schema as schema


def test_load_default_config_has_expected_shape():
    cfg = config.load_default_config()
    assert cfg["segments"] == ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
    assert cfg["tree_pass"]["threshold"] == 90
    assert cfg["nextclade"]["ha_only"] is True
    assert cfg["nextclade"]["legacy_clade"] is False
    assert cfg["expected_columns"] == schema.CANONICAL_COLUMNS


def test_overrides_deep_merge_without_clobbering_sibling_keys():
    cfg = config.load_config(overrides={"nextclade": {"ha_only": False}})
    assert cfg["nextclade"]["ha_only"] is False
    # legacy_clade wasn't touched by the override -- deep merge must not
    # wipe it out the way the old flat dict.update() would have.
    assert cfg["nextclade"]["legacy_clade"] is False


def test_legacy_flat_keys_are_normalized_with_deprecation_warning(caplog):
    import logging
    with caplog.at_level(logging.WARNING):
        cfg = config.load_config(overrides={"legacy-clade": True, "nextclade-ha-only": False})
    assert cfg["nextclade"]["legacy_clade"] is True
    assert cfg["nextclade"]["ha_only"] is False
    assert "legacy-clade" not in cfg
    assert "nextclade-ha-only" not in cfg
    events = [json.loads(r.message)["event_type"] for r in caplog.records]
    assert events.count("deprecated_config_key") == 2


def test_explicit_nested_key_wins_over_legacy_flat_key():
    cfg = config.load_config(overrides={"legacy-clade": True, "nextclade": {"legacy_clade": False}})
    assert cfg["nextclade"]["legacy_clade"] is False


def test_load_config_from_json_file(tmp_path):
    path = tmp_path / "user_config.json"
    path.write_text(json.dumps({"tree_pass": {"threshold": 80}}))
    cfg = config.load_config(user_config_path=str(path))
    assert cfg["tree_pass"]["threshold"] == 80
    # untouched defaults survive the merge
    assert cfg["segments"] == ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]


def test_load_config_from_yaml_file(tmp_path):
    path = tmp_path / "user_config.yaml"
    path.write_text("tree_pass:\n  threshold: 75\n")
    cfg = config.load_config(user_config_path=str(path))
    assert cfg["tree_pass"]["threshold"] == 75


def test_overrides_apply_after_and_on_top_of_config_file(tmp_path):
    path = tmp_path / "user_config.json"
    path.write_text(json.dumps({"tree_pass": {"threshold": 80}}))
    cfg = config.load_config(user_config_path=str(path), overrides={"tree_pass": {"threshold": 70}})
    assert cfg["tree_pass"]["threshold"] == 70
