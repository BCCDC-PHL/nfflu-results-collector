import json
import logging
import os

import yaml

import nfflu_results_collector.schema as schema

_PACKAGE_DEFAULTS_PATH = os.path.join(os.path.dirname(__file__), "config", "defaults.json")

# legacy flat key -> (nested section, nested key). Kept for backward
# compatibility with existing caller configs (e.g. auto-nfflu's
# submodule_params.nfflu_results_collector block) that predate the nested
# config layout.
_LEGACY_KEY_MAP = {
    "nextclade-ha-only": ("nextclade", "ha_only"),
    "legacy-clade": ("nextclade", "legacy_clade"),
}


def _load_packaged_defaults():
    with open(_PACKAGE_DEFAULTS_PATH, "r") as f:
        config = json.load(f)
    config.setdefault("expected_columns", schema.CANONICAL_COLUMNS)
    return config


def _deep_merge(base, override):
    result = dict(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def _normalize_legacy_keys(layer):
    """Map legacy flat keys onto their new nested locations, with a
    deprecation warning. An explicit new-style nested key within the SAME
    layer always wins if both are present.

    Must be applied to each config layer (packaged defaults, user file,
    overrides) individually, BEFORE that layer is deep-merged onto the
    others. Applying it once to the final merged config would be wrong:
    the packaged defaults always define the nested keys, so by the time
    a later layer's legacy flat key is normalized against the
    already-merged dict, the "new_key not in section_dict" check would
    always find the default there first and the legacy key would never
    take effect.
    """
    layer = dict(layer)
    for legacy_key, (section, new_key) in _LEGACY_KEY_MAP.items():
        if legacy_key not in layer:
            continue
        section_dict = dict(layer.get(section, {}))
        if new_key not in section_dict:
            section_dict[new_key] = layer[legacy_key]
        logging.warning(json.dumps({
            "event_type": "deprecated_config_key",
            "legacy_key": legacy_key,
            "new_key": f"{section}.{new_key}",
        }))
        layer[section] = section_dict
        del layer[legacy_key]
    return layer


def load_default_config():
    """Backward-compatible alias for the packaged defaults (no user
    overrides applied)."""
    return _load_packaged_defaults()


def load_config(user_config_path=None, overrides=None):
    """Build a config dict: packaged defaults, deep-merged with an
    optional user config file (YAML or JSON, selected by extension),
    deep-merged with an optional `overrides` dict (e.g. CLI flags or a
    caller-supplied dict, as in `Nfflu_Results_Collector(overrides)`).
    Legacy flat keys (auto-nfflu's pre-existing submodule_params, etc.)
    are normalized onto their nested locations within each layer before
    that layer is merged in, so a legacy key in a higher-precedence layer
    (e.g. overrides) still takes effect over a nested key set by a
    lower-precedence layer (e.g. packaged defaults).
    """
    config = _load_packaged_defaults()

    if user_config_path:
        with open(user_config_path, "r") as f:
            if str(user_config_path).endswith((".yaml", ".yml")):
                user_config = yaml.safe_load(f) or {}
            else:
                user_config = json.load(f)
        config = _deep_merge(config, _normalize_legacy_keys(user_config))

    if overrides:
        config = _deep_merge(config, _normalize_legacy_keys(overrides))

    return config
