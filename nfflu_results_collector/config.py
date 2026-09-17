import json
import os

import yaml

import nfflu_results_collector.schema as schema

_PACKAGE_DEFAULTS_PATH = os.path.join(os.path.dirname(__file__), "config", "defaults.json")


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


def load_default_config():
    """The packaged defaults, with no user overrides applied."""
    return _load_packaged_defaults()


def load_config(user_config_path=None, overrides=None):
    """Build a config dict: packaged defaults, deep-merged with an optional user
    config file (YAML or JSON, selected by extension), deep-merged with an
    optional `overrides` dict (CLI flags, or the dict passed to
    `Nfflu_Results_Collector`)."""
    config = _load_packaged_defaults()

    if user_config_path:
        with open(user_config_path, "r") as f:
            if str(user_config_path).endswith((".yaml", ".yml")):
                user_config = yaml.safe_load(f) or {}
            else:
                user_config = json.load(f)
        config = _deep_merge(config, user_config)

    if overrides:
        config = _deep_merge(config, overrides)

    return config
