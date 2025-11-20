import json
import os

def load_default_config():
    """Load default config."""
    config = {}

    # load defaults first
    default_path = os.path.join(os.path.dirname(__file__), 'config', 'defaults.json')

    if not os.path.exists(default_path):
        return config

    with open(default_path, 'r') as infile:
        config = json.load(infile)

    return config
