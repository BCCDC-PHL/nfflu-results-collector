import pytest

from tests.fixture_builder import build_fixture, RUN_ID, SAMPLE_IDS


@pytest.fixture
def analysis_dir(tmp_path):
    """Path to a synthetic nf-flu output directory (see fixture_builder.py)."""
    return build_fixture(tmp_path)


@pytest.fixture
def run_id():
    return RUN_ID


@pytest.fixture
def sample_ids():
    return list(SAMPLE_IDS)
