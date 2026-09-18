import os

import pytest

import nfflu_results_collector.layouts as layouts


def test_both_layouts_declare_the_same_outputs():
    """The archived per-sample layout must stay complete, or the eventual
    migration silently loses whichever paths were left behind."""
    assert set(layouts.STAGE_CENTRIC) == set(layouts.SAMPLE_CENTRIC)


@pytest.mark.parametrize("layout,output_name,expected", [
    ('stage',  'mapping_dir', '/out/mapping/*'),
    ('sample', 'mapping_dir', '/out/*/mapping'),
    ('stage',  'read_counts', '/out/irma/*/tables/READ_COUNTS.txt'),
    ('sample', 'read_counts', '/out/*/irma/tables/READ_COUNTS.txt'),
    ('stage',  'mixtures_csv', '/out/mixtures/*/*_mixtures.csv'),
])
def test_output_path_renders_glob_for_each_layout(layout, output_name, expected):
    assert layouts.output_path('/out', output_name, layout=layout) == expected


def test_output_path_renders_a_single_sample():
    assert layouts.output_path('/out', 'mixtures_csv', sample='S1', layout='stage') == '/out/mixtures/S1/S1_mixtures.csv'


@pytest.mark.parametrize("platform,expected", [
    ('illumina', '/out/blast/blastn/irma/S1*blastn.txt'),
    ('nanopore', '/out/blast/irma/S1*blastn.txt'),
])
def test_blast_paths_differ_by_platform(platform, expected):
    assert layouts.output_path('/out', 'blastn_ref', sample='S1', platform=platform, layout='stage') == expected


def _touch(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, 'a').close()


@pytest.mark.parametrize("layout,rel", [
    ('stage',  'irma/{s}/tables/READ_COUNTS.txt'),
    ('sample', '{s}/irma/tables/READ_COUNTS.txt'),
])
def test_find_by_sample_recovers_sample_names_under_either_layout(tmp_path, layout, rel):
    for sample in ('S1', 'S2'):
        _touch(os.path.join(str(tmp_path), rel.format(s=sample)))

    found = layouts.find_by_sample(str(tmp_path), 'read_counts', layout=layout)

    assert [sample for sample, _ in found] == ['S1', 'S2']
    assert all(os.path.exists(path) for _, path in found)


def test_find_by_sample_recovers_sample_name_when_template_repeats_it(tmp_path):
    _touch(os.path.join(str(tmp_path), 'mixtures', 'S1', 'S1_mixtures.csv'))

    assert layouts.find_by_sample(str(tmp_path), 'mixtures_csv', layout='stage')[0][0] == 'S1'


def test_find_by_sample_refuses_templates_where_the_sample_name_is_ambiguous(tmp_path):
    """'{sample}*blastn.txt' cannot be split reliably, so find_by_sample must refuse
    rather than return a truncated name."""
    with pytest.raises(ValueError):
        layouts.find_by_sample(str(tmp_path), 'blastn_ref', layout='stage')
