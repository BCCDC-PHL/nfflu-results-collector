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
def test_blast_paths_differ_by_platform_under_the_stage_layout(platform, expected):
    assert layouts.output_path('/out', 'blastn_ref', sample='S1', platform=platform, layout='stage') == expected


@pytest.mark.parametrize("platform", ['illumina', 'nanopore'])
def test_blast_paths_are_the_same_on_both_platforms_under_the_sample_layout(platform):
    """The per-sample nextflow config is one file covering both platforms and
    publishes BLAST to a single path, unlike upstream's two per-platform configs."""
    path = layouts.output_path('/out', 'blastn_ref', sample='S1', platform=platform, layout='sample')
    assert path == '/out/S1/blast/irma/S1*blastn.txt'


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


@pytest.mark.parametrize("layout,rel", [
    ('stage',  'mapping/S1'),
    ('sample', 'S1/mapping'),
])
def test_detect_layout_reads_the_layout_off_the_directory(tmp_path, layout, rel):
    os.makedirs(os.path.join(str(tmp_path), rel))
    assert layouts.detect_layout(str(tmp_path)) == layout


def test_detect_layout_returns_none_when_nothing_matches(tmp_path):
    assert layouts.detect_layout(str(tmp_path)) is None


@pytest.mark.parametrize("analysis_type", ['short', 'long'])
@pytest.mark.parametrize("trailing", ['', os.sep])
def test_detect_analysis_type_reads_the_run_tree(tmp_path, analysis_type, trailing):
    outdir = os.path.join(str(tmp_path), 'RUN-1', analysis_type, 'nf-flu-3.10-output')
    os.makedirs(outdir)
    assert layouts.detect_analysis_type(outdir + trailing) == analysis_type


def test_detect_analysis_type_returns_none_outside_the_run_tree(tmp_path):
    """A standalone `-d /some/path` run isn't under analysis_output/<run>/<type>,
    so the column stays empty rather than filling with a stray directory name."""
    outdir = os.path.join(str(tmp_path), 'somewhere', 'nf-flu-output')
    os.makedirs(outdir)
    assert layouts.detect_analysis_type(outdir) is None
