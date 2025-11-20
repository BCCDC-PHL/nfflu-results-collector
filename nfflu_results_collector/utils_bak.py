import yaml 
import glob
import pandas as pd 
import numpy as np

def parse_provenance(provenance_path, prefix):
    with open(provenance_path, 'r') as f:
        try:
            parsed_provenance = yaml.safe_load(f)
        except yaml.YAMLError as e:
            print(e)
            exit(-1)

    final_provenance = {}
    for provenance_item in parsed_provenance:
        if 'pipeline_name' in provenance_item: 
            final_provenance[f'{prefix}_pipeline_name'] = provenance_item['pipeline_name']
        if 'pipeline_version' in provenance_item:
            final_provenance[f'{prefix}_pipeline_version'] = provenance_item['pipeline_version']
    return final_provenance

def glob_search_single(glob_expr):
    file_list = glob.glob(glob_expr)

    if len(file_list) == 0:
        print('WARN: No paths found for glob search')
        return None
    elif len(file_list) > 1:
        print('WARN: Multiple paths found for glob search')
        return None
    return file_list[0]

def compute_numeric_pass_fail(df, pass_fail_config):
    reverse_map = {0: 'PASS', 1: 'WARN', 2: 'FAIL'}

    pass_fail_series_list = []
    for column, vals in pass_fail_config.items():
        if vals['direction'] == 'ascending':
            bins = [0, vals['warn'], vals['fail'], np.inf]
            numeric_label = [0, 1, 2]
            pass_fail_series_list.append(pd.cut(df[column], bins=bins, labels=numeric_label))
        else:
            bins = [0, vals['fail'], vals['warn'], np.inf]
            numeric_label = [2, 1, 0]
            pass_fail_series_list.append(pd.cut(df[column], bins=bins, labels=numeric_label))
    
    result = pass_fail_series_list[0]

    for pf_series in pass_fail_series_list[1:]:
        result = result.combine(pf_series, max)

    return result.map(reverse_map)

def add_global_pass_fail(df):
    reverse_map = {'PASS': 0, 'WARN': 1, 'FAIL': 2}
    pass_fail = df[df.columns[df.columns.str.endswith("pass_fail")]]

    pass_fail = pass_fail.drop(pass_fail.columns[pass_fail.isna().all()],axis=1)

    pass_fail = pass_fail.apply(lambda col: col.map(reverse_map))

    df['overall_pass_fail'] = pass_fail.max(axis=1).map({y:x for x, y in reverse_map.items()})
    # df['overall_pass_fail'] = pd.cut(df['overall_pass_fail'], bins=[-0.5, 0.5, 1.5, np.inf], labels=['PASS','WARN','FAIL'])

    return df