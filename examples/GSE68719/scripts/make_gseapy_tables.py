#!/usr/bin/env python3
"""Run enabled named analyses, with paths relative to the selected config."""
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import _gseapy_pre_rank_wrap as rgw
from workflow_config import load_config, resolve_path, resolve_collections, setup_logging


def main(argv=None):
    config, config_dir = load_config(argv)
    control = config.get('make_gseapy_tables_run_control', {})
    for mode in ('run_4_all_nf_differentialabundance_runs', 'run_4_selected_nf_differentialabundance_runs'):
        if control.get(mode, False):
            raise ValueError(f'{mode} is unsupported; configure defined_gseapy_runs instead.')
    if not control.get('run_defined_gseapy_runs', False):
        logging.info('Skipping defined runs: run_defined_gseapy_runs=false')
        return
    defaults = config.get('default_params', {})
    plans = []
    for run_name, overrides in config.get('defined_gseapy_runs', {}).items():
        params = {**defaults, **overrides}
        if not params.get('run', True):
            logging.info('Skipping %s: run=false', run_name)
            continue
        # Resolve aliases per layer: a legacy run key still overrides a default.
        permutation_num = overrides.get('permutation_num', overrides.get('gsea_permutation_num',
                              defaults.get('permutation_num', defaults.get('gsea_permutation_num', 1000))))
        if 'list_of_term_collections' in params:
            explicit = params['list_of_term_collections']
            explicit = [explicit] if isinstance(explicit, str) else explicit
            if not isinstance(explicit, list) or not explicit or not all(isinstance(item, str) for item in explicit):
                raise ValueError('list_of_term_collections must be a nonempty list of GMT file paths.')
            paths = [resolve_path(item, config_dir) for item in explicit]
            missing = [str(path) for path in paths if not path.is_file()]
            if missing:
                raise FileNotFoundError(f'Explicit GMT files not found: {missing}')
            collections = list(map(str, paths))
            names = params.get('list_of_term_collection_names', [path.stem.replace('.symbols', '') for path in paths])
            names = [names] if isinstance(names, str) else names
            if not isinstance(names, list) or len(names) != len(paths) or not all(isinstance(name, str) and name.strip() for name in names):
                raise ValueError('list_of_term_collection_names must provide one nonempty name per GMT path.')
        else:
            collections, names = resolve_collections(params, config_dir)
            names = params.get('list_of_term_collection_names', names)
            names = [names] if isinstance(names, str) else names
            if not isinstance(names, list) or len(names) != len(collections) or not all(isinstance(name, str) and name.strip() for name in names):
                raise ValueError('list_of_term_collection_names must provide one nonempty name per resolved GMT file.')
        input_path = resolve_path(params['data_csv_path'], config_dir)
        if not input_path.is_file():
            raise FileNotFoundError(f'{run_name}: input file not found: {input_path}')
        output = resolve_path(params.get('results_dir', params.get('repo_results_dir', config['repo_results_dir'])), config_dir)
        if 'results_dir' not in params:
            output = output / params['DA_run_name'] / params['gseapy_run_name']
        fallback = dict(data_csv_path_sepstr=',', feature_name_column_name='var_names',
                        rank_metric_calc_flavor='signed_neg_log10_pvalue', l2fc_column_name='log2FoldChange',
                        pvalue_column_name='pvalue', existing_rank_metric_column_name=None,
                        min_size=10, max_size=3000, results_table_min_genes_matched=1,
                        threads=10, seed=6, GSEA_run_filename_tag=None)
        kwargs = {key: params.get(key, value) for key, value in fallback.items()}
        for key in ('make_feature_names_unique', 'keep_first_occurrence_for_rows_with_multiple_feature_names',
                    'explode_rows_with_multiple_feature_names', 'rows_with_multiple_feature_names_delimiter',
                    'new_rank_metric_column_name', 'verbose', 'save_table_output'):
            if key in params:
                kwargs[key] = params[key]
        kwargs.update(data_csv_path=str(input_path), results_dir=str(output),
                      file_base_name=params.get('file_base_name', input_path.stem),
                      permutation_num=permutation_num, list_of_term_collections=collections,
                      list_of_term_collection_names=names)
        plans.append((run_name, kwargs))
    # Preflight every enabled input/reference before analysis or output writes.
    for run_name, kwargs in plans:
        setup_logging(Path(kwargs['results_dir']), __file__)
        logging.info('Running %s with %d resolved collections', run_name, len(kwargs['list_of_term_collections']))
        rgw.run_gseapy_prerank_multiple_term_collections(**kwargs)


if __name__ == '__main__':
    main()
