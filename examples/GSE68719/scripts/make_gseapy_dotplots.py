#!/usr/bin/env python3
"""Render configured saved tables without relying on checkout name or cwd."""
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import _gseapy_dot_plots as rgw
from workflow_config import load_config, resolve_path, setup_logging


def main(argv=None):
    config, config_dir = load_config(argv)
    params = config['gseapy_dotplot_params']
    if not params.get('run', True):
        logging.info('Skipping dotplots: run=false')
        return
    directories = [resolve_path(item, config_dir) for item in params['gseapy_run_dir_list']]
    for directory in directories:
        if not directory.is_dir():
            raise FileNotFoundError(f'Dotplot input directory not found: {directory}')
    setup_logging(resolve_path(config['repo_results_dir'], config_dir), __file__)
    import matplotlib.pyplot as plt
    for directory in directories:
        csv_files = [path for path in sorted(directory.glob('*.csv')) if not path.name.endswith('.rnk_df.csv')]
        if not csv_files:
            logging.warning('No result CSV files found in %s', directory)
            continue
        output_dir = directory / 'dotplots'
        output_dir.mkdir(parents=True, exist_ok=True)
        for csv_path in csv_files:
            try:
                rgw.dotplot_from_table(csv_path, title=csv_path.stem.replace('.', '\n', 1),
                                      top_term=params.get('top_terms_2_plot', 20),
                                      figsize=params.get('figsize') or (6, 8),
                                      ofname=str(output_dir / f'{csv_path.stem}.dotplot.png'))
            finally:
                plt.close('all')


if __name__ == '__main__':
    main()
