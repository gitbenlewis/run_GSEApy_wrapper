"""Shared path and logging setup for the example CLI scripts."""
import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path

import yaml


def load_config(argv=None):
    parser = argparse.ArgumentParser(description='Run a configured GSEApy example workflow.')
    parser.add_argument('--config', type=Path, default=Path(__file__).resolve().parents[1] / 'config' / 'config.yaml')
    path = parser.parse_args(argv).config.resolve()
    with path.open() as handle:
        return yaml.safe_load(handle)['GSEApy_params'], path.parent


def resolve_path(value, config_dir):
    path = Path(value).expanduser()
    return (config_dir / path).resolve() if not path.is_absolute() else path.resolve()


def setup_logging(output_dir, script):
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    paths = [output_dir / 'logs', Path(script).resolve().parent / 'logs']
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',
                        handlers=[*(logging.FileHandler(path / f'{Path(script).stem}_{stamp}.log') for path in paths),
                                  logging.StreamHandler(sys.stdout)], force=True)


def resolve_collections(params, config_dir):
    """Resolve requested exact GMT stems once each, retaining request order."""
    directories = params['gene_term_collections_gmt_file_directory_list']
    if isinstance(directories, dict):
        directories = directories.values()
    inventory = {}
    for directory in dict.fromkeys(resolve_path(item, config_dir) for item in directories):
        for path in sorted(directory.glob('*.gmt')):
            inventory.setdefault(path.stem, set()).add(path.resolve())
    names = list(dict.fromkeys(params['gene_term_collections_2_use']))
    if not names:
        raise ValueError('At least one GMT collection must be requested.')
    missing = [name for name in names if name not in inventory]
    ambiguous = {name: sorted(map(str, inventory[name])) for name in names if len(inventory.get(name, ())) > 1}
    if missing or ambiguous:
        raise ValueError(f'GMT collection resolution failed: missing={missing}; ambiguous={ambiguous}')
    return [str(next(iter(inventory[name]))) for name in names], names
