import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'examples/GSE68719/scripts'
spec = importlib.util.spec_from_file_location('portable_config', SCRIPTS / 'workflow_config.py')
workflow = importlib.util.module_from_spec(spec)
spec.loader.exec_module(workflow)


class PortableWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name).resolve()
        self.checkout = self.base / 'renamed-checkout'
        self.scripts = self.checkout / 'examples/GSE68719/scripts'
        shutil.copytree(SCRIPTS, self.scripts, ignore=shutil.ignore_patterns('old', 'logs', '__pycache__'))
        self.config = self.base / 'settings' / 'chosen.yaml'
        self.config.parent.mkdir()
        (self.config.parent / 'refs').mkdir()
        (self.config.parent / 'refs/A.gmt').write_text('T\tna\tA\n')
        (self.config.parent / 'input.tsv').write_text('symbol\tscore\nA\t1\n')
        self.calls = self.base / 'calls.jsonl'
        self.params = {
            'repo_results_dir': 'results',
            'default_params': {'data_csv_path': 'input.tsv', 'permutation_num': 1000,
                'gene_term_collections_gmt_file_directory_list': {'human': 'refs'},
                'gene_term_collections_2_use': ['A', 'A']},
            'make_gseapy_tables_run_control': {'run_defined_gseapy_runs': True},
            'defined_gseapy_runs': {'one': {'run': True, 'gseapy_run_name': 'score', 'DA_run_name': 'comparison'}},
            'gseapy_dotplot_params': {'run': False, 'gseapy_run_dir_list': ['absent']},
        }

    def launch(self, script='make_gseapy_tables.py'):
        self.config.write_text(yaml.safe_dump({'GSEApy_params': self.params}))
        # Modules are substituted before script execution: no network or GSEA.
        code = '''import sys,types,runpy,json
from pathlib import Path
calls=Path(sys.argv[3])
def capture(*args,**kwargs):
 with calls.open('a') as f: f.write(json.dumps({'args':args,'kwargs':kwargs},default=str)+'\\n')
for name in ('_gseapy_pre_rank_wrap','_gseapy_dot_plots'):
 module=types.ModuleType(name)
 module.run_gseapy_prerank_multiple_term_collections=capture
 module.dotplot_from_table=capture
 sys.modules[name]=module
script=sys.argv[1];config=sys.argv[2]
sys.path.insert(0,str(Path(script).parent))
sys.argv=[script,'--config',config]
runpy.run_path(script,run_name='__main__')
'''
        return subprocess.run([sys.executable, '-c', code, str(self.scripts / script), str(self.config), str(self.calls)],
                              cwd=self.base, text=True, capture_output=True,
                              env={**os.environ, 'MPLBACKEND': 'Agg', 'PYTHONDONTWRITEBYTECODE': '1'})

    def test_renamed_checkout_arbitrary_cwd_overrides_and_disabled_run(self):
        run = self.params['defined_gseapy_runs']['one']
        run.update(gsea_permutation_num=37, min_size=2, results_dir='override', seed=19)
        self.params['defined_gseapy_runs']['disabled'] = {'run': False, 'data_csv_path': 'missing'}
        result = self.launch()
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = [json.loads(line)['kwargs'] for line in self.calls.read_text().splitlines()]
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]['permutation_num'], 37)
        self.assertEqual(calls[0]['seed'], 19)
        self.assertEqual(calls[0]['min_size'], 2)
        self.assertEqual(calls[0]['list_of_term_collection_names'], ['A'])
        self.assertEqual(calls[0]['data_csv_path'], str(self.config.parent / 'input.tsv'))
        self.assertEqual(calls[0]['results_dir'], str(self.config.parent / 'override'))
        self.assertTrue(list((self.config.parent / 'override/logs').glob('*.log')))
        self.assertTrue(list((self.scripts / 'logs').glob('*.log')))

    def test_missing_later_run_prevents_any_analysis_or_output(self):
        self.params['defined_gseapy_runs']['two'] = {
            'gene_term_collections_2_use': ['absent'], 'gseapy_run_name': 'second', 'DA_run_name': 'comparison'}
        result = self.launch()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing=['absent']", result.stderr)
        self.assertFalse(self.calls.exists())
        self.assertFalse((self.config.parent / 'results').exists())

    def test_explicit_gmt_paths_and_name_overrides(self):
        run = self.params['defined_gseapy_runs']['one']
        run.update(list_of_term_collections=['refs/A.gmt'], list_of_term_collection_names=['custom'])
        self.params['default_params']['gene_term_collections_2_use'] = ['absent']
        result = self.launch()
        self.assertEqual(result.returncode, 0, result.stderr)
        call = json.loads(self.calls.read_text())['kwargs']
        self.assertEqual(call['list_of_term_collections'], [str(self.config.parent / 'refs/A.gmt')])
        self.assertEqual(call['list_of_term_collection_names'], ['custom'])
        self.calls.unlink()
        run['list_of_term_collection_names'] = []
        self.assertNotEqual(self.launch().returncode, 0)
        self.assertFalse(self.calls.exists())

    def test_ambiguous_collection_rejected(self):
        other = self.config.parent / 'other'; other.mkdir()
        (other / 'A.gmt').write_text('T\tna\tB\n')
        self.params['default_params']['gene_term_collections_gmt_file_directory_list']['other'] = 'other'
        result = self.launch()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('ambiguous=', result.stderr)
        self.assertFalse(self.calls.exists())

    def test_disabled_and_unsupported_modes(self):
        self.params['make_gseapy_tables_run_control']['run_defined_gseapy_runs'] = False
        self.assertEqual(self.launch().returncode, 0)
        self.assertFalse(self.calls.exists())
        self.assertFalse((self.config.parent / 'results').exists())
        for mode in ('run_4_all_nf_differentialabundance_runs', 'run_4_selected_nf_differentialabundance_runs'):
            self.params['make_gseapy_tables_run_control'][mode] = True
            result = self.launch()
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('unsupported', result.stderr)
            self.params['make_gseapy_tables_run_control'][mode] = False

    def test_dotplot_toggle_and_config_paths(self):
        self.assertEqual(self.launch('make_gseapy_dotplots.py').returncode, 0)
        self.assertFalse((self.config.parent / 'results').exists())
        tables = self.config.parent / 'tables'; tables.mkdir()
        (tables / 'example.csv').write_text('Term,fdr\nT,.01\n')
        (tables / 'example.rnk_df.csv').write_text('feature_name,rank\nA,1\n')
        self.params['gseapy_dotplot_params'].update(run=True, gseapy_run_dir_list=['tables'], top_terms_2_plot=7)
        result = self.launch('make_gseapy_dotplots.py')
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = [json.loads(line) for line in self.calls.read_text().splitlines()]
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]['args'], [str(tables / 'example.csv')])
        self.assertEqual(calls[0]['kwargs']['top_term'], 7)
        self.assertEqual(calls[0]['kwargs']['ofname'], str(tables / 'dotplots/example.dotplot.png'))

    def test_bundled_inventory_and_effective_default(self):
        config = ROOT / 'examples/GSE68719/config/config.yaml'
        params = yaml.safe_load(config.read_text())['GSEApy_params']['default_params']
        files, names = workflow.resolve_collections(params, config.parent)
        self.assertEqual(names, list(dict.fromkeys(params['gene_term_collections_2_use'])))
        self.assertTrue(all(Path(path).is_file() for path in files))
        self.assertIn('MGI_Mammalian_Phenotype_Level_4_2024', names)
        self.assertEqual(params['permutation_num'], 1000)
        self.assertEqual({Path(path).parent.name for path in files}, {'enrichr_Hu', 'msigdb_Hu'})

    def test_downloader_import_and_mocked_main(self):
        config_dir = self.checkout / 'config'; config_dir.mkdir(parents=True)
        script = config_dir / 'download_gseapy_gmt_files.py'
        shutil.copy(ROOT / 'config/download_gseapy_gmt_files.py', script)
        code = '''import sys,types,runpy
module=types.ModuleType('_gseapy_pre_rank_wrap')
calls=[]
for name in ('make_gmts_from_enrichr_libraries','make_gmts_from_Msigdb','convert_gmt_homologs_directory'):
 setattr(module,name,lambda **kwargs:calls.append(kwargs))
for name in ('basic_set_of_enrichr_term_collections','expansion_set_of_enrichr_term_collections','all_enrichr_libraries_2025_sin_archive'):
 setattr(module,name,[])
sys.modules[module.__name__]=module
runpy.run_path(sys.argv[1],run_name='__main__')
assert len(calls)==6,calls
'''
        result = subprocess.run([sys.executable, '-c', code, str(script)], cwd=self.base, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
