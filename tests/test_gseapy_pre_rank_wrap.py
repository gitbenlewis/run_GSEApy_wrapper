import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd

import _gseapy_pre_rank_wrap as wrapper


class ResultsSerializationTests(unittest.TestCase):
    def setUp(self):
        self.result = SimpleNamespace(
            gene_sets={'pathway': ['A', 'B']},
            results={'pathway': {
                'name': 'prerank', 'es': 0.5, 'nes': 1.5, 'pval': 0.02,
                'fdr': 0.04, 'fwerp': 0.1, 'matched_genes': 'A;B',
                'tag %': '1/2', 'RES': [0.5, 0.0],
            }},
        )

    def test_omitted_names_write_csv_and_string_metadata(self):
        cases = [
            ({'pathway': ['A', 'B']}, 'InMemory_TermCollection'),
            ('/references/collection.symbols.gmt', 'collection.symbols'),
            (Path('/references/collection.symbols.GMT'), 'collection.symbols'),
            ('MSigDB_Hallmark_2020', 'MSigDB_Hallmark_2020'),
            (None, 'Unknown_TermCollection'),
        ]
        for source, expected in cases:
            with self.subTest(source=source), tempfile.TemporaryDirectory() as directory:
                self.result.gene_sets = source
                frame = wrapper.make_results_dataframe(
                    self.result, results_table_min_genes_matched=0, results_dir=directory,
                )
                output = Path(directory) / f'datasetID.GSEA_prerank_run.{expected}.csv'
                saved = pd.read_csv(output)
                self.assertEqual(saved['TermCollection'].tolist(), [expected])
                self.assertEqual(frame['TermCollection_Term'].tolist(), [expected + '__pathway'])
                self.assertEqual(saved['nes'].tolist(), [1.5])
                self.assertEqual(saved['fdr'].tolist(), [0.04])
                self.assertNotIn('RES', saved.columns)

    def test_explicit_label_preserved_with_safe_filename(self):
        with tempfile.TemporaryDirectory() as directory:
            frame = wrapper.make_results_dataframe(
                self.result, collection_name='custom/path\\name',
                results_table_min_genes_matched=0, results_dir=directory,
            )
            self.assertEqual(frame['TermCollection'].tolist(), ['custom/path\\name'])
            self.assertTrue((Path(directory) / 'datasetID.GSEA_prerank_run.custom_path_name.csv').is_file())

    def test_single_run_resolves_metadata_before_analysis(self):
        with patch.object(wrapper.gp, 'prerank', return_value=self.result) as prerank:
            result = wrapper.run_single_gsea_prerank(pd.DataFrame(), self.result.gene_sets)
        self.assertEqual(result.results['pathway']['TermCollection'], 'InMemory_TermCollection')
        self.assertIs(prerank.call_args.kwargs['gene_sets'], self.result.gene_sets)

    def test_invalid_explicit_name_fails_before_analysis(self):
        with patch.object(wrapper.gp, 'prerank') as prerank:
            with self.assertRaisesRegex(ValueError, 'nonempty string'):
                wrapper.run_single_gsea_prerank(pd.DataFrame(), {}, collection_name=' ')
        prerank.assert_not_called()

    def test_schema_error_propagates(self):
        del self.result.results['pathway']['matched_genes']
        with self.assertRaises(KeyError):
            wrapper.make_results_dataframe(self.result, save_table_output=False)

    def test_write_failure_prevents_wrapper_success_logs(self):
        with patch.object(wrapper, 'process_input_table', return_value=pd.DataFrame()), \
                patch.object(wrapper.gp, 'prerank', return_value=self.result), \
                patch.object(pd.DataFrame, 'to_csv', side_effect=OSError('disk full')), \
                self.assertLogs(wrapper.LOGGER, level='INFO') as logs:
            with self.assertRaisesRegex(OSError, 'disk full'):
                wrapper.run_gseapy_prerank_multiple_term_collections(
                    data_df=pd.DataFrame(), list_of_term_collections=[self.result.gene_sets],
                )
        self.assertFalse(any('Results saved' in line or 'Processed term collection:' in line
                             or 'analysis completed' in line for line in logs.output))

    def test_batch_collisions_warn_and_save_distinct_files(self):
        sets = self.result.gene_sets
        cases = [
            ([sets, sets], None, ['InMemory_TermCollection', 'InMemory_TermCollection__2']),
            ([sets, sets], ['same', 'same'], ['same', 'same__2']),
            (['/first/library.gmt', '/second/library.gmt'], None, ['library', 'library__2']),
            ([sets, sets], ['A/B', 'A_B'], ['A_B', 'A_B__2']),
            ([sets, sets], ['Case', 'case'], ['Case', 'case__2']),
            ([sets], ['rnk_df'], ['rnk_df__2']),
            ([sets, sets, sets], ['first', 'later', 'later'], ['first', 'later', 'later__2']),
            ([sets] * 4, ['A', 'A', 'A__2', 'A'], ['A', 'A__3', 'A__2', 'A__4']),
        ]
        data = pd.DataFrame({'feature_name': ['A', 'B'], 'score': [2.0, -1.0]})
        for collections, names, components in cases:
            with self.subTest(names=names, components=components), \
                    tempfile.TemporaryDirectory() as directory, \
                    patch.object(wrapper.gp, 'prerank', return_value=self.result) as prerank, \
                    self.assertLogs(wrapper.LOGGER, level='INFO'), \
                    warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                for repeat in range(2):
                    wrapper.run_gseapy_prerank_multiple_term_collections(
                        data_df=data, list_of_term_collections=collections,
                        list_of_term_collection_names=names, results_dir=directory,
                        rank_metric_calc_flavor='use_column_key', existing_rank_metric_column_name='score',
                    )
                    current = {p.name: p.read_bytes() for p in Path(directory).iterdir()}
                    if repeat == 0:
                        first_outputs = current
                    else:
                        self.assertEqual(current, first_outputs)
                prefix = 'datasetID.GSEA_prerank_run.'
                self.assertEqual(set(current), {prefix + name + '.csv' for name in components + ['rnk_df']} | {prefix + 'run_metadata.json'})
                self.assertEqual(prerank.call_count, len(collections) * 2)
                expected_labels = names or [
                    'InMemory_TermCollection' if isinstance(source, dict) else 'library'
                    for source in collections
                ]
                for component, label in zip(components, expected_labels):
                    saved = pd.read_csv(Path(directory) / (prefix + component + '.csv'))
                    self.assertEqual(saved['TermCollection'].tolist(), [label])
                    self.assertEqual(saved['TermCollection_Term'].tolist(), [label + '__pathway'])
                    self.assertNotIn('RES', saved.columns)
                    for metric in ('es', 'nes', 'pval', 'fdr', 'fwerp'):
                        self.assertEqual(saved[metric].tolist(), [self.result.results['pathway'][metric]])
                self.assertEqual(pd.read_csv(Path(directory) / (prefix + 'rnk_df.csv'))['rank'].tolist(), [2.0, -1.0])
                self.assertTrue(caught)
                for warning in caught:
                    self.assertIs(warning.category, RuntimeWarning)
                    self.assertIn('saving as', str(warning.message))
                    self.assertIn('Collection metadata is unchanged', str(warning.message))

    def test_batch_invalid_names_fail_before_processing(self):
        for names in (['only_one'], ['valid', ' ']):
            with self.subTest(names=names), \
                    patch.object(wrapper, 'process_input_table') as prepare, \
                    patch.object(wrapper.gp, 'prerank') as prerank:
                with self.assertRaises(ValueError):
                    wrapper.run_gseapy_prerank_multiple_term_collections(
                        data_df=pd.DataFrame(),
                        list_of_term_collections=[self.result.gene_sets] * 2,
                        list_of_term_collection_names=names,
                    )
                prepare.assert_not_called()
                prerank.assert_not_called()

    def test_batch_valid_names_preserve_files_labels_order_and_values(self):
        data = pd.DataFrame({'feature_name': ['A', 'B'], 'score': [2.0, -1.0]})
        cases = [(None, ['InMemory_TermCollection']), (['Second', 'First'], ['Second', 'First'])]
        for names, expected_names in cases:
            with self.subTest(names=names), tempfile.TemporaryDirectory() as directory, \
                    patch.object(wrapper.gp, 'prerank', return_value=self.result) as prerank, \
                    patch.object(wrapper, 'make_results_dataframe', wraps=wrapper.make_results_dataframe) as serialize:
                wrapper.run_gseapy_prerank_multiple_term_collections(
                    data_df=data, list_of_term_collections=[self.result.gene_sets] * len(expected_names),
                    list_of_term_collection_names=names, results_dir=directory,
                    rank_metric_calc_flavor='use_column_key', existing_rank_metric_column_name='score',
                )
                self.assertEqual([call.args[1] for call in serialize.call_args_list], expected_names)
                self.assertEqual(prerank.call_count, len(expected_names))
                prefix = 'datasetID.GSEA_prerank_run.'
                self.assertEqual(
                    {path.name for path in Path(directory).iterdir()},
                    {prefix + name + '.csv' for name in expected_names + ['rnk_df']} | {prefix + 'run_metadata.json'},
                )
                saved_ranks = pd.read_csv(Path(directory) / (prefix + 'rnk_df.csv'))
                self.assertEqual(saved_ranks.to_dict('list'), {'feature_name': ['A', 'B'], 'rank': [2.0, -1.0]})
                for name in expected_names:
                    saved = pd.read_csv(Path(directory) / (prefix + name + '.csv'))
                    self.assertEqual(saved['TermCollection'].tolist(), [name])
                    for metric in ('es', 'nes', 'pval', 'fdr', 'fwerp'):
                        self.assertEqual(saved[metric].tolist(), [self.result.results['pathway'][metric]])

    def test_batch_no_save_permits_duplicate_names(self):
        for names in (None, ['same', 'same']):
            with self.subTest(names=names), \
                    patch.object(wrapper, 'process_input_table', return_value=pd.DataFrame()), \
                    patch.object(wrapper.gp, 'prerank', return_value=self.result) as prerank, \
                    patch.object(pd.DataFrame, 'to_csv') as write, \
                    warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                wrapper.run_gseapy_prerank_multiple_term_collections(
                    data_df=pd.DataFrame(), list_of_term_collections=[self.result.gene_sets] * 2,
                    list_of_term_collection_names=names, save_table_output=False,
                )
                self.assertEqual(prerank.call_count, 2)
                write.assert_not_called()
                self.assertEqual(caught, [])

    def test_no_save_retains_results_without_writing(self):
        with patch.object(pd.DataFrame, 'to_csv') as write:
            frame = wrapper.make_results_dataframe(
                self.result, collection_name='Explicit',
                results_table_min_genes_matched=0, save_table_output=False,
            )
        write.assert_not_called()
        self.assertEqual(frame['TermCollection'].tolist(), ['Explicit'])
        self.assertIn('RES', frame.columns)


if __name__ == '__main__':
    unittest.main()
