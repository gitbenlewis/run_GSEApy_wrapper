import copy
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

import _gseapy_pre_rank_wrap as wrapper


class RankInputTests(unittest.TestCase):
    def setUp(self):
        self.data = pd.DataFrame({'feature_name': ['A', 'B', 'C'],
                                  'log2FoldChange': [2., -1., 0.], 'pvalue': [.01, .1, 1.]})

    def test_valid_signed_ranks_and_input_unchanged(self):
        before = self.data.copy(deep=True)
        result = wrapper.process_input_table(self.data, save_table_output=False)
        self.assertEqual(result.feature_name.tolist(), ['A', 'C', 'B'])
        np.testing.assert_array_equal(result['rank'], [2., 0., -1.])
        pd.testing.assert_frame_equal(self.data, before)

    def test_zero_and_duplicate_policy_preserved(self):
        data = pd.DataFrame({'feature_name': ['A', 'A', 'B'],
                             'log2FoldChange': [1., -1., 1.], 'pvalue': [.01, .001, 0.]})
        result = wrapper.process_input_table(data, save_table_output=False).set_index('feature_name')
        self.assertEqual(result.loc['A', 'rank'], -3.)
        self.assertEqual(result.loc['B', 'rank'], -np.log10(np.nextafter(0, 1)))

    def test_invalid_identifiers_rejected_before_output(self):
        for value in (None, np.nan, pd.NA, '', '   '):
            with self.subTest(value=value), patch.object(pd.DataFrame, 'to_csv') as write:
                data = self.data.copy()
                data.loc[0, 'feature_name'] = value
                with self.assertRaisesRegex(ValueError, 'identifiers'):
                    wrapper.process_input_table(data)
                write.assert_not_called()

    def test_invalid_metrics_and_probability_domains(self):
        for column in ('pvalue', 'log2FoldChange'):
            for value in (None, np.nan, np.inf, -np.inf, 'bad'):
                with self.subTest(column=column, value=value):
                    data = self.data.astype({column: object})
                    data.loc[0, column] = value
                    with self.assertRaisesRegex(ValueError, column):
                        wrapper.process_input_table(data, save_table_output=False)
        for value in (-.01, 1.01):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, r'\[0, 1\]'):
                wrapper.process_input_table(self.data.assign(pvalue=value), save_table_output=False)

    def test_precomputed_validation_and_numeric_strings(self):
        for value in (None, np.nan, np.inf, 'bad'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                wrapper.process_input_table(self.data.assign(score=[value, 1, 2]),
                    rank_metric_calc_flavor='use_column_key', existing_rank_metric_column_name='score', save_table_output=False)
        result = wrapper.process_input_table(self.data.assign(score=['3', '-2', '1']),
            rank_metric_calc_flavor='use_column_key', existing_rank_metric_column_name='score', save_table_output=False)
        self.assertEqual(result['rank'].tolist(), [3, 1, -2])
        with self.assertRaisesRegex(ValueError, 'Missing required'):
            wrapper.process_input_table(self.data.drop(columns='pvalue'), save_table_output=False)

    def test_split_empty_identifiers_and_flag_precedence(self):
        for flags, identifier in [({'keep_first_occurrence_for_rows_with_multiple_feature_names': True}, ';A'),
                                  ({'explode_rows_with_multiple_feature_names': True}, 'A;;B')]:
            with self.subTest(flags=flags), self.assertRaisesRegex(ValueError, 'after splitting'):
                wrapper.process_input_table(self.data.iloc[:1].assign(feature_name=identifier),
                    rows_with_multiple_feature_names_delimiter=';', save_table_output=False, **flags)
        result = wrapper.process_input_table(self.data.iloc[:1].assign(feature_name='A;B'),
            keep_first_occurrence_for_rows_with_multiple_feature_names=True,
            explode_rows_with_multiple_feature_names=True,
            rows_with_multiple_feature_names_delimiter=';', save_table_output=False)
        self.assertEqual(result.feature_name.tolist(), ['A'])


class RunContractTests(unittest.TestCase):
    def setUp(self):
        self.data = pd.DataFrame({'feature_name': ['A', 'B'], 'rank': [2., -1.]})
        self.result = SimpleNamespace(weight=1., ascending=False, gene_sets={'T': ['A', 'B']},
            results={'T': {'name': 'prerank', 'es': .5, 'nes': 1.2, 'pval': .02,
                           'fdr': .05, 'fwerp': .1, 'matched_genes': 'A;B'}})
        self.args = dict(data_df=self.data, list_of_term_collections=[{'T': ['A', 'B']}],
                         rank_metric_calc_flavor='use_column_key', existing_rank_metric_column_name='rank')

    def test_supported_options_forward_and_unknown_rejected(self):
        with patch.object(wrapper.gp, 'prerank', return_value=self.result) as run:
            wrapper.run_single_gsea_prerank(self.data, {}, weight=0, ascending=None, no_plot=True)
            self.assertEqual(run.call_args.kwargs['weight'], 0)
            self.assertIsNone(run.call_args.kwargs['ascending'])
            self.assertTrue(run.call_args.kwargs['no_plot'])
        with patch.object(wrapper.gp, 'prerank') as run:
            with self.assertRaisesRegex(TypeError, 'Unsupported'):
                wrapper.run_single_gsea_prerank(self.data, {}, wieght=0)
            run.assert_not_called()

    def test_default_call_and_positive_count_forms_preserved(self):
        for count in (100, 1, np.int64(2), 2., '2'):
            with self.subTest(count=count), patch.object(wrapper.gp, 'prerank', return_value=self.result) as run:
                wrapper.run_single_gsea_prerank(self.data, {}, permutation_num=count)
                self.assertEqual(run.call_args.kwargs['permutation_num'], count)
                self.assertNotIn('weight', run.call_args.kwargs)
                self.assertIsNone(run.call_args.kwargs['ascending'])

    def test_invalid_permutations_fail_before_processing(self):
        for count in (0, -1, .5, 1.5, np.inf, None, True, 'bad'):
            with self.subTest(count=count), patch.object(wrapper, 'process_input_table') as prepare, \
                    patch.object(wrapper.gp, 'prerank') as run:
                with self.assertRaisesRegex(ValueError, 'positive integer'):
                    wrapper.run_single_gsea_prerank(self.data, {}, permutation_num=count)
                with self.assertRaisesRegex(ValueError, 'positive integer'):
                    wrapper.run_gseapy_prerank_multiple_term_collections(**self.args, permutation_num=count)
                prepare.assert_not_called()
                run.assert_not_called()

    def test_return_modes_and_no_save(self):
        for mode, count in [('all', 2), ('single', 2), ('default', 1)]:
            first, second = copy.deepcopy(self.result), copy.deepcopy(self.result)
            args = {**self.args, 'list_of_term_collections': self.args['list_of_term_collections'] * count}
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory, \
                    patch.object(wrapper.gp, 'prerank', side_effect=[first, second]):
                actual = wrapper.run_gseapy_prerank_multiple_term_collections(**args,
                    return_all_pre_res=mode == 'all', return_single_pre_res=mode == 'single',
                    save_table_output=False, results_dir=directory)
                if mode == 'all':
                    self.assertEqual(actual, [first, second])
                elif mode == 'single':
                    self.assertIs(actual, second)
                else:
                    self.assertIsNone(actual)
                self.assertEqual(list(Path(directory).iterdir()), [])

    def test_empty_and_conflicting_return_modes_fail_early(self):
        for options in ({'list_of_term_collections': []}, {'return_single_pre_res': True, 'return_all_pre_res': True}):
            with self.subTest(options=options), patch.object(wrapper, 'process_input_table') as prepare:
                with self.assertRaises(ValueError):
                    wrapper.run_gseapy_prerank_multiple_term_collections(**{**self.args, **options})
                prepare.assert_not_called()

    def test_metadata_families_match_colliding_outputs(self):
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(wrapper.gp, 'prerank', return_value=self.result), \
                self.assertWarns(RuntimeWarning), self.assertLogs(wrapper.LOGGER) as logs:
            wrapper.run_gseapy_prerank_multiple_term_collections(**{
                **self.args, 'list_of_term_collections': self.args['list_of_term_collections'] * 2}, results_dir=directory)
            metadata = json.loads((Path(directory) / 'datasetID.GSEA_prerank_run.run_metadata.json').read_text())
            self.assertEqual(metadata['fdr_scope'], 'per_collection')
            self.assertEqual(metadata['parameters']['permutation_num'], 100)
            self.assertEqual([c['adjustment_family'] for c in metadata['collections']], ['collection_1', 'collection_2'])
            self.assertEqual(metadata['collections'][0]['weight'], 1.)
            self.assertFalse(metadata['collections'][0]['ascending'])
            for collection in metadata['collections']:
                self.assertTrue((Path(directory) / collection['result_file']).exists())
            self.assertTrue(any('reported terms' in line for line in logs.output))
            self.assertFalse(any('significant terms' in line for line in logs.output))

    def test_metabolite_nulls_filtered_without_losing_literal_nan(self):
        data = pd.DataFrame({'SUB_PATHWAY': ['P'] * 5, 'CHEM_ID': [None, ' ', 'nan', 'A', 'A']})
        with patch.object(wrapper.pd, 'read_excel', return_value=data):
            result = wrapper.metabolon_chemical_annotation_to_pathway_table('unused.xlsx', columns_to_include=['CHEM_ID'])
        self.assertEqual(result.CHEM_ID.tolist(), [['nan', 'A', 'A']])

    def test_metadata_preserves_numpy_scalar_values(self):
        self.result.weight = np.float64(.5)
        self.result.ascending = np.bool_(False)
        with tempfile.TemporaryDirectory() as directory, patch.object(wrapper.gp, 'prerank', return_value=self.result):
            wrapper.run_gseapy_prerank_multiple_term_collections(**self.args, results_dir=directory,
                results_table_min_genes_matched=np.int64(1), make_feature_names_unique=np.bool_(True),
                min_size=np.float64(1.5), max_size=np.int64(100), threads=np.int64(2), seed=np.int64(6))
            metadata = json.loads((Path(directory) / 'datasetID.GSEA_prerank_run.run_metadata.json').read_text())
        self.assertEqual(metadata['parameters']['min_size'], 1.5)
        self.assertEqual(metadata['parameters']['results_table_min_genes_matched'], 1)
        self.assertIs(metadata['parameters']['make_feature_names_unique'], True)
        self.assertEqual(metadata['collections'][0]['weight'], .5)
        self.assertIs(metadata['collections'][0]['ascending'], False)


if __name__ == '__main__':
    unittest.main()
