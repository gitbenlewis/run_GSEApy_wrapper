import json
import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

import _gseapy_pre_rank_wrap as wrapper


class RankingPolicyTests(unittest.TestCase):
    def prepare(self, names, values, **kwargs):
        return wrapper.process_input_table(
            pd.DataFrame({'feature_name': names, 'score': values}),
            rank_metric_calc_flavor='use_column_key',
            existing_rank_metric_column_name='score', save_table_output=False, **kwargs)

    def test_ties_are_lexical_without_numeric_perturbation(self):
        ranks = self.prepare(['a', 'B', 'A', 'low'], [2., 2., 2., -1.])
        self.assertEqual(ranks.feature_name.tolist(), ['A', 'B', 'a', 'low'])
        np.testing.assert_array_equal(ranks['rank'], [2., 2., 2., -1.])

    def test_duplicate_max_absolute_and_counts(self):
        with self.assertWarnsRegex(RuntimeWarning, '1 duplicate feature IDs'):
            ranks = self.prepare(['A', 'A', 'B', 'B'], [2., -3., 1., 2.])
        self.assertEqual(dict(zip(ranks.feature_name, ranks['rank'])), {'B': 2., 'A': -3.})
        metadata = ranks.attrs['ranking_metadata']
        self.assertEqual(metadata['duplicate_id_count'], 2)
        self.assertEqual(metadata['duplicate_row_count'], 4)
        self.assertEqual(metadata['opposing_sign_duplicate_id_count'], 1)

    def test_opposite_maximum_ties_fail_before_writing(self):
        with warnings.catch_warnings(), patch.object(pd.DataFrame, 'to_csv') as save:
            warnings.simplefilter('ignore', RuntimeWarning)
            with self.assertRaisesRegex(ValueError, 'equal-magnitude opposing-sign maximum'):
                wrapper.process_input_table(
                    pd.DataFrame({'feature_name': ['A', 'A'], 'score': [2., -2.]}),
                    rank_metric_calc_flavor='use_column_key', existing_rank_metric_column_name='score')
            save.assert_not_called()

    def test_submaximum_opposite_ties_are_not_ambiguous(self):
        with self.assertWarns(RuntimeWarning):
            ranks = self.prepare(['A', 'A', 'A'], [2., -2., 3.])
        self.assertEqual(ranks['rank'].tolist(), [3.])

    def test_explicit_duplicate_retention(self):
        with self.assertWarnsRegex(RuntimeWarning, 'Duplicates are retained'):
            ranks = self.prepare(['A', 'A'], [2., -2.], make_feature_names_unique=False)
        self.assertEqual(ranks['rank'].tolist(), [2., -2.])
        self.assertEqual(ranks.attrs['ranking_metadata']['ambiguous_opposite_sign_maximum_id_count'], 1)

    def test_zero_floor_count_before_expansion_and_signed_zero_normalization(self):
        data = pd.DataFrame({'feature_name': ['B;A', 'C', 'C'],
                             'log2FoldChange': [1., 1., -1.], 'pvalue': [0., 1., 1.]})
        outputs = []
        for seed in range(3):
            with self.assertWarnsRegex(RuntimeWarning, '1 input zero p-values'):
                ranks = wrapper.process_input_table(
                    data.sample(frac=1, random_state=seed), save_table_output=False,
                    explode_rows_with_multiple_feature_names=True, rows_with_multiple_feature_names_delimiter=';')
            outputs.append((ranks.to_csv(index=False), ranks.attrs))
            np.testing.assert_array_equal(ranks['rank'], [-np.log10(np.nextafter(0, 1))] * 2 + [0.])
            self.assertFalse(np.signbit(ranks['rank'].iloc[-1]))
        self.assertEqual(outputs[0], outputs[1])
        self.assertEqual(outputs[0], outputs[2])

    def test_direct_call_preserves_order_and_explicit_overrides(self):
        ranks = pd.DataFrame({'feature_name': ['B', 'A'], 'rank': [1., 2.]})
        for options in ({}, {'ascending': True}, {'ascending': False}):
            with patch.object(wrapper.gp, 'prerank', return_value=SimpleNamespace(results={})) as run:
                wrapper.run_single_gsea_prerank(ranks, {}, **options)
            self.assertIs(run.call_args.kwargs['rnk'], ranks)
            self.assertIs(run.call_args.kwargs['ascending'], options.get('ascending'))

    def test_real_signed_gsea_shuffles_match_results_and_all_saved_bytes(self):
        # Deterministic synthetic fixture, including ties, zero probabilities,
        # signed zero ranks, and a duplicate whose maximum has a clear direction.
        names = [f'G{i:03d}' for i in range(120)]
        data = pd.DataFrame({'feature_name': names + ['G050'],
                             'log2FoldChange': [1.] * 60 + [-1.] * 60 + [-1.],
                             'padj': [0.] * 2 + [.01] * 38 + [.1] * 40 + [1.] * 40 + [.5]})
        sets = {'upper': names[:30], 'mixed': names[10:110:3]}
        snapshots = []
        metrics = []
        for seed in (1, 2, 3):
            with tempfile.TemporaryDirectory() as directory, warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                results = wrapper.run_gseapy_prerank_multiple_term_collections(
                    data_df=data.sample(frac=1, random_state=seed), list_of_term_collections=[sets],
                    min_size=5, max_size=100, permutation_num=100, threads=1,
                    verbose=False, results_dir=directory, return_all_pre_res=True)
                snapshots.append({p.name: p.read_bytes() for p in Path(directory).iterdir()})
                metrics.append({term: [value[k] for k in ('es', 'nes', 'pval', 'fdr', 'fwerp')]
                                for term, value in results[0].results.items()})
        self.assertEqual(metrics[0], metrics[1])
        self.assertEqual(metrics[0], metrics[2])
        self.assertEqual(snapshots[0], snapshots[1])
        self.assertEqual(snapshots[0], snapshots[2])
        metadata = json.loads(snapshots[0]['datasetID.GSEA_prerank_run.run_metadata.json'])
        self.assertEqual(metadata['ranking']['input_zero_pvalue_count'], 2)
        self.assertEqual(metadata['ranking']['opposing_sign_duplicate_id_count'], 1)
        self.assertIsNone(metadata['collections'][0]['ascending'])


if __name__ == '__main__':
    unittest.main()
