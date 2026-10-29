import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import _gseapy_dot_plots as plots


class DotplotTests(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def capture(self, data, **kwargs):
        with patch.object(plots.gp, 'dotplot', return_value=None) as draw:
            plots.dotplot_from_table(data, **kwargs)
        return draw.call_args.args[0], draw.call_args.kwargs

    def test_requested_nominal_and_explicit_adjusted_precedence(self):
        data = pd.DataFrame({'Term': ['nominal', 'adjusted'], 'pval': ['.001', '.5'], 'fdr': ['.8', '.002']})
        before = data.copy(deep=True)
        nominal, options = self.capture(data, pvalue_col='pval', top_term=1)
        adjusted, _ = self.capture(data, pvalue_col='pval', column='FDR q-val', top_term=1)
        self.assertEqual(nominal.Term.tolist(), ['nominal'])
        self.assertEqual(adjusted.Term.tolist(), ['adjusted'])
        self.assertEqual(options['column'], '-log10(pval)')
        pd.testing.assert_frame_equal(data, before)

    def test_zero_selection_raw_cutoff_and_floor_cannot_reverse_rank(self):
        data = pd.DataFrame({'Term': ['zero', 'positive'], 'fdr': [0., 1e-8]})
        selected, options = self.capture(data, cutoff=0, top_term=1, n_perm_4_zero_pavlues=10)
        self.assertEqual(selected.Term.tolist(), ['zero'])
        self.assertEqual(selected.fdr.tolist(), [0.])
        self.assertEqual(selected[options['column']].tolist(), [1.])
        selected, options = self.capture(data, n_perm_4_zero_pavlues=10)
        scores = selected.set_index('Term')[options['column']]
        self.assertGreater(scores['zero'], scores['positive'])
        self.assertAlmostEqual(scores['zero'], -np.log10(5e-9))

    def test_allzero_and_none_floor(self):
        data = pd.DataFrame({'Term': ['A', 'B'], 'fdr': [0., 0.]})
        ax = plots.dotplot_from_table(data, cutoff=0)
        np.testing.assert_array_equal(ax.collections[0].get_array(), [4., 4.])
        with self.assertRaisesRegex(ValueError, 'only zeros'):
            plots.dotplot_from_table(data, n_perm_4_zero_pavlues=None)
        data.loc[1, 'fdr'] = .02
        selected, options = self.capture(data, n_perm_4_zero_pavlues=None)
        self.assertEqual(selected.loc[selected.Term == 'A', options['column']].iloc[0], 2.)

    def test_subnormal_display_and_unrepresentable_floor(self):
        tiny = np.nextafter(0., 1.)
        ax = plots.dotplot_from_table(pd.DataFrame({'Term': ['A', 'B'], 'fdr': [tiny, 1e-310]}))
        self.assertGreater(max(ax.collections[0].get_array()), 323.)
        with self.assertRaisesRegex(ValueError, 'not representable'):
            plots.dotplot_from_table(pd.DataFrame({'Term': ['A', 'B'], 'fdr': [0., tiny]}))

    def test_mixed_aliases_missing_cells_and_groupwise_selection(self):
        a = pd.DataFrame({'Term': ['A', 'B'], 'pval': ['.01', '.2'], 'fdr': [.4, .3]})
        b = pd.DataFrame({'Term': ['C', 'D'], 'NOM p-val': ['.03', None], 'FDR q-val': [.5, .4]})
        with patch.object(plots.gp, 'dotplot', return_value=None) as draw:
            plots.dotplot_from_tables([a, b], labels=['one', 'two'], pvalue_col='P-value', top_term=1)
        data = draw.call_args.args[0]
        self.assertEqual(data.Term.tolist(), ['A', 'C'])
        self.assertEqual(data.Dataset.tolist(), ['one', 'two'])
        with self.assertRaisesRegex(ValueError, 'unavailable'):
            self.capture(a[['Term', 'fdr']], pvalue_col='pval')

    def test_full_identity_custom_y_and_saved_short_ticks(self):
        data = pd.DataFrame({'Identity': ['same_prefix_A', 'same_prefix_B'], 'fdr': [.01, .02]})
        with tempfile.TemporaryDirectory() as directory:
            fig, ax = plt.subplots()
            result = plots.dotplot_from_table(data, term_label_col='Identity', term_label_char_limit=4,
                                             ordered_terms_list=['same_prefix_B', 'same_prefix_A'],
                                             ax=ax, ofname=str(Path(directory) / 'plot.png'))
            self.assertIsNone(result)
            self.assertGreater((Path(directory) / 'plot.png').stat().st_size, 0)
            self.assertEqual([t.get_text() for t in ax.get_yticklabels()], ['same', 'same'])
            self.assertEqual(len(set(ax.collections[0].get_offsets()[:, 1])), 2)
        selected, options = self.capture(data, term_label_col='Identity', term_label_char_limit=4)
        self.assertEqual(options['y'], 'Identity')
        self.assertEqual(set(selected.Identity), {'same_prefix_A', 'same_prefix_B'})

    def test_explicit_effect_column_retains_native_behavior(self):
        selected, options = self.capture(pd.DataFrame({'Term': ['A'], 'NES': [2.]}), column='NES')
        self.assertEqual(options['column'], 'NES')
        self.assertEqual(selected.NES.tolist(), [2.])


if __name__ == '__main__':
    unittest.main()
