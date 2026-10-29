import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.legend import Legend
import numpy as np
import pandas as pd

import _gseapy_dot_plots as plots


class LollipopTests(unittest.TestCase):
    def setUp(self):
        self.table = pd.DataFrame({
            "Term": list("ABCDE"), "Label": ["Shared label prefix one", "Shared label prefix two", "C", "D", "E"],
            "nes": [2.5, -1.5, 0, 1.2, -2.0], "fdr": [0, 0.01, 0.05, 0.2, np.nan],
            "fraction": [1, 0.5, 0, 0.25, 0.75],
        })

    def tearDown(self):
        plt.close("all")

    def test_signed_stems_identity_order_and_shortened_labels(self):
        ax = plots.lollipop_from_table(self.table, term_label_col="Label", term_label_char_limit=6,
                                      ordered_terms_list=["B", "A", "C"])
        segments = ax.collections[0].get_segments()
        np.testing.assert_allclose(segments, [[[0, 0], [-1.5, 0]], [[0, 1], [2.5, 1]], [[0, 2], [0, 2]]])
        np.testing.assert_allclose(ax.collections[1].get_offsets(), [[-1.5, 0], [2.5, 1], [0, 2]])
        np.testing.assert_allclose(ax.collections[1].get_sizes(), [64, 64, 64])
        self.assertEqual([tick.get_text() for tick in ax.get_yticklabels()], ["Shared", "Shared", "C"])
        self.assertIsNone(ax.get_legend())

    def test_missing_terms_warn_without_substitution(self):
        with self.assertWarnsRegex(UserWarning, "missing"):
            ax = plots.lollipop_from_table(self.table, ordered_terms_list=["B", "missing", "A"])
        self.assertEqual([tick.get_text() for tick in ax.get_yticklabels()], ["B", "A"])
        with self.assertWarns(UserWarning), self.assertRaisesRegex(ValueError, "No terms"):
            plots.lollipop_from_table(self.table, ordered_terms_list=["missing"])

    def test_identity_validation(self):
        for values in (["A"] * 5, [None, "B", "C", "D", "E"]):
            with self.subTest(values=values), self.assertRaises(ValueError):
                plots.lollipop_from_table(self.table.assign(Term=values))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            plots.lollipop_from_table(self.table, ordered_terms_list=["A", "A"])

    def test_probability_bars_floor_missing_and_alignment(self):
        colors = {"A": "red", "B": "blue", "C": "green", "D": "orange", "E": "black"}
        left, right = plots.lollipop_from_table(self.table, significance_panel=True, zero_pvalue_floor=1e-4,
                                              term_colors=colors, text_kwargs={"significance": "numeric"})
        np.testing.assert_allclose([bar.get_width() for bar in right.patches], [4, 2, -np.log10(.05), -np.log10(.2)])
        np.testing.assert_allclose([bar.get_y() + bar.get_height()/2 for bar in right.patches], [0, 1, 2, 3])
        np.testing.assert_allclose(left.get_ylim(), right.get_ylim())
        np.testing.assert_allclose(left.collections[1].get_facecolors()[:4], [bar.get_facecolor() for bar in right.patches])
        self.assertEqual([text.get_text() for text in right.texts], ["0", "0.01", "0.05", "0.2", "NA"])
        self.assertIn("display only", right.get_xlabel())
        self.assertIn("Adjusted P-value", right.get_xlabel())

    def test_floor_and_probability_validation(self):
        for floor in (None, 0, -1, 1, np.inf, 0.02):
            with self.subTest(floor=floor), self.assertRaises(ValueError):
                plots.lollipop_from_table(self.table, significance_panel=True, zero_pvalue_floor=floor)
        for invalid in (-0.1, 1.1, np.inf, "bad"):
            table = self.table.copy()
            table["fdr"] = [invalid, .01, .05, .2, np.nan]
            with self.subTest(invalid=invalid), self.assertRaises((ValueError, TypeError)):
                plots.lollipop_from_table(table, text_kwargs={"significance": "numeric"})

    def test_thresholds_inclusive_and_annotations_without_panel(self):
        ax = plots.lollipop_from_table(self.table, significance_thresholds={.05: "*", .01: "**"},
                                      text_kwargs={"significance": "symbols"})
        self.assertEqual([text.get_text() for text in ax.texts], ["**", "**", "*", "ns", "NA"])
        self.assertIn("Adjusted P-value", ax.get_xlabel())
        labels = [text.get_text() for text in ax.get_legend().get_texts()]
        self.assertEqual(labels, ["**: Adjusted P-value ≤ 0.01", "*: Adjusted P-value ≤ 0.05"])
        with self.assertRaisesRegex(ValueError, "require significance_thresholds"):
            plots.lollipop_from_table(self.table, text_kwargs={"significance": "symbols"})

    def test_custom_columns_and_no_nominal_fallback(self):
        table = self.table.rename(columns={"Term": "id", "nes": "effect", "fdr": "q"})
        ax = plots.lollipop_from_table(table, term_id_col="id", term_label_col="Label", nes_col="effect", pvalue_col="q",
                                      text_kwargs={"significance": "numeric"})
        self.assertIn("Annotations: q", ax.get_xlabel())
        with self.assertRaisesRegex(ValueError, "No nominal-P fallback"):
            plots.lollipop_from_table(self.table.drop(columns="fdr").assign(pval=.01), significance_panel=True)

    def test_fraction_area_reference_and_size_legend(self):
        ax = plots.lollipop_from_table(self.table, fraction_col="fraction",
                                      marker_kwargs={"fraction_scale": 400, "reference_fraction": 1})
        np.testing.assert_allclose(ax.collections[1].get_sizes(), [400])
        np.testing.assert_allclose(ax.collections[2].get_sizes(), [400, 200, 0, 100, 300])
        handles = ax.get_legend().legend_handles
        np.testing.assert_allclose([handle.get_markersize()**2 for handle in handles], [100, 200, 400])
        with self.assertRaisesRegex(ValueError, "fraction_scale"):
            plots.lollipop_from_table(self.table, fraction_col="fraction")

    def test_legends_independent_visibility_order_and_composition(self):
        fig, (ax, unrelated) = plt.subplots(1, 2)
        unrelated.plot([1, 2], [3, 4])
        result = plots.lollipop_from_table(
            self.table, ax=ax, fraction_col="fraction", marker_kwargs={"fraction_scale": 200},
            significance_thresholds={.05: "*", .01: "**"}, text_kwargs={"significance": "symbols"},
            legend_kwargs={"size": {"visible": False}, "significance": {"order": [.05, .01], "ncol": 2, "fontsize": 7},
                           "color": {"visible": True, "order": ["B", "A"], "title": "Selected", "marker_size": 5}},
        )
        self.assertIs(result, ax)
        self.assertEqual(len(fig.axes), 2)
        self.assertEqual(len(unrelated.lines), 1)
        legends = [child for child in ax.get_children() if isinstance(child, Legend)]
        self.assertEqual(len(legends), 2)
        self.assertEqual([text.get_text() for text in legends[0].get_texts()], ["*: Adjusted P-value ≤ 0.05", "**: Adjusted P-value ≤ 0.01"])
        self.assertEqual([text.get_text() for text in legends[1].get_texts()], ["B", "A"])

    def test_large_marker_labels_and_legend_do_not_overlap(self):
        ax = plots.lollipop_from_table(
            self.table, fraction_col="fraction", figsize=(10, 7),
            marker_kwargs={"fraction_scale": 1600, "reference_fraction": 1},
            text_kwargs={"nes_labels": True, "offset": 7},
        )
        ax.figure.canvas.draw()
        renderer = ax.figure.canvas.get_renderer()
        radius_pixels = np.sqrt(1600) / 2 * ax.figure.dpi / 72
        for row, (value, annotation) in enumerate(zip(self.table.nes, ax.texts)):
            center_x = ax.transData.transform((value, row))[0]
            bounds = annotation.get_window_extent(renderer)
            if value >= 0:
                self.assertGreater(bounds.x0, center_x + radius_pixels)
            else:
                self.assertLess(bounds.x1, center_x - radius_pixels)
        bounds = [handle.get_window_extent(renderer) for handle in ax.get_legend().legend_handles]
        for upper, lower in zip(bounds, bounds[1:]):
            self.assertGreater(upper.y0, lower.y1)

    def test_input_and_settings_unchanged(self):
        before = self.table.copy(deep=True)
        settings = {"marker_kwargs": {"fraction_scale": 100}, "text_kwargs": {"significance": "symbols", "nes_labels": True},
                    "legend_kwargs": {"size": {"order": [1., .5]}}}
        previous = copy.deepcopy(settings)
        plots.lollipop_from_table(self.table, fraction_col="fraction", significance_thresholds={.05: "*"}, **settings)
        pd.testing.assert_frame_equal(self.table, before)
        self.assertEqual(settings, previous)

    def test_supplied_paired_axes_and_exports(self):
        with tempfile.TemporaryDirectory() as directory:
            csv_path = Path(directory) / "table.csv"
            self.table.to_csv(csv_path, index=False)
            for extension in ("png", "pdf", "svg"):
                fig, axes = plt.subplots(1, 2)
                path = Path(directory) / f"plot.{extension}"
                returned = plots.lollipop_from_table(csv_path, ax=axes, significance_panel=True, zero_pvalue_floor=1e-4, ofname=path)
                self.assertIs(returned[0], axes[0])
                self.assertIs(returned[1], axes[1])
                self.assertGreater(path.stat().st_size, 100)

    def test_long_labels_and_annotations_fit_figure(self):
        ax = plots.lollipop_from_table(self.table, term_label_col="Label", text_kwargs={"nes_labels": True, "wrap_width": 15},
                                      figsize=(8, 5))
        ax.figure.canvas.draw()
        renderer = ax.figure.canvas.get_renderer()
        for artist in ax.get_yticklabels() + ax.texts:
            bounds = artist.get_window_extent(renderer)
            self.assertGreaterEqual(bounds.x0, 0)
            self.assertLessEqual(bounds.x1, ax.figure.bbox.x1)
            self.assertGreaterEqual(bounds.y0, 0)
            self.assertLessEqual(bounds.y1, ax.figure.bbox.y1)

    def test_legacy_dotplot_forwarding_unchanged(self):
        table = self.table.assign(fdr=.05)
        with patch.object(plots.gp, "dotplot", return_value="sentinel") as draw:
            result = plots.dotplot_from_table(table, ordered_terms_list=["B", "A"], figsize=(3, 4))
        self.assertEqual(result, "sentinel")
        self.assertEqual(draw.call_args.kwargs["y_order"], ["B", "A"])
        self.assertEqual(draw.call_args.kwargs["column"], "-log10(fdr)")
        self.assertEqual(draw.call_args.kwargs["top_term"], 2)
        self.assertEqual(draw.call_args.kwargs["figsize"], (3, 4))


if __name__ == "__main__":
    unittest.main()
