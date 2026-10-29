# GSEA lollipop plots

`lollipop_from_table` plots existing enrichment results; it never reruns enrichment,
recomputes adjusted probabilities, or modifies the input table. Positive NES values
extend right from zero and negative NES values extend left. Interpret direction
relative to the ranking and contrast used to generate your results, not as a new
test of activation or inhibition.

## Complete examples

Run from the parent directory of this repository, or add that directory to
`PYTHONPATH`, using the existing environment with GSEApy and Matplotlib:

```python
import pandas as pd
from run_GSEApy_wrapper import lollipop_from_table

results = pd.DataFrame({
    "Term": list("ABCDE"),
    "Label": ["Immune response", "Energy metabolism", "Protein transport",
              "Cell cycle", "Extracellular matrix"],
    "nes": [2.4, -1.8, 1.3, -0.8, 0.45],
    "fdr": [0, 0.01, 0.05, 0.2, None],
    "fraction": [0.8, 0.6, 0.35, 0.2, 0.45],
})
colors = {"A": "#4477aa", "B": "#ee6677", "C": "#228833",
          "D": "#aa3377", "E": "#66ccee"}
common = dict(
    term_label_col="Label", ordered_terms_list=list("ABCDE"), term_colors=colors,
    significance_thresholds={0.001: "***", 0.01: "**", 0.05: "*"},
    text_kwargs={"nes_labels": True, "significance": "symbols"},
    legend_kwargs={"significance": {"loc": "upper left", "bbox_to_anchor": (0, -0.22),
                                    "ncol": 3, "fontsize": 8, "handlelength": 0}},
)

# Equal-size circles; symbols use the supplied probability, including zero.
ax = lollipop_from_table(results, **common, ofname="equal_size.png")

# The rendering floor applies only to bars for reported zeros.
nes_ax, probability_ax = lollipop_from_table(
    results, **common, significance_panel=True, zero_pvalue_floor=1e-4,
    ofname="significance_panel.pdf",
)

# Circle AREA in points squared is fraction * 500.
axes = lollipop_from_table(
    results, **common, fraction_col="fraction",
    marker_kwargs={"fraction_scale": 500, "reference_fraction": 1},
    significance_panel=True, zero_pvalue_floor=1e-4, ofname="fraction_sizes.svg",
)
```

The checked-in [render script](../examples/lollipop/render_examples.py) uses the
same public API and a synthetic A–E table with long labels. It can run from any
working directory:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/gsea_lollipop_mpl python /absolute/path/to/run_GSEApy_wrapper/examples/lollipop/render_examples.py
```

It creates three layouts, each in PNG/PDF/SVG, under
`examples/lollipop/rendered/`. No input download or enrichment run is needed.

## Input and selection

`table_path` accepts a DataFrame or a CSV path; `sep` defaults to `","`.
`term_id_col="Term"` selects identity, `term_label_col=None` uses that identity
for labels, and `ordered_terms_list=None` preserves input row order. An explicit
order selects those IDs in that order. Missing requested IDs produce a warning;
no substitute rows are selected. An all-missing selection raises an error.

IDs are compared as strings and must be unique, nonblank and nonmissing. Use
`term_id_col="TermCollection_Term"` when a combined table repeats term names
across collections. Duplicate requested IDs raise an error. Distinct IDs can
have identical labels: numeric row coordinates keep them separate.
`term_label_char_limit` truncates display text only. `text_kwargs['wrap_width']`
wraps labels (default 38 characters; `None` disables wrapping).

`nes_col="NES"` and `pvalue_col="Adjusted P-value"` use the existing aliases:
lowercase `nes` supplies `NES`, and `fdr` or `FDR q-val` supplies
`Adjusted P-value` when absent. Existing explicit columns take precedence;
columns are never coalesced row by row. For a specific adjusted statistic,
pass its actual column name. Nominal P values are never selected as a fallback.
Probability columns are required only when bars or probability annotations
are requested. Nonfinite NES values raise an error; missing probabilities are
allowed, but nonnumeric, infinite, negative and above-one values raise.

## Scientific display rules

The bar width is `-log10(supplied adjusted probability)`. Missing probabilities
have no bar and show `NA`, never zero or a significance symbol. A probability of
one has zero bar width and remains a valid measured value.

For reported zeros, bar mode requires `zero_pvalue_floor` strictly between zero
and one, and no greater than the smallest positive selected probability. Only
zeros are replaced in the private display array. This prevents a zero from
looking less significant than a positive value. The bar-axis label identifies
the selected column and rendering floor. The floor is **not** a measured value,
an estimate, a permutation bound, or a new significance threshold. Numeric
annotations display the original reported `0`. All-zero data can be drawn with
an explicit floor; entirely missing data shows `NA` without invented bars.

Symbols require `significance_thresholds={upper_probability: symbol}`. Thresholds
are sorted numerically; the smallest threshold satisfying `p <= threshold` wins.
For `{0.01: "**", 0.05: "*"}`, exactly 0.01 is `**`, exactly 0.05 is `*`, above
0.05 is `ns`, and missing is `NA`. Symbols use original probabilities, not the
bar floor. No significance threshold is imposed by default. The legend is
derived from the same mapping; annotations-only plots label their probability
column beneath the NES axis.

## Styling dictionaries

| Dictionary | Supported semantic keys and defaults |
| --- | --- |
| `marker_kwargs` | `size=64` constant area in points²; `color="#4477aa"` fallback color; `fraction_scale=None` required with `fraction_col`; `reference_fraction=None`; `legend_fractions=[0.25, 0.5, 1.0]`. Other options go to Matplotlib `scatter`, e.g. `edgecolors`, `linewidths`, `alpha`, `marker`, `zorder`. Use `size`, not `s`; use `term_colors`, not `c`. |
| `stem_kwargs` | Passed to `hlines`; default `linewidths=1.8`. Supports e.g. `linestyles`, `alpha`, `zorder`. Row colors come from `term_colors`. |
| `text_kwargs` | `nes_labels=False`; `significance="none"` (`"numeric"` or `"symbols"`); `nes_format=".2f"`; `pvalue_format=".3g"`; `offset=7` points beyond the marker edge (or bar endpoint); `wrap_width=38`. Remaining text properties such as `fontsize=9`, `color`, and `fontweight` style annotations. Alignment/coordinates are owned by the renderer. |
| `panel_kwargs` | `title=None`; `nes_xlabel="Normalized enrichment score (NES)"`; `significance_xlabel` defaults to the selected probability column's transform; `bar_height=0.55`; `width_ratios=(1.2, 1)`; `wspace=0.12`; `xmargin` defaults to 0.25 with annotations, otherwise 0.08; `zero_line_kwargs` passes style to the vertical zero reference. |
| `legend_kwargs` | Independent dictionaries named `size`, `significance`, `color`; see below. |

Fractions must be finite and in `[0,1]`. Area is exactly `fraction * fraction_scale`,
so a zero fraction has zero marker area. No minimum size or data-dependent
rescaling is added. `reference_fraction` draws an outline reference circle behind
each point using the same area scale. It requires fraction mode. Size legends
use the same scale and appear only in fraction mode, even if requested in
constant-size mode. Caller `term_colors` maps stable term IDs to colors; unspecified
IDs use the fallback marker color. Both panels use the same row colors.

Each legend supports `visible`, `order`, `loc`, `bbox_to_anchor`, `ncol`, `title`,
`fontsize`, `title_fontsize`, `marker_size`, and all standard Matplotlib legend
spacing options (`labelspacing`, `columnspacing`, `handlelength`, `handletextpad`,
`borderpad`, etc.). `order` is an ordered subset of fraction values for `size`,
threshold numbers for `significance`, or string term IDs for `color`. It may also
omit entries. Unknown or duplicate order keys raise an error. `marker_size`
overrides legend marker diameter in points; leave it unset in a size legend to
preserve the plotted area scale. `markerscale` likewise changes displayed legend
marker size only. Significance legend entries are text, so marker sizing has no
visual effect there.

Size/significance legends default visible when their encoding is active; term
color legends default hidden. Each has its own anchor and survives the creation
of the others. Use explicit anchors for dense compositions, as in the example.
For example, hide only size and reorder only significance:

```python
legend_kwargs={
    "size": {"visible": False},
    "significance": {"order": [0.05, 0.01], "ncol": 2, "loc": "lower left"},
    "color": {"visible": True, "order": ["B", "A"], "title": "Terms"},
}
```

## Figure composition and saving

Single-panel mode accepts one `ax` and returns that Axes. Two-panel mode accepts
two distinct Axes from the same figure and returns `(nes_ax, probability_ax)`.
The helper does not close supplied figures or create another figure when axes
are provided. `figsize` affects only internally created figures; internal figures
use constrained layout and a row-aware default height. Caller-owned layouts remain
the caller's responsibility.

```python
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 2, figsize=(10, 5), layout="constrained")
lollipop_from_table(results, ax=axes, significance_panel=True, zero_pvalue_floor=1e-4)
fig.savefig("composed.svg", bbox_inches="tight")
```

`ofname` delegates to Matplotlib `savefig` with `bbox_inches="tight"`; PNG, PDF,
and SVG are supported. You can also save or further style the returned figure.
For unusually long labels or large legends, increase `figsize`, adjust wrapping
and anchors, and use tight export bounds. Existing dotplot APIs and defaults are
unchanged.
