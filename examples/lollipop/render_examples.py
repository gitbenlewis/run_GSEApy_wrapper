"""Render deterministic synthetic GSEA lollipop examples; no enrichment is run."""
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Allow this checked-in example to run from any working directory.
REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY.parent))
from run_GSEApy_wrapper import lollipop_from_table


def main():
    output = Path(__file__).resolve().parent / "rendered"
    output.mkdir(exist_ok=True)
    table = pd.DataFrame({
        "Term": list("ABCDE"),
        "Label": [
            "A · Interferon response and regulation of immune cell activation",
            "B · Mitochondrial respiratory chain and cellular energy production",
            "C · Protein transport and intracellular localization",
            "D · Cell cycle checkpoint signaling",
            "E · Extracellular matrix organization",
        ],
        "nes": [2.4, -1.8, 1.3, -0.8, 0.45],
        "fdr": [0, 0.01, 0.05, 0.2, None],
        "fraction": [0.8, 0.6, 0.35, 0.2, 0.45],
    })
    colors = dict(zip("ABCDE", ["#4477aa", "#ee6677", "#228833", "#aa3377", "#66ccee"]))
    common = dict(
        term_label_col="Label", term_colors=colors, pvalue_col="fdr",
        significance_thresholds={0.001: "***", 0.01: "**", 0.05: "*"},
        text_kwargs={"nes_labels": True, "significance": "symbols", "wrap_width": 34},
        legend_kwargs={"significance": {"loc": "upper left", "bbox_to_anchor": (1.12, 1), "ncol": 1,
                                        "fontsize": 8, "title_fontsize": 9, "handlelength": 0}},
    )
    layouts = {
        "equal_size": dict(figsize=(9, 4.8), panel_kwargs={"title": "Synthetic enrichment · equal-size markers"}),
        "significance_panel": dict(figsize=(11, 4.8), significance_panel=True, zero_pvalue_floor=1e-4,
                                   panel_kwargs={"title": "Synthetic enrichment · aligned significance bars"}),
        "fraction_sizes": dict(figsize=(11.5, 5), significance_panel=True, zero_pvalue_floor=1e-4,
                               fraction_col="fraction", marker_kwargs={"fraction_scale": 500, "reference_fraction": 1},
                               panel_kwargs={"title": "Synthetic enrichment · area proportional to fraction"}),
    }
    for name, options in layouts.items():
        axes = lollipop_from_table(table, **common, **options)
        fig = axes[0].figure if isinstance(axes, tuple) else axes.figure
        for extension in ("png", "pdf", "svg"):
            fig.savefig(output / f"{name}.{extension}", dpi=160, bbox_inches="tight")
        plt.close(fig)
    print(f"Rendered three synthetic layouts in PNG, PDF and SVG: {output}")


if __name__ == "__main__":
    main()
