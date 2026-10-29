"""Helpers to plot GSEApy dotplots from saved result tables."""
# version 1 2026-01-14
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional, Sequence
import textwrap
import warnings

import gseapy as gp
import numpy as np
import pandas as pd


def _add_gseapy_aliases(df: pd.DataFrame) -> pd.DataFrame:
    if "FDR q-val" not in df.columns and "fdr" in df.columns:
        df["FDR q-val"] = df["fdr"]
    if "Adjusted P-value" not in df.columns:
        if "FDR q-val" in df.columns:
            df["Adjusted P-value"] = df["FDR q-val"]
        elif "fdr" in df.columns:
            df["Adjusted P-value"] = df["fdr"]
    if "P-value" not in df.columns and "pval" in df.columns:
        df["P-value"] = df["pval"]
    if "NES" not in df.columns and "nes" in df.columns:
        df["NES"] = df["nes"]
    if "ES" not in df.columns and "es" in df.columns:
        df["ES"] = df["es"]
    if "Tag %" not in df.columns and "tag %" in df.columns:
        df["Tag %"] = df["tag %"]
    return df


def _apply_pvalue_zero_floor(
    df: pd.DataFrame,
    pvalue_col: str,
    n_perm_4_zero_pavlues: Optional[int],
) -> pd.DataFrame:
    if pvalue_col not in df.columns:
        raise ValueError(f"pvalue_col '{pvalue_col}' not found in table columns.")
    series = pd.to_numeric(df[pvalue_col], errors="raise")
    zero_mask = series == 0
    if not zero_mask.any():
        return df
    nonzero = series[series > 0]
    if n_perm_4_zero_pavlues is not None:
        if not np.isfinite(n_perm_4_zero_pavlues) or n_perm_4_zero_pavlues <= 0:
            raise ValueError("n_perm_4_zero_pavlues must be > 0.")
        floor = 1.0 / float(n_perm_4_zero_pavlues)
        if not nonzero.empty:
            floor = min(floor, nonzero.min() / 2.0)
    else:
        nonzero = series[series > 0]
        if nonzero.empty:
            raise ValueError(
                f"{pvalue_col} has only zeros; set n_perm_4_zero_pavlues to use 1/N."
            )
        floor = nonzero.min() / 2.0
    df = df.copy()
    if floor <= 0 or not np.isfinite(floor):
        raise ValueError("The display zero floor is not representable as a positive finite float.")
    df[pvalue_col] = series.mask(zero_mask, floor)
    return df


def _dotplot_prepared(df, term_label_char_limit, term_label_col, ordered_terms_list,
                      pvalue_col, n_perm_4_zero_pavlues, options):
    """Select raw probabilities, then render scores and shorten ticks only."""
    nominal = ('pval', 'P-value', 'NOM p-val')
    adjusted = ('fdr', 'FDR q-val', 'Adjusted P-value')
    requested = options.get('column', pvalue_col)
    family = nominal if requested in nominal else adjusted if requested in adjusted else ()
    if family:
        available = [name for name in (requested, *family) if name in df.columns]
        if not available:
            raise ValueError(f"Requested probability column '{requested}' is unavailable.")
        values = pd.to_numeric(df[available[0]], errors='raise')
        for name in available[1:]:
            values = values.fillna(pd.to_numeric(df[name], errors='raise'))
        df[requested] = values
    elif requested not in df.columns:
        raise ValueError(f"Requested column '{requested}' is unavailable.")
    df = _add_gseapy_aliases(df)
    y_col = term_label_col or options.get('y', 'Term')
    options['y'] = y_col
    if ordered_terms_list is not None:
        ordered = [str(term) for term in ordered_terms_list if str(term) in set(df[y_col].astype(str))]
        if not ordered:
            raise ValueError('ordered_terms_list had no matches in the table.')
        df = df[df[y_col].astype(str).isin(ordered)].copy()
        options.setdefault('y_order', ordered)
        options.setdefault('top_term', len(ordered))
    options.setdefault('top_term', 15)
    options.setdefault('cutoff', 1.0)
    if family or 'column' not in options:
        values = pd.to_numeric(df[requested], errors='raise')
        invalid = values.notna() & (~np.isfinite(values) | (values < 0) | (values > 1))
        if invalid.any():
            raise ValueError(f"{requested} must contain probabilities in [0, 1] or missing values.")
        df[requested] = values
        df = df.loc[values <= options['cutoff']].copy()
        if df.empty:
            raise ValueError(f"No enrichment terms meet cutoff={options['cutoff']} for {requested}.")
        n_terms = int(options['top_term'])
        x = options.get('x')
        if x is not None and x in df and not pd.to_numeric(df[x], errors='coerce').notna().all():
            df = pd.concat([group.sort_values(requested, kind='stable').head(n_terms)
                            for _, group in df.groupby(x)], ignore_index=True)
        else:
            df = df.sort_values(requested, kind='stable').head(n_terms)
        display = _apply_pvalue_zero_floor(df, requested, n_perm_4_zero_pavlues)
        # Pretransform to avoid GSEApy clipping subnormal probabilities or
        # reapplying the statistical cutoff to a display-only replacement.
        display_column = f'-log10({requested})'
        df[display_column] = -np.log10(display[requested])
        options['column'] = display_column
        options['cutoff'] = np.inf
        options['top_term'] = len(df)
    else:
        options['column'] = requested
    output = options.pop('ofname', None)
    ax = gp.dotplot(df, **options)
    if ax is not None and term_label_char_limit is not None:
        ax.set_yticks(ax.get_yticks(), [tick.get_text()[:int(term_label_char_limit)] for tick in ax.get_yticklabels()])
    if output is not None and ax is not None:
        ax.figure.savefig(output, bbox_inches='tight', dpi=300)
        return None
    return ax


def load_gseapy_table(
    table_path: str | Path,
    sep: str = ",",
) -> pd.DataFrame:
    """Load a GSEApy output table for dotplotting."""
    return pd.read_csv(table_path, sep=sep)


def load_gseapy_tables(
    table_paths: Iterable[str | Path | pd.DataFrame],
    labels: Optional[Sequence[str]] = None,
    x_label: str = "Dataset",
    sep: str = ",",
) -> pd.DataFrame:
    """Load multiple GSEApy tables and annotate them with a dataset column."""
    tables = list(table_paths)
    paths = [Path(p) for p in tables if not isinstance(p, pd.DataFrame)]
    if labels is None:
        labels = []
        for idx, item in enumerate(tables):
            if isinstance(item, pd.DataFrame):
                labels.append(f"table_{idx}")
            else:
                labels.append(Path(item).stem)
    if len(labels) != len(tables):
        raise ValueError("labels length must match table_paths length")
    frames = []
    for item, label in zip(tables, labels):
        if isinstance(item, pd.DataFrame):
            df = item.copy()
        else:
            df = pd.read_csv(item, sep=sep)
        df = df.assign(**{x_label: label})
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def dotplot_from_table(
    table_path: str | Path | pd.DataFrame,
    sep: str = ",",
    term_label_char_limit: Optional[int] = None,
    term_label_col: Optional[str] = None,
    ordered_terms_list: Optional[Sequence[str]] = None,
    pvalue_col: str = "fdr",
    n_perm_4_zero_pavlues: Optional[int] = 10000,
    **dotplot_kwargs,
):
    """Plot a GSEApy dotplot from a single saved table or DataFrame."""
    if isinstance(table_path, pd.DataFrame):
        df = table_path.copy()
    else:
        df = load_gseapy_table(table_path, sep=sep)
    return _dotplot_prepared(df, term_label_char_limit, term_label_col, ordered_terms_list,
                             pvalue_col, n_perm_4_zero_pavlues, dotplot_kwargs)


def dotplot_from_tables(
    table_paths: Iterable[str | Path],
    labels: Optional[Sequence[str]] = None,
    x_label: str = "Dataset",
    sep: str = ",",
    term_label_char_limit: Optional[int] = None,
    term_label_col: Optional[str] = None,
    ordered_terms_list: Optional[Sequence[str]] = None,
    pvalue_col: str = "fdr",
    n_perm_4_zero_pavlues: Optional[int] = 10000,
    **dotplot_kwargs,
):
    """Plot a GSEApy dotplot across multiple saved tables."""
    df = load_gseapy_tables(
        table_paths,
        labels=labels,
        x_label=x_label,
        sep=sep,
    )
    dotplot_kwargs = {"x": x_label, **dotplot_kwargs}
    return _dotplot_prepared(df, term_label_char_limit, term_label_col, ordered_terms_list,
                             pvalue_col, n_perm_4_zero_pavlues, dotplot_kwargs)


def _lollipop_legend(ax, entries, settings, defaults):
    """Keep independent legends alive and preserve explicit entry ordering."""
    from matplotlib.lines import Line2D

    options = {**defaults, **settings}
    if not options.pop("visible", True) or not entries:
        return
    order = options.pop("order", [key for key, _, _ in entries])
    available = {key: (handle, label) for key, handle, label in entries}
    if len(order) != len(set(order)) or any(key not in available for key in order):
        raise ValueError("Legend order must contain distinct available entry keys.")
    marker_size = options.pop("marker_size", None)
    handles, labels = [], []
    for key in order:
        handle, label = available[key]
        if marker_size is not None and isinstance(handle, Line2D):
            handle.set_markersize(marker_size)
        handles.append(handle)
        labels.append(label)
    if handles:
        previous = ax.get_legend()
        if previous is not None:
            ax.add_artist(previous)
        ax.legend(handles, labels, **options)


def lollipop_from_table(
    table_path: str | Path | pd.DataFrame,
    sep: str = ",",
    term_id_col: str = "Term",
    term_label_col: Optional[str] = None,
    ordered_terms_list: Optional[Sequence[str]] = None,
    term_label_char_limit: Optional[int] = None,
    nes_col: str = "NES",
    pvalue_col: str = "Adjusted P-value",
    term_colors: Optional[dict] = None,
    fraction_col: Optional[str] = None,
    significance_panel: bool = False,
    significance_thresholds: Optional[dict[float, str]] = None,
    zero_pvalue_floor: Optional[float] = None,
    ax=None,
    figsize=None,
    ofname: Optional[str | Path] = None,
    marker_kwargs: Optional[dict] = None,
    stem_kwargs: Optional[dict] = None,
    text_kwargs: Optional[dict] = None,
    panel_kwargs: Optional[dict] = None,
    legend_kwargs: Optional[dict] = None,
):
    """Plot signed NES values without recalculating enrichment statistics.

    Returns one Axes, or (NES Axes, probability Axes) when significance_panel
    is true. Input order is preserved unless ordered_terms_list supplies IDs.
    Missing requested IDs warn; duplicate IDs are ambiguous and raise.

    text_kwargs['significance'] is 'none', 'numeric', or 'symbols'. Symbols
    require significance_thresholds={upper_probability: symbol}; the smallest
    matching inclusive threshold wins. Missing probabilities are 'NA'. Bars
    for reported zeros require zero_pvalue_floor, a display convention only.

    marker_kwargs['size'] is constant area in points squared (default 64).
    With fraction_col, supply 'fraction_scale': area = fraction * scale.
    'reference_fraction' optionally draws an outline behind each marker.
    legend_kwargs has independent 'size', 'significance', and 'color' dictionaries.
    Full styling keys and examples are documented in docs/lollipop.md.
    """
    import matplotlib.pyplot as plt
    from matplotlib.axes import Axes
    from matplotlib.lines import Line2D

    df = table_path.copy() if isinstance(table_path, pd.DataFrame) else load_gseapy_table(table_path, sep=sep)
    df = _add_gseapy_aliases(df)
    marker = dict(marker_kwargs or {})
    stem = {"linewidths": 1.8, **(stem_kwargs or {})}
    text = dict(text_kwargs or {})
    panel = dict(panel_kwargs or {})
    legends = dict(legend_kwargs or {})
    if set(legends) - {"size", "significance", "color"}:
        raise ValueError("legend_kwargs keys must be size, significance, or color.")
    show_nes = text.pop("nes_labels", False)
    annotation = text.pop("significance", "none")
    if annotation not in {"none", "numeric", "symbols"}:
        raise ValueError("text_kwargs['significance'] must be none, numeric, or symbols.")
    nes_format = text.pop("nes_format", ".2f")
    pvalue_format = text.pop("pvalue_format", ".3g")
    offset = text.pop("offset", 7)
    wrap_width = text.pop("wrap_width", 38)
    text = {"fontsize": 9, **text}
    needed = [term_id_col, term_label_col or term_id_col, nes_col]
    use_probability = significance_panel or annotation != "none"
    if use_probability:
        needed.append(pvalue_col)
    if fraction_col is not None:
        needed.append(fraction_col)
    missing_columns = [name for name in dict.fromkeys(needed) if name not in df.columns]
    if missing_columns:
        raise ValueError(f"Required columns not found: {missing_columns}. No nominal-P fallback is used.")
    if df[term_id_col].isna().any():
        raise ValueError(f"{term_id_col} contains missing term IDs.")
    ids = df[term_id_col].astype(str)
    if ids.str.strip().eq("").any() or ids.duplicated().any():
        raise ValueError(f"{term_id_col} must contain unique, nonblank IDs; use TermCollection_Term for multiple collections.")
    df.index = ids
    if ordered_terms_list is not None:
        requested = [str(value) for value in ordered_terms_list]
        if len(requested) != len(set(requested)):
            raise ValueError("ordered_terms_list contains duplicate IDs.")
        absent = [value for value in requested if value not in df.index]
        if absent:
            warnings.warn(f"Requested term IDs not found: {absent}", UserWarning, stacklevel=2)
        df = df.loc[[value for value in requested if value in df.index]]
    if df.empty:
        raise ValueError("No terms remain to plot.")
    ids = df.index.tolist()
    nes = pd.to_numeric(df[nes_col], errors="raise").to_numpy(dtype=float, na_value=np.nan)
    if not np.isfinite(nes).all():
        raise ValueError(f"{nes_col} must contain finite NES values.")
    labels = df[term_label_col or term_id_col].fillna("").astype(str).tolist()
    if term_label_char_limit is not None:
        if term_label_char_limit < 1:
            raise ValueError("term_label_char_limit must be positive.")
        labels = [label[:term_label_char_limit] for label in labels]
    if wrap_width is not None:
        labels = [textwrap.fill(label, width=wrap_width) for label in labels]
    color_default = marker.pop("color", "#4477aa")
    colors = [(term_colors or {}).get(term, color_default) for term in ids]
    if "s" in marker or "c" in marker:
        raise ValueError("Use marker_kwargs['size'] and term_colors instead of s or c.")
    size = marker.pop("size", 64)
    scale = marker.pop("fraction_scale", None)
    reference = marker.pop("reference_fraction", None)
    size_values = marker.pop("legend_fractions", [0.25, 0.5, 1.0])
    if fraction_col is not None:
        if scale is None or not np.isfinite(scale) or scale <= 0:
            raise ValueError("fraction_col requires a positive finite marker_kwargs['fraction_scale'].")
        fractions = pd.to_numeric(df[fraction_col], errors="raise").to_numpy(dtype=float, na_value=np.nan)
        if not np.isfinite(fractions).all() or ((fractions < 0) | (fractions > 1)).any():
            raise ValueError(f"{fraction_col} must contain finite fractions in [0, 1].")
        if any(not np.isfinite(value) or value < 0 or value > 1 for value in size_values):
            raise ValueError("legend_fractions must be finite fractions in [0, 1].")
        sizes = fractions * scale
    else:
        if not np.isscalar(size) or not np.isfinite(size) or size <= 0:
            raise ValueError("marker_kwargs['size'] must be a positive finite scalar area.")
        if reference is not None:
            raise ValueError("reference_fraction requires fraction_col.")
        sizes = np.full(len(df), size)
    if reference is not None and (not np.isfinite(reference) or reference < 0 or reference > 1):
        raise ValueError("reference_fraction must be in [0, 1].")
    probabilities = np.full(len(df), np.nan)
    if use_probability:
        probabilities = pd.to_numeric(df[pvalue_col], errors="raise").to_numpy(dtype=float, na_value=np.nan)
        available = ~np.isnan(probabilities)
        if not np.isfinite(probabilities[available]).all() or ((probabilities[available] < 0) | (probabilities[available] > 1)).any():
            raise ValueError(f"{pvalue_col} must contain probabilities in [0, 1] or missing values.")
    thresholds = []
    if annotation == "symbols":
        if not significance_thresholds:
            raise ValueError("Symbol annotations require significance_thresholds={upper_probability: symbol}.")
        thresholds = sorted(significance_thresholds.items())
        if any(not np.isfinite(value) or value < 0 or value > 1 for value, _ in thresholds):
            raise ValueError("Significance thresholds must be finite probabilities in [0, 1].")
    display_probability = probabilities.copy()
    zero_mask = probabilities == 0
    if significance_panel and zero_mask.any():
        if zero_pvalue_floor is None or not np.isfinite(zero_pvalue_floor) or not 0 < zero_pvalue_floor < 1:
            raise ValueError("Zero probabilities require an explicit zero_pvalue_floor in (0, 1) for finite bars.")
        positive = probabilities[probabilities > 0]
        if positive.size and zero_pvalue_floor > positive.min():
            raise ValueError("zero_pvalue_floor must not exceed the smallest positive selected probability.")
        display_probability[zero_mask] = zero_pvalue_floor

    title = panel.pop("title", None)
    nes_xlabel = panel.pop("nes_xlabel", "Normalized enrichment score (NES)")
    probability_xlabel = panel.pop("significance_xlabel", f"−log10({pvalue_col})")
    bar_height = panel.pop("bar_height", 0.55)
    width_ratios = panel.pop("width_ratios", (1.2, 1))
    wspace = panel.pop("wspace", 0.12)
    xmargin = panel.pop("xmargin", 0.25 if show_nes or annotation != "none" else 0.08)
    zero_line = {"color": "0.5", "linewidth": 0.8, "zorder": 0, **panel.pop("zero_line_kwargs", {})}
    if panel:
        raise ValueError(f"Unknown panel_kwargs: {sorted(panel)}")
    if ax is None:
        height = max(2.4, 0.48 * sum(max(1, label.count('\n') + 1) for label in labels) + 1.4)
        fig, axes = plt.subplots(
            1, 2 if significance_panel else 1,
            figsize=figsize or ((10 if significance_panel else 7), height),
            layout="constrained",
            gridspec_kw={"width_ratios": width_ratios, "wspace": wspace} if significance_panel else None,
        )
        nes_ax, probability_ax = axes if significance_panel else (axes, None)
    elif significance_panel:
        if isinstance(ax, Axes) or len(ax) != 2 or not all(isinstance(item, Axes) for item in ax):
            raise ValueError("significance_panel=True requires two Axes from the same figure.")
        nes_ax, probability_ax = ax
        if nes_ax is probability_ax or nes_ax.figure is not probability_ax.figure:
            raise ValueError("Provide two distinct Axes from the same figure.")
        fig = nes_ax.figure
    else:
        if not isinstance(ax, Axes):
            raise ValueError("Single-panel plotting requires one Axes.")
        nes_ax, probability_ax, fig = ax, None, ax.figure
    y = np.arange(len(df))
    nes_ax.hlines(y, 0, nes, colors=colors, **stem)
    nes_ax.axvline(0, **zero_line)
    if reference is not None:
        nes_ax.scatter(nes, y, s=reference * scale, facecolors="none", edgecolors=colors, linewidths=0.8, zorder=2)
    nes_ax.scatter(nes, y, s=sizes, c=colors, **{"zorder": 3, **marker})
    nes_ax.set_yticks(y, labels)
    nes_ax.set_xlabel(nes_xlabel + (f"\nAnnotations: {pvalue_col}" if use_probability and not significance_panel else ""))
    if title is not None:
        nes_ax.set_title(title)
    nes_ax.margins(x=xmargin)
    if significance_panel:
        available = np.isfinite(display_probability)
        bar_values = -np.log10(display_probability[available])
        probability_ax.barh(y[available], bar_values, height=bar_height,
                            color=[color for color, keep in zip(colors, available) if keep])
        probability_ax.set_yticks(y)
        probability_ax.tick_params(axis="y", left=False, labelleft=False)
        probability_ax.set_xlim(left=0)
        probability_ax.margins(x=max(xmargin, 0.15))
        if zero_mask.any():
            probability_xlabel += f"\nReported zeros drawn at {zero_pvalue_floor:g} (display only)"
        probability_ax.set_xlabel(probability_xlabel)
    for target in (nes_ax, probability_ax):
        if target is not None:
            target.set_ylim(len(df) - 0.5, -0.5)
            target.spines[["top", "right"]].set_visible(False)
    for row, (value, probability) in enumerate(zip(nes, probabilities)):
        sign = 1 if value >= 0 else -1
        if show_nes:
            # Offset from the circle edge, including an optional reference outline.
            radius = np.sqrt(max(sizes[row], reference * scale if reference is not None else 0)) / 2
            nes_ax.annotate(format(value, nes_format), (value, row), xytext=(sign * (offset + radius), 0),
                            textcoords="offset points", ha="left" if sign > 0 else "right", va="center", **text)
        if annotation != "none" or (significance_panel and np.isnan(probability)):
            if np.isnan(probability):
                label = "NA"
            elif annotation == "numeric":
                label = format(probability, pvalue_format)
            elif annotation == "symbols":
                label = next((str(symbol) for threshold, symbol in thresholds if probability <= threshold), "ns")
            else:
                continue
            if significance_panel:
                x = 0 if np.isnan(probability) else -np.log10(display_probability[row])
                probability_ax.annotate(label, (x, row), xytext=(offset, 0), textcoords="offset points",
                                        ha="left", va="center", **text)
            else:
                # Put probability annotations in a separate column to avoid NES-label overlap.
                nes_ax.text(1.02, row, label, transform=nes_ax.get_yaxis_transform(), va="center", **text)

    legend_ax = probability_ax or nes_ax
    if fraction_col is not None:
        from matplotlib.font_manager import FontProperties

        entries = [(value, Line2D([], [], linestyle="none", marker="o", color=color_default,
                                 markersize=np.sqrt(value * scale)), f"{value:g}") for value in size_values]
        size_options = legends.get("size", {})
        legend_fontsize = FontProperties(size=size_options.get("fontsize", plt.rcParams["legend.fontsize"])).get_size_in_points()
        diameter = size_options.get("marker_size", np.sqrt(max(size_values, default=0) * scale))
        diameter *= size_options.get("markerscale", plt.rcParams["legend.markerscale"])
        _lollipop_legend(legend_ax, entries, legends.get("size", {}),
                         {"title": fraction_col, "loc": "upper left", "bbox_to_anchor": (1.12, 0.5),
                          "labelspacing": max(0.5, diameter / legend_fontsize), "frameon": False})
    if annotation == "symbols":
        entries = [(threshold, Line2D([], [], linestyle="none", marker="", color="none"),
                    f"{symbol}: {pvalue_col} ≤ {threshold:g}") for threshold, symbol in thresholds]
        _lollipop_legend(legend_ax, entries, legends.get("significance", {}),
                         {"title": "Significance", "loc": "upper left", "bbox_to_anchor": (1.12, 1), "frameon": False})
    color_entries = [(term, Line2D([], [], linestyle="none", marker="o", color=color), label)
                     for term, color, label in zip(ids, colors, labels)]
    _lollipop_legend(legend_ax, color_entries, legends.get("color", {}),
                     {"visible": False, "title": "Terms", "loc": "upper left", "bbox_to_anchor": (1.02, 0), "frameon": False})
    if ofname is not None:
        fig.savefig(ofname, bbox_inches="tight")
    return (nes_ax, probability_ax) if significance_panel else nes_ax


__all__ = [
    "load_gseapy_table",
    "load_gseapy_tables",
    "dotplot_from_table",
    "dotplot_from_tables",
    "lollipop_from_table",
]
