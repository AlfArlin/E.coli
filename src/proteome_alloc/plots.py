"""Figures. One function per figure; each returns the Figure and, when a path
is given, saves it (PNG, 200 dpi) so results never depend on manual export.

Encoding rules used throughout
------------------------------
* dataset identity = fixed colour AND fixed marker shape (never colour alone);
* one y-axis per panel; magnitudes on a single-hue sequential scale;
* sign of an effect = blue (rises with growth) / orange (falls), grey = not significant.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from .ingest import GROWTH

SURFACE, INK, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#898781", "#e1e0d9"
UP, DOWN, NEUTRAL = "#2a78d6", "#eb6834", "#898781"
DATASET_STYLE = {            # fixed order: colour follows the dataset, never its rank
    "schmidt_2016": ("#2a78d6", "o"), "mori_2021": ("#eb6834", "s"),
    "wu_2023": ("#1baf7a", "^"), "zhu_2025": ("#eda100", "D"),
}
DIVERGING = LinearSegmentedColormap.from_list("orange_grey_blue", ["#c4501f", "#eb6834", "#f0efec", "#2a78d6", "#184f95"])
SEQUENTIAL = LinearSegmentedColormap.from_list("blues", ["#f3f7fd", "#86b6ef", "#2a78d6", "#0d366b"])
MU_LABEL = "Specific growth rate (h$^{-1}$)"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "axes.axisbelow": True, "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "legend.frameon": False, "figure.dpi": 110,
})


def _style(dataset: str):
    return DATASET_STYLE.get(str(dataset), (MUTED, "X"))


def _save(fig, path):
    if path is not None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=200, bbox_inches="tight")
    return fig


def _scatter_by_dataset(ax, x, y, datasets, **kwargs):
    for name in pd.unique(datasets):
        colour, marker = _style(name)
        mask = np.asarray(datasets) == name
        ax.scatter(np.asarray(x)[mask], np.asarray(y)[mask], c=colour, marker=marker, s=46,
                   edgecolor=INK, linewidth=0.6, label=name, zorder=3, **kwargs)


# --------------------------------------------------------------------------- #
def plot_design(obs: pd.DataFrame, path=None):
    """Growth rate of every sample per dataset; filled = in scope, open = excluded."""
    datasets = list(pd.unique(obs["Dataset_Source"]))
    fig, ax = plt.subplots(figsize=(7.5, 0.75 * len(datasets) + 1.6))
    rng = np.random.default_rng(0)                      # display jitter only
    for row, name in enumerate(datasets):
        sub = obs[(obs["Dataset_Source"] == name) & obs[GROWTH].notna()]
        colour, marker = _style(name)
        jitter = rng.uniform(-0.18, 0.18, len(sub))
        inside = sub["In_Scope"].to_numpy(bool)
        ax.scatter(sub[GROWTH][inside], row + jitter[inside], c=colour, marker=marker, s=46,
                   edgecolor=INK, linewidth=0.6, zorder=3)
        ax.scatter(sub[GROWTH][~inside], row + jitter[~inside], facecolors="none", marker=marker,
                   s=46, edgecolor=MUTED, linewidth=0.9, zorder=2)
        n_in, n_all = int(obs[obs["Dataset_Source"] == name]["In_Scope"].sum()), int((obs["Dataset_Source"] == name).sum())
        ax.text(1.01, row, f"{n_in} of {n_all} in scope", transform=ax.get_yaxis_transform(),
                va="center", color="#52514e")
    ax.set_yticks(range(len(datasets)), datasets)
    ax.set_ylim(-0.6, len(datasets) - 0.4)
    ax.invert_yaxis()
    ax.set_xlabel(MU_LABEL)
    ax.set_title("Design check: growth-rate coverage per dataset (filled = nutrient-quality scope)")
    ax.grid(axis="y", visible=False)
    return _save(fig, path)


def plot_distributions(X: pd.DataFrame, obs: pd.DataFrame, path=None):
    """Box plot per sample of log10 mass fractions (quantified proteins only).
    Shows dynamic range and depth of each platform before any modelling."""
    order = obs.sort_values(["Dataset_Source"], kind="stable").index
    fig, ax = plt.subplots(figsize=(max(7.5, 0.2 * len(order) + 2), 3.6))
    data = [np.log10(X.loc[s][X.loc[s] > 0].to_numpy(float)) for s in order]
    boxes = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False,
                       medianprops={"color": INK, "linewidth": 1.2}, whiskerprops={"color": MUTED},
                       capprops={"color": MUTED})
    seen = {}
    for patch, sample in zip(boxes["boxes"], order):
        name = obs.loc[sample, "Dataset_Source"]
        patch.set(facecolor=_style(name)[0], edgecolor=INK, linewidth=0.6)
        seen.setdefault(name, patch)
    ax.set_xticks([])
    ax.set_xlabel("Samples (grouped by dataset)")
    ax.set_ylabel("log$_{10}$ protein mass fraction")
    ax.set_title("Distribution of quantified protein mass fractions per sample")
    ax.legend(seen.values(), seen.keys(), ncol=len(seen), loc="lower center", bbox_to_anchor=(0.5, -0.32))
    ax.grid(axis="x", visible=False)
    return _save(fig, path)


def plot_coverage(coverage: pd.DataFrame, obs: pd.DataFrame, path=None):
    """Proteins detected and share of proteome mass with a COG category, per sample."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
    datasets = obs.loc[coverage.index, "Dataset_Source"]
    for ax, column, label in zip(axes, ["n_proteins_detected", "mass_fraction_with_COG"],
                                 ["Proteins quantified per sample", "Proteome mass with a COG category"]):
        positions, ticks = [], []
        for i, name in enumerate(pd.unique(datasets)):
            values = coverage.loc[datasets == name, column]
            colour, marker = _style(name)
            x = i + np.linspace(-0.25, 0.25, len(values)) if len(values) > 1 else np.array([i])
            ax.scatter(x, values, c=colour, marker=marker, s=40, edgecolor=INK, linewidth=0.6, zorder=3)
            positions.append(i), ticks.append(name)
        ax.set_xticks(positions, ticks, rotation=20, ha="right")
        ax.set_title(label)
        ax.grid(axis="x", visible=False)
    axes[1].set_ylim(min(0.5, coverage["mass_fraction_with_COG"].min() - 0.05), 1.0)
    return _save(fig, path)


def plot_composition(Fc: pd.DataFrame, info: pd.DataFrame, labels: pd.Series, path=None, top: int = 14):
    """Heatmap of category mass fractions; conditions ordered by growth within dataset."""
    order = info.sort_values(["Dataset_Source", GROWTH]).index
    categories = Fc.mean().sort_values(ascending=False).index[:top]
    data = Fc.loc[order, categories].T * 100
    fig, ax = plt.subplots(figsize=(0.33 * len(order) + 4.5, 0.3 * len(categories) + 1.8))
    image = ax.imshow(data.to_numpy(), aspect="auto", cmap=SEQUENTIAL, vmin=0)
    ax.set_yticks(range(len(categories)), [f"{c}  {labels.get(c, '')[:38]}" for c in categories])
    ax.set_xticks(range(len(order)), [f"{info.loc[c, GROWTH]:.2f}" for c in order], rotation=90)
    ax.set_xlabel(MU_LABEL + " - conditions sorted within dataset")
    boundaries = np.flatnonzero(info.loc[order, "Dataset_Source"].to_numpy()[1:] !=
                                info.loc[order, "Dataset_Source"].to_numpy()[:-1])
    for b in boundaries:
        ax.axvline(b + 0.5, color=SURFACE, linewidth=3)
    starts = np.r_[0, boundaries + 1]
    for start, end in zip(starts, np.r_[boundaries, len(order) - 1]):
        ax.text((start + end) / 2, -0.8, info.loc[order[start], "Dataset_Source"], ha="center", va="bottom")
    ax.grid(False)
    fig.colorbar(image, ax=ax, label="% of proteome mass", fraction=0.03, pad=0.02)
    ax.set_title(f"COG allocation per condition (top {len(categories)} categories)", pad=22)
    return _save(fig, path)


def plot_growth_laws(Fc: pd.DataFrame, info: pd.DataFrame, table: pd.DataFrame, labels: pd.Series,
                     path=None, n_panels: int = 8, batch_key: str = "Dataset_Source"):
    """Small multiples: category fraction vs growth rate with the common-slope fit."""
    tested = table[table["tested"]]
    chosen = tested.reindex(tested["slope_per_h-1"].abs().sort_values(ascending=False).index).index[:n_panels]
    columns = 4
    rows = int(np.ceil(len(chosen) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(3.1 * columns, 2.7 * rows), squeeze=False)
    mu, groups = info[GROWTH].to_numpy(float), info[batch_key].astype(str).to_numpy()
    for ax, category in zip(axes.ravel(), chosen):
        y = Fc[category].to_numpy(float) * 100
        slope = table.loc[category, "slope_per_h-1"] * 100
        for g in np.unique(groups):                       # one intercept per dataset, shared slope
            mask = groups == g
            xs = np.array([mu[mask].min(), mu[mask].max()])
            ax.plot(xs, y[mask].mean() + slope * (xs - mu[mask].mean()), color=_style(info.loc[mask, "Dataset_Source"].iloc[0])[0],
                    linewidth=2, zorder=2)
        _scatter_by_dataset(ax, mu, y, info["Dataset_Source"].to_numpy())
        q = table.loc[category, "q_BH"]
        ax.set_title(f"{category}  {labels.get(category, '')[:34]}\n", fontsize=8.5)
        ax.text(0.0, 1.03, f"slope {slope:+.2f} points per h$^{{-1}}$,  q = {q:.3g}", transform=ax.transAxes,
                va="bottom", fontsize=8, color="#52514e")
    for ax in axes.ravel()[len(chosen):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel(MU_LABEL)
    for ax in axes[:, 0]:
        ax.set_ylabel("% of proteome mass")
    handles, names = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, names, loc="upper center", ncol=len(names), bbox_to_anchor=(0.5, 1.03))
    fig.tight_layout()
    return _save(fig, path)


def plot_slopes(table: pd.DataFrame, labels: pd.Series, path=None):
    """Forest plot: slope and bootstrap CI for every tested category."""
    data = table[table["tested"]].sort_values("slope_per_h-1")
    fig, ax = plt.subplots(figsize=(7.5, 0.3 * len(data) + 1.5))
    y = np.arange(len(data))
    slope, low, high = (data[c].to_numpy() * 100 for c in ("slope_per_h-1", "ci_low", "ci_high"))
    colours = [NEUTRAL if not sig else (UP if s > 0 else DOWN) for s, sig in zip(slope, data["significant"])]
    ax.hlines(y, low, high, color=colours, linewidth=2, zorder=2)
    ax.scatter(slope, y, c=colours, s=44, edgecolor=INK, linewidth=0.6, zorder=3)
    ax.axvline(0, color="#52514e", linewidth=1)
    ax.set_yticks(y, [f"{c}  {labels.get(c, '')[:42]}" for c in data.index])
    ax.set_xlabel("Change in proteome share per unit growth rate (percentage points per h$^{-1}$)")
    ax.set_title("Growth-rate dependence of COG categories (95% bootstrap CI; coloured = FDR-significant)")
    ax.grid(axis="y", visible=False)
    return _save(fig, path)


def plot_asca(result, path=None, top_n: int = 20, batch_key: str = "Dataset_Source"):
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8))
    vt = result.variance_table.drop(index="Total (centred)")
    ax = axes[0, 0]
    colours = {"Dataset": NEUTRAL, "Growth rate | dataset": UP, "Dataset x growth rate": DOWN, "Residual": "#c3c2b7"}
    bars = ax.bar([i.replace(" | ", " |\n").replace(" x ", " x\n") for i in vt.index], vt["SS_fraction"] * 100,
                  color=[colours.get(i, NEUTRAL) for i in vt.index], width=0.55)
    ax.bar_label(bars, fmt="%.1f%%", padding=2)
    ax.set_ylabel("% of total centred sum of squares")
    p_inter = getattr(result, "p_interaction", np.nan)
    ax.set_title(f"A. Variance partition (growth p = {result.p_growth:.4f}; interaction p = {p_inter:.4f})")
    ax.grid(axis="x", visible=False)

    ax = axes[0, 1]
    ax.hist(result.null_distribution, bins=40, color="#c3c2b7")
    ax.axvline(result.observed_ss_growth, color=UP, linewidth=2)
    ax.text(result.observed_ss_growth, ax.get_ylim()[1] * 0.95, " observed", color=INK, va="top")
    ax.set_xlabel("SS(growth | dataset) under within-dataset permutation")
    ax.set_title("B. Permutation null distribution")

    ax = axes[1, 0]
    s = result.scores
    _scatter_by_dataset(ax, s[GROWTH], s["effect_plus_residual_score"], s["Dataset_Source"].to_numpy())
    ax.set_xlabel(MU_LABEL)
    ax.set_ylabel("Score on the growth component (effect + residual)")
    ax.set_title("C. Conditions along the growth component")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    table = result.protein_table
    extreme = pd.concat([table.head(top_n - top_n // 2), table.tail(top_n // 2)])
    values = extreme["slope_log_per_h-1"]
    ax.barh([f"{g} ({b})" for g, b in zip(extreme["gene"].astype(str), extreme.index)], values,
            color=[UP if v > 0 else DOWN for v in values], height=0.65)
    ax.invert_yaxis()
    ax.axvline(0, color="#52514e", linewidth=1)
    ax.set_xlabel("Within-dataset slope (log$_{10}$ mass fraction per h$^{-1}$)")
    ax.set_title(f"D. {top_n} most growth-dependent proteins")
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    return _save(fig, path)


def plot_prediction(result: dict, path=None):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    table = result["predictions"]
    ax = axes[0]
    _scatter_by_dataset(ax, table["observed"], table["predicted"], table["Dataset_Source"].to_numpy())
    limits = [min(table["observed"].min(), table["predicted"].min()) - 0.1,
              max(table["observed"].max(), table["predicted"].max()) + 0.1]
    ax.plot(limits, limits, color="#52514e", linewidth=1, linestyle="--", zorder=1)
    ax.set_xlim(limits), ax.set_ylim(limits)
    ax.set_xlabel("Observed growth rate (h$^{-1}$)")
    ax.set_ylabel("Predicted when the whole dataset is held out (h$^{-1}$)")
    ax.set_title(f"A. Leave-one-dataset-out prediction (RMSE {result['rmse_model']:.2f} vs "
                 f"baseline {result['rmse_baseline']:.2f})")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.hist(result["null_rmse"], bins=25, color="#c3c2b7")
    ax.axvline(result["rmse_model"], color=UP, linewidth=2)
    ax.axvline(result["rmse_baseline"], color=DOWN, linewidth=2, linestyle="--")
    top = ax.get_ylim()[1]
    ax.text(result["rmse_model"], top * 0.97, " model", va="top")
    ax.text(result["rmse_baseline"], top * 0.85, " mean baseline", va="top")
    ax.set_xlabel("Pooled RMSE (h$^{-1}$), growth rates permuted within dataset")
    ax.set_title(f"B. Permutation null (p = {result['p_permutation']:.3f})")
    fig.tight_layout()
    return _save(fig, path)


# --------------------------------------------------------------------------- #
# combined views (all datasets together, batch offsets removed for display)
# --------------------------------------------------------------------------- #
def plot_condition_heatmap(F_adj: pd.DataFrame, info: pd.DataFrame, labels: pd.Series, cond_labels: pd.Series,
                           path=None, top: int = 16, min_share: float = 0.001):
    """All conditions of all datasets on ONE axis, ordered by growth rate.
    Values = batch-adjusted share of each category, standardised per category (z-score),
    so small categories are as visible as the ribosome. Row label gives the mean share;
    the percentage-point changes are in the proteome-shift table."""
    order = info.sort_values(GROWTH).index
    means = F_adj.mean()
    categories = means[means >= min_share].sort_values(ascending=False).index[:top]   # skip near-empty categories
    values = F_adj.loc[order, categories]
    data = ((values - values.mean()) / values.std(ddof=1).replace(0, np.nan)).T.fillna(0)   # z-score per category
    limit = float(np.nanmax(np.abs(data.to_numpy()))) or 1.0
    fig, (strip, ax) = plt.subplots(2, 1, figsize=(0.34 * len(order) + 5, 0.32 * len(categories) + 3.2),
                                    gridspec_kw={"height_ratios": [0.6, len(categories)], "hspace": 0.04},
                                    sharex=True)
    for i, name in enumerate(info.loc[order, "Dataset_Source"]):
        strip.add_patch(plt.Rectangle((i - 0.5, 0), 1, 1, facecolor=_style(name)[0], edgecolor=SURFACE, linewidth=1.5))
    strip.set_xlim(-0.5, len(order) - 0.5), strip.set_ylim(0, 1)
    strip.set_yticks([]), strip.grid(False)
    for side in strip.spines.values():
        side.set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=_style(n)[0]) for n in pd.unique(info.loc[order, "Dataset_Source"])]
    strip.legend(handles, pd.unique(info.loc[order, "Dataset_Source"]), ncol=4, loc="lower left",
                 bbox_to_anchor=(0, 1.05), fontsize=8)
    strip.set_title("Proteome allocation per condition - all datasets combined, sorted by growth rate", pad=26)
    image = ax.imshow(data.to_numpy(), aspect="auto", cmap=DIVERGING, vmin=-limit, vmax=limit)
    ax.set_yticks(range(len(categories)),
                  [f"{c}  {labels.get(c, '')[:34]}  ({values[c].mean() * 100:.1f}%)" for c in categories])
    ax.set_xticks(range(len(order)), [cond_labels.loc[c] for c in order], rotation=90)
    ax.set_xlabel("Condition (specific growth rate, h$^{-1}$)")
    ax.grid(False)
    fig.colorbar(image, ax=[strip, ax], label="z-score within category (blue = above its mean share)",
                 fraction=0.025, pad=0.01)
    return _save(fig, path)


def plot_combined_growth_laws(F_adj: pd.DataFrame, info: pd.DataFrame, table: pd.DataFrame, labels: pd.Series,
                              path=None, n_panels: int = 8):
    """Pooled view: batch-adjusted fraction vs growth rate, ONE line for all datasets."""
    tested = table[table["tested"]]
    chosen = tested.reindex(tested["slope_per_h-1"].abs().sort_values(ascending=False).index).index[:n_panels]
    columns = 4
    rows = int(np.ceil(len(chosen) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(3.1 * columns, 2.7 * rows), squeeze=False)
    mu = info[GROWTH].to_numpy(float)
    xs = np.array([mu.min(), mu.max()])
    for ax, category in zip(axes.ravel(), chosen):
        y = F_adj[category].to_numpy(float) * 100
        slope = table.loc[category, "slope_per_h-1"] * 100
        ax.plot(xs, y.mean() + slope * (xs - mu.mean()), color=INK, linewidth=2, zorder=2)
        _scatter_by_dataset(ax, mu, y, info["Dataset_Source"].to_numpy())
        ax.set_title(f"{category}  {labels.get(category, '')[:34]}\n", fontsize=8.5)
        ax.text(0.0, 1.03, f"slope {slope:+.2f} points per h$^{{-1}}$,  q = {table.loc[category, 'q_BH']:.3g}",
                transform=ax.transAxes, va="bottom", fontsize=8, color="#52514e")
    for ax in axes.ravel()[len(chosen):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel(MU_LABEL)
    for ax in axes[:, 0]:
        ax.set_ylabel("% of proteome (offset removed)")
    handles, names = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, names, loc="upper center", ncol=len(names), bbox_to_anchor=(0.5, 1.03))
    fig.tight_layout()
    return _save(fig, path)


def plot_pca_before_after(before, after, ratio_before, ratio_after, info: pd.DataFrame, path=None):
    """2 x 2: PCA of the combined protein matrix before / after removing dataset offsets.
    Top row coloured by dataset, bottom row by growth rate (marker = dataset in both)."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    mu = info.loc[before.index, GROWTH].to_numpy(float)
    datasets = info.loc[before.index, "Dataset_Source"].to_numpy()
    for column, (scores, ratio, title) in enumerate([(before, ratio_before, "Before: log10 proteome, all datasets"),
                                                     (after, ratio_after, "After: dataset offsets removed")]):
        ax = axes[0, column]
        _scatter_by_dataset(ax, scores["PC1"], scores["PC2"], datasets)
        ax.set_title(f"{'AB'[column]}. {title} - colour = dataset")
        ax.legend(fontsize=8)
        ax = axes[1, column]
        for name in pd.unique(datasets):
            mask = datasets == name
            sc = ax.scatter(scores["PC1"][mask], scores["PC2"][mask], c=mu[mask], cmap=SEQUENTIAL,
                            vmin=mu.min(), vmax=mu.max(), marker=_style(name)[1], s=55, edgecolor=INK, linewidth=0.6)
        ax.set_title(f"{'CD'[column]}. {title.split(':')[0]} - colour = growth rate")
        for row in (0, 1):
            axes[row, column].set_xlabel(f"PC1 ({ratio[0]:.1%} of variance)")
            axes[row, column].set_ylabel(f"PC2 ({ratio[1]:.1%} of variance)")
    fig.colorbar(sc, ax=axes[1, :], label=MU_LABEL, fraction=0.03, pad=0.02)
    return _save(fig, path)


def plot_condition_dendrograms(before: pd.DataFrame, after: pd.DataFrame, info: pd.DataFrame,
                               cond_labels: pd.Series, path=None):
    """Ward clustering of conditions before / after removing dataset offsets.
    Leaf label colour = dataset. Before: conditions cluster by study; after: by growth/medium."""
    from scipy.cluster.hierarchy import dendrogram, linkage
    fig, axes = plt.subplots(1, 2, figsize=(12, 0.26 * len(before) + 2.2))
    for ax, data, title in ((axes[0], before, "A. Before (log10 proteome)"),
                            (axes[1], after, "B. After removing dataset offsets")):
        centred = data - data.mean(axis=0)
        tree = linkage(centred.to_numpy(float), method="ward")
        result = dendrogram(tree, orientation="left", labels=[cond_labels.loc[i] for i in data.index],
                            ax=ax, color_threshold=0, above_threshold_color="#52514e", leaf_font_size=8)
        lookup = dict(zip(cond_labels.loc[data.index], info.loc[data.index, "Dataset_Source"]))
        for tick in ax.get_yticklabels():
            tick.set_color(_style(lookup[tick.get_text()])[0])
            tick.set_fontweight("bold")
        ax.set_title(title)
        ax.set_xlabel("Ward linkage distance")
        ax.grid(False)
        for side in ("left", "top", "right"):
            ax.spines[side].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=_style(n)[0]) for n in pd.unique(info["Dataset_Source"])]
    fig.legend(handles, pd.unique(info["Dataset_Source"]), ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout()
    return _save(fig, path)
