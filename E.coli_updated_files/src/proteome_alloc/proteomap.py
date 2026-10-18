"""Proteomap-style treemaps: where the proteome mass goes in two contrasting conditions.

Each panel is one growth condition; area = share of proteome mass. Tiles are
nested  COG super-group -> COG category -> protein  (squarified treemap,
Bruls et al. 2000; same visual idea as Liebermeister et al. 2014 proteomaps,
with rectangles instead of Voronoi cells). Colour families follow the four COG
super-groups: information processing (blue), cellular processes (green),
metabolism (orange/brown), poorly characterised (grey), not in COG (dark grey).

Which conditions are shown
--------------------------
Two conditions are only comparable when they come from the SAME measurement
batch; otherwise part of the visible difference would be platform, not biology.
With `proteomap.conditions: auto` the batch with the widest growth-rate range is
chosen and its slowest and fastest conditions are drawn. Any pair can be fixed
in config.yml instead. The maps illustrate; the statistical evidence for each
category is the slope analysis of notebook 03 (marked with * in the bar chart).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

from .config import Config
from .ingest import GROWTH
from .plots import INK, SURFACE, UP, DOWN, MU_LABEL, _save

SUPERGROUPS = {
    "Information storage and processing": list("JAKLB"),
    "Cellular processes and signalling": list("DYVTMNZWUOX"),
    "Metabolism": list("CGEFHIPQ"),
    "Poorly characterised": list("RS"),
    "Not in COG": ["-"],
}
COLOURS = {
    "J": "#2a78d6", "K": "#5598e7", "L": "#86b6ef", "A": "#1c5cab", "B": "#184f95",
    "O": "#1baf7a", "M": "#199e70", "U": "#4cc79a", "T": "#7dd6b4", "V": "#0f7a55",
    "N": "#a6e3cb", "D": "#2e8b6e", "Y": "#5fb79a", "Z": "#5fb79a", "W": "#5fb79a", "X": "#5fb79a",
    "E": "#eda100", "G": "#eb6834", "C": "#c4501f", "F": "#f2b84b", "H": "#d98a1c",
    "I": "#e8935e", "P": "#a8642a", "Q": "#f5c98a",
    "R": "#b9b8b2", "S": "#9c9b95", "-": "#5c5b57",
}
SHORT = {
    "J": "Translation & ribosome", "K": "Transcription", "L": "Replication & repair",
    "A": "RNA processing", "B": "Chromatin", "O": "Chaperones & turnover", "M": "Cell envelope",
    "U": "Secretion", "T": "Signal transduction", "V": "Defence", "N": "Motility",
    "D": "Cell cycle", "C": "Energy production", "G": "Carbohydrate metabolism",
    "E": "Amino acid metabolism", "F": "Nucleotide metabolism", "H": "Coenzyme metabolism",
    "I": "Lipid metabolism", "P": "Inorganic ion transport", "Q": "Secondary metabolites",
    "R": "General function", "S": "Function unknown", "-": "Not in COG",
}


# --------------------------------------------------------------------------- #
# squarified treemap (Bruls, Huizing & van Wijk 2000)
# --------------------------------------------------------------------------- #
def _worst(row, side):
    """Worst aspect ratio of a row of rectangles (squarified treemap criterion)."""
    total = sum(row)
    return max(max(side ** 2 * r / total ** 2, total ** 2 / (side ** 2 * r)) for r in row)


def squarify(values, x, y, w, h):
    """Rectangles (x, y, w, h) for positive `values` (already in display order,
    ideally descending) filling the box. Returns a list aligned with `values`."""
    values = [float(v) for v in values]
    total = sum(values)
    if total <= 0 or w <= 0 or h <= 0:
        return [(x, y, 0.0, 0.0) for _ in values]
    areas = [v * w * h / total for v in values]
    rects, row, start = [], [], 0
    i = 0
    while i < len(areas):
        side = min(w, h)
        candidate = row + [areas[i]]
        if not row or _worst(candidate, side) <= _worst(row, side):
            row = candidate
            i += 1
            continue
        x, y, w, h = _lay_row(row, x, y, w, h, rects)
        row = []
    if row:
        _lay_row(row, x, y, w, h, rects)
    return rects


def _lay_row(row, x, y, w, h, rects):
    """Place one row of treemap rectangles and return the remaining free area."""
    total = sum(row)
    if w >= h:                                   # row is a column on the left
        width = total / h
        offset = y
        for area in row:
            height = area / width
            rects.append((x, offset, width, height))
            offset += height
        return x + width, y, w - width, h
    height = total / w                           # row along the bottom
    offset = x
    for area in row:
        width = area / height
        rects.append((offset, y, width, height))
        offset += width
    return x, y + height, w, h - height


# --------------------------------------------------------------------------- #
# data
# --------------------------------------------------------------------------- #
def primary_category(annotation: pd.DataFrame, loci: pd.Index) -> pd.Series:
    """One tile per protein: its first (most relevant) COG letter, '-' if none."""
    letters = annotation.set_index("b_number")["COG_Categories"].reindex(loci).fillna("")
    return letters.map(lambda s: s[0] if s else "-")


def condition_profiles(X: pd.DataFrame, obs: pd.DataFrame, cfg: Config):
    """Protein mass fractions per condition (replicates averaged, each row sums to 1)."""
    from .allocation import aggregate_conditions
    profiles, info = aggregate_conditions(X.fillna(0.0), obs, cfg)
    return profiles.div(profiles.sum(axis=1), axis=0), info


def select_conditions(info: pd.DataFrame, cfg: Config) -> list[str]:
    """Pick the two conditions shown in the proteomap.

    'auto' = slowest and fastest condition of the dataset with the widest
    growth-rate range, so both maps come from one study (no platform offset).
    """
    spec = cfg.get("proteomap", {}).get("conditions", "auto")
    if spec != "auto":
        missing = [c for c in spec if c not in info.index]
        if missing:
            raise KeyError(f"proteomap.conditions not found: {missing}. Valid names: {list(info.index)}")
        return list(spec)
    batch_key = cfg["allocation"]["batch_key"]
    spread = info.groupby(batch_key)[GROWTH].agg(["count", "min", "max"])
    spread = spread[spread["count"] >= 2]
    best = (spread["max"] - spread["min"]).idxmax()
    within = info[info[batch_key] == best]
    return [within[GROWTH].idxmin(), within[GROWTH].idxmax()]


def category_change(profiles: pd.DataFrame, category: pd.Series, conditions: list[str]) -> pd.DataFrame:
    """Percentage of the proteome in each category for the chosen conditions,
    and the change in percentage points between slow and fast growth.
    """
    shares = profiles.loc[conditions].T.groupby(category).sum() * 100
    shares.columns = ["slow", "fast"] if len(conditions) == 2 else conditions
    if len(conditions) == 2:
        shares["change_points"] = shares["fast"] - shares["slow"]
    return shares


# --------------------------------------------------------------------------- #
# drawing
# --------------------------------------------------------------------------- #
def _text_colour(hex_colour):
    """Dark or light label text, whichever is readable on the given fill colour."""
    r, g, b = to_rgb(hex_colour)
    return INK if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.55 else SURFACE


def _fit_text(ax, text, rect, max_size, colour, weight="bold"):
    """Write `text` centred in `rect` with the largest font (<= max_size) that fits; skip if < 5.5 pt."""
    x, y, w, h = rect
    box = ax.get_position()
    fig_w, fig_h = ax.figure.get_size_inches()
    width_pt, height_pt = w * box.width * fig_w * 72, h * box.height * fig_h * 72
    lines = text.split("\n")
    longest = max(len(line) for line in lines)
    size = min(max_size, 0.92 * width_pt / (0.6 * longest), 0.85 * height_pt / (1.25 * len(lines)))
    if size >= 5.5:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size, color=colour,
                fontweight=weight, linespacing=1.05)


def _layout(profile: pd.Series, category: pd.Series, group_order, category_order):
    """Nested rectangles in the unit square: {category: rect}, {protein: (rect, category)}."""
    cat_rects, prot_rects = {}, {}
    group_mass = {g: profile[category.isin(c)].sum() for g, c in SUPERGROUPS.items()}
    groups = [g for g in group_order if group_mass[g] > 0]
    for g, (gx, gy, gw, gh) in zip(groups, squarify([group_mass[g] for g in groups], 0, 0, 1, 1)):
        cats = [c for c in category_order if c in SUPERGROUPS[g] and profile[category == c].sum() > 0]
        masses = [profile[category == c].sum() for c in cats]
        for c, rect in zip(cats, squarify(masses, gx, gy, gw, gh)):
            cat_rects[c] = rect
            proteins = profile[(category == c) & (profile > 0)].sort_values(ascending=False)
            for p, sub in zip(proteins.index, squarify(proteins.to_numpy(), *rect)):
                prot_rects[p] = (sub, c)
    return cat_rects, prot_rects


def plot_proteomaps(profiles: pd.DataFrame, info: pd.DataFrame, category: pd.Series, conditions: list[str],
                    gene_names: pd.Series, change: pd.DataFrame, significant: pd.Series | None = None,
                    path=None, label_min_share: float = 0.006):
    """Top row: categories. Bottom row: proteins (same layout). Right: change per category."""
    selected = profiles.loc[conditions]
    mean_share = selected.mean().groupby(category).sum()
    group_order = sorted(SUPERGROUPS, key=lambda g: -mean_share.reindex(SUPERGROUPS[g]).fillna(0).sum())
    category_order = list(mean_share.sort_values(ascending=False).index)   # same order in every panel

    n = len(conditions)
    fig = plt.figure(figsize=(5.2 * n + 6.2, 10.4))
    grid = fig.add_gridspec(2, n + 2, width_ratios=[1] * n + [0.42, 0.75], wspace=0.06, hspace=0.12)
    for column, condition in enumerate(conditions):
        cat_rects, prot_rects = _layout(selected.loc[condition], category, group_order, category_order)
        top, bottom = fig.add_subplot(grid[0, column]), fig.add_subplot(grid[1, column])
        for ax in (top, bottom):
            ax.set_xlim(0, 1), ax.set_ylim(0, 1), ax.set_aspect("equal"), ax.axis("off")
        for c, (x, y, w, h) in cat_rects.items():
            top.add_patch(plt.Rectangle((x, y), w, h, facecolor=COLOURS.get(c, "#9c9b95"), edgecolor=SURFACE, linewidth=2))
            share = selected.loc[condition][category == c].sum()
            name = SHORT.get(c, c)
            if len(name) > 14 and " " in name:                     # wrap long names on two lines
                cut = name.rfind(" ", 0, len(name) // 2 + 4)
                name = name[:cut] + "\n" + name[cut + 1:] if cut > 0 else name
            _fit_text(top, f"{name}\n{share:.0%}", (x, y, w, h), 16, _text_colour(COLOURS.get(c, "#9c9b95")))
        for p, ((x, y, w, h), c) in prot_rects.items():
            bottom.add_patch(plt.Rectangle((x, y), w, h, facecolor=COLOURS.get(c, "#9c9b95"), edgecolor=SURFACE, linewidth=0.35))
            if w * h >= label_min_share:
                _fit_text(bottom, str(gene_names.get(p, p)), (x, y, w, h), 11, _text_colour(COLOURS.get(c, "#9c9b95")))
        for c, (x, y, w, h) in cat_rects.items():
            bottom.add_patch(plt.Rectangle((x, y), w, h, facecolor="none", edgecolor=SURFACE, linewidth=2))
        name = condition.split("|", 1)[-1]
        top.set_title(f"{name}  (μ = {info.loc[condition, GROWTH]:.2f} h$^{{-1}}$)\n{info.loc[condition, 'Dataset_Source']}",
                      fontsize=12, loc="center")

    ax = fig.add_subplot(grid[:, n + 1])                       # column n is a spacer for the bar labels
    data = change.reindex(category_order).dropna()
    data = data[(data[["slow", "fast"]].max(axis=1) >= 0.5)].sort_values("change_points")
    colours = [UP if v > 0 else DOWN for v in data["change_points"]]
    bars = ax.barh([SHORT.get(c, c) for c in data.index], data["change_points"], color=colours, height=0.65)
    marks = [("*" if significant is not None and bool(significant.get(c, False)) else "") for c in data.index]
    ax.bar_label(bars, labels=[f"{v:+.1f}{m}" for v, m in zip(data["change_points"], marks)], padding=3, fontsize=8)
    ax.axvline(0, color="#52514e", linewidth=1)
    ax.set_xlabel("Change in share, fast − slow (percentage points)")
    ax.set_title("Change per COG category\n(* = significant growth slope, notebook 03)", fontsize=10)
    ax.grid(axis="y", visible=False)
    low, high = ax.get_xlim()
    ax.set_xlim(low - 0.15 * (high - low), high + 0.15 * (high - low))
    fig.suptitle("Proteome allocation in the two most contrasting conditions of one batch "
                 "(area = share of proteome mass)", fontsize=13, fontweight="bold", y=0.99)
    return _save(fig, path)
