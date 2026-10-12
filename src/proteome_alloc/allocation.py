"""Step 3 - primary analysis: COG-category mass fractions versus growth rate.

Question answered: *which functional groups take a larger or smaller share of
the proteome as nutrient-limited growth gets faster, and by how much?*

Statistical design
------------------
* Unit of analysis = growth CONDITION, not MS run: technical/biological
  replicates of one condition are averaged first, otherwise nine glucose
  replicates would count as nine independent conditions (pseudo-replication).
* Model per category:  phi = a_dataset + b * mu.  Each dataset keeps its own
  intercept (absorbs additive platform/quantification offsets); the slope b
  is shared and is estimated ONLY from within-dataset contrasts, so a
  difference between datasets can never masquerade as a growth effect.
  phi is analysed on its natural (linear) scale because bacterial growth laws
  are linear in mass fraction and because fractions are additive.
* Inference without distributional assumptions and with few conditions:
  - permutation test: growth rates shuffled WITHIN each dataset;
  - bootstrap CI: conditions resampled WITHIN each dataset;
  - Benjamini-Hochberg FDR across categories;
  - leave-one-dataset-out (jackknife) slopes: is the result carried by one study?
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import false_discovery_control

from .annotation import category_weights
from .config import Config
from .ingest import GROWTH


# --------------------------------------------------------------------------- #
# protein -> category mass fractions
# --------------------------------------------------------------------------- #
def cog_fractions(X: pd.DataFrame, annotation: pd.DataFrame, cfg: Config):
    """Sum protein mass fractions per COG category.

    NaN (protein not in a dataset's table) and 0 (not detected) both contribute
    zero mass - the correct value for a mass balance, so nothing is imputed.
    Returns (fractions samples x categories, per-sample coverage table)."""
    W = category_weights(X.columns, annotation, cfg)
    mass = X.fillna(0.0)
    raw_sum = mass.sum(axis=1)
    F = mass @ W
    if cfg["allocation"]["renormalise"]:
        F = F.div(raw_sum, axis=0)
    unassigned = cfg["annotation"]["unassigned_label"]
    coverage = pd.DataFrame({
        "mass_fraction_sum_before_renormalisation": raw_sum,
        "n_proteins_detected": (mass > 0).sum(axis=1),
        "mass_fraction_with_COG": 1.0 - F[unassigned] / F.sum(axis=1),
        "n_detected_without_COG": ((mass > 0) & (W[unassigned] == 1.0).values).sum(axis=1),
    })
    return F, coverage


def aggregate_conditions(values: pd.DataFrame, obs: pd.DataFrame, cfg: Config):
    """Average replicates: one row per (batch, condition). Returns (values, obs)."""
    a = cfg["allocation"]
    keys = [a["batch_key"], a["condition_key"]]
    obs = obs.loc[values.index].assign(_dataset=lambda d: d["Dataset_Source"])
    if obs[keys].isna().any().any():
        raise ValueError(f"Missing {keys} for some in-scope samples.")
    mean = values.groupby([obs[k] for k in keys], sort=False).mean()
    info = obs.groupby(keys, sort=False).agg(
        n_replicates=(GROWTH, "size"), growth=(GROWTH, "mean"), growth_sd=(GROWTH, "std"),
        dataset=("_dataset", "first"),
        samples=("Original_Sample_ID", lambda s: ";".join(map(str, s))))
    info = info.rename(columns={"growth": GROWTH, "dataset": "Dataset_Source"}).loc[mean.index]
    labels = [f"{b}|{c}" for b, c in mean.index]
    mean.index = info.index = pd.Index(labels, name="Condition")
    info[a["batch_key"]] = [label.split("|", 1)[0] for label in labels]
    return mean, info


# --------------------------------------------------------------------------- #
# pooled within-dataset slope
# --------------------------------------------------------------------------- #
def _center_within(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    out = np.asarray(values, dtype=float).copy()
    for g in np.unique(groups):
        mask = groups == g
        out[mask] -= out[mask].mean(axis=0)
    return out


def within_slopes(Y: np.ndarray, mu: np.ndarray, groups: np.ndarray):
    """Common slope of every column of Y on mu with one intercept per group."""
    mu_c, Y_c = _center_within(mu, groups), _center_within(Y, groups)
    ss = float(mu_c @ mu_c)
    if ss <= 0:
        raise ValueError("No within-dataset variation in growth rate: slope not estimable.")
    slope = (mu_c @ Y_c) / ss
    residual = Y_c - np.outer(mu_c, slope)
    ss_within = (Y_c ** 2).sum(axis=0)
    r2 = np.where(ss_within > 0, 1.0 - (residual ** 2).sum(axis=0) / np.where(ss_within > 0, ss_within, 1), np.nan)
    return slope, r2


def _bootstrap_chunk(seed, n_draws, Y, mu, groups):
    rng = np.random.default_rng(seed)
    index_by_group = [np.where(groups == g)[0] for g in np.unique(groups)]
    draws = np.full((n_draws, Y.shape[1]), np.nan)
    for i in range(n_draws):
        pick = np.concatenate([rng.choice(idx, size=idx.size, replace=True) for idx in index_by_group])
        try:
            draws[i] = within_slopes(Y[pick], mu[pick], groups[pick])[0]
        except ValueError:
            continue                         # degenerate resample (no growth-rate spread)
    return draws


def fit_growth_laws(Fc: pd.DataFrame, info: pd.DataFrame, cfg: Config):
    """Slope of each COG category against growth rate. Returns (table, details)."""
    a = cfg["allocation"]
    mu = info[GROWTH].to_numpy(float)
    groups = info[a["batch_key"]].astype(str).to_numpy()
    Y = Fc.to_numpy(float)
    slope, r2 = within_slopes(Y, mu, groups)

    # ---- permutation test: growth rate shuffled within dataset ------------- #
    rng = cfg.rng("allocation.permutation")
    mu_c, Y_c = _center_within(mu, groups), _center_within(Y, groups)
    ss = float(mu_c @ mu_c)
    index_by_group = [np.where(groups == g)[0] for g in np.unique(groups)]
    n_perm = int(a["n_permutations"])
    exceed = np.zeros(Y.shape[1])
    for _ in range(n_perm):
        shuffled = mu_c.copy()
        for idx in index_by_group:
            shuffled[idx] = mu_c[rng.permutation(idx)]
        exceed += np.abs(shuffled @ Y_c / ss) >= np.abs(slope) - 1e-15
    p_perm = (1.0 + exceed) / (1.0 + n_perm)

    # ---- bootstrap CI (parallel, worker-count independent) ----------------- #
    n_boot, n_chunks = int(a["n_bootstrap"]), 8
    seeds = cfg.seed_sequence("allocation.bootstrap").spawn(n_chunks)
    sizes = [n_boot // n_chunks + (i < n_boot % n_chunks) for i in range(n_chunks)]
    chunks = Parallel(n_jobs=cfg["n_jobs"])(
        delayed(_bootstrap_chunk)(seed, size, Y, mu, groups) for seed, size in zip(seeds, sizes))
    boot = np.vstack(chunks)
    alpha = 1.0 - float(a["ci_level"])
    ci_low, ci_high = np.nanquantile(boot, [alpha / 2, 1 - alpha / 2], axis=0)

    # ---- per-dataset slopes and leave-one-dataset-out stability ------------- #
    per_dataset, leave_out = {}, {}
    for g in np.unique(groups):
        inside = groups == g
        if inside.sum() >= 2 and np.ptp(mu[inside]) > 0:
            per_dataset[g] = within_slopes(Y[inside], mu[inside], groups[inside])[0]
        try:
            leave_out[g] = within_slopes(Y[~inside], mu[~inside], groups[~inside])[0]
        except ValueError:
            pass
    per_dataset = pd.DataFrame(per_dataset, index=Fc.columns).add_prefix("slope_")
    leave_out = pd.DataFrame(leave_out, index=Fc.columns).add_prefix("slope_without_")

    table = pd.DataFrame({
        "mean_fraction": Fc.mean(axis=0),
        "slope_per_h-1": slope,
        "ci_low": ci_low, "ci_high": ci_high,
        "relative_slope": slope / Fc.mean(axis=0).replace(0, np.nan),
        "within_dataset_R2": r2,
        "p_permutation": p_perm,
    }, index=Fc.columns)
    tested = table["mean_fraction"] >= float(a["min_mean_fraction"])
    table["tested"] = tested
    table["q_BH"] = np.nan
    table.loc[tested, "q_BH"] = false_discovery_control(table.loc[tested, "p_permutation"], method="bh")
    table["significant"] = table["q_BH"] <= float(a["fdr_alpha"])
    signs = np.sign(leave_out.to_numpy())
    table["sign_stable_leave_one_dataset_out"] = (signs == np.sign(slope)[:, None]).all(axis=1)
    table = pd.concat([table, per_dataset, leave_out], axis=1)
    table.index.name = "COG_Category"
    details = {"n_conditions": int(len(mu)), "conditions_per_dataset": pd.Series(groups).value_counts().to_dict(),
               "n_permutations": n_perm, "n_bootstrap": n_boot,
               "bootstrap_failed_draws": int(np.isnan(boot[:, 0]).sum()),
               "smallest_attainable_p": 1.0 / (1.0 + n_perm)}
    return table.sort_values("slope_per_h-1", ascending=False), details
