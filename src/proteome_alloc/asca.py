"""Step 4 - secondary analysis: ASCA at protein level.

ANOVA-Simultaneous Component Analysis splits the centred data matrix into
additive, mutually orthogonal effect matrices

    X_c = X_dataset + X_growth + E

X_dataset is the projection onto the dataset indicators (entered FIRST), and
X_growth the projection onto growth rate AFTER the dataset space has been
removed (sequential sums of squares). Growth is therefore credited only with
variation that dataset membership cannot explain - necessary because the
design is unbalanced (datasets cover different growth-rate ranges).

Growth rate is used as a continuous covariate (1 degree of freedom). The
growth effect matrix then has rank 1, so its single component IS the vector of
per-protein within-dataset slopes; ASCA adds the multivariate significance
test and the variance partition, not extra dimensions. Quantile-binning the
growth rate was dropped: it discards information and its bin edges depend on
which samples are in the matrix.

ASCA is descriptive (no prediction); the only data-driven filter (detection
completeness) never looks at growth rate, so there is no label leakage.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, mannwhitneyu

from .config import Config
from .ingest import GROWTH


def prepare_protein_matrix(X: pd.DataFrame, obs: pd.DataFrame, cfg: Config):
    """In-scope samples -> filtered, log-transformed protein matrix (NOT yet centred).

    1. keep proteins quantified (> 0) in at least `min_detect_fraction_per_dataset`
       of the in-scope samples of EVERY dataset (logical AND). With the default
       1.0 the matrix is complete and nothing is imputed.
    2. only if the threshold is < 1: per-sample half-minimum for the remaining
       zeros (left-censored values; uses that sample alone, so no leakage).
    3. log transform (base from config), no offset.
    Returns (log matrix, report dict)."""
    p = cfg["asca"]
    threshold = float(p["min_detect_fraction_per_dataset"])
    detected = X.fillna(0.0) > 0
    per_dataset = detected.groupby(obs.loc[X.index, "Dataset_Source"]).mean()
    keep = (per_dataset >= threshold).all(axis=0)
    kept = X.loc[:, keep].fillna(0.0)
    n_imputed = int((kept <= 0).sum().sum())
    if n_imputed:
        half_min = X.where(X > 0).min(axis=1) * float(p["half_min_factor"])
        kept = kept.where(kept > 0, np.nan).apply(lambda col: col.fillna(half_min))
    logged = np.log(kept) / np.log(float(p["log_base"]))
    report = {"proteins_in_merged_matrix": int(X.shape[1]), "proteins_kept": int(keep.sum()),
              "threshold_per_dataset": threshold, "values_imputed_half_minimum": n_imputed,
              "share_imputed": n_imputed / float(kept.size) if kept.size else 0.0,
              "log_base": p["log_base"],
              "mass_fraction_covered_by_kept_proteins_min": float(X.loc[:, keep].sum(axis=1).min())}
    return logged, report


def _onehot(labels) -> np.ndarray:
    return pd.get_dummies(pd.Series(np.asarray(labels).astype(str)), dtype=float).to_numpy()


def _basis(D: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """Orthonormal basis of the column space of D (safe for rank-deficient D)."""
    U, s, _ = np.linalg.svd(D, full_matrices=False)
    if s.size == 0 or s[0] <= 0:
        return np.zeros((D.shape[0], 0))
    return U[:, : int((s > tol * s[0]).sum())]


@dataclass
class ASCAResult:
    variance_table: pd.DataFrame
    protein_table: pd.DataFrame
    scores: pd.DataFrame
    p_growth: float
    null_distribution: np.ndarray
    observed_ss_growth: float
    params: dict = field(default_factory=dict)
    p_interaction: float = np.nan
    null_interaction: np.ndarray | None = None
    observed_f_interaction: float = np.nan


def run_asca(L: pd.DataFrame, info: pd.DataFrame, cfg: Config, gene_names: pd.Series | None = None) -> ASCAResult:
    """L: conditions x proteins, log scale, complete. info: condition metadata.

    Model (sequential sums of squares, each term orthogonal to the ones before it):
        X_c = X_batch + X_growth + X_batch:growth + E
    X_batch        : one mean level per batch
    X_growth       : one growth slope shared by all batches (after batch)
    X_batch:growth : how much each batch's own slope deviates from the shared slope
    E              : residual
    """
    p, batch_key = cfg["asca"], cfg["allocation"]["batch_key"]
    include_interaction = bool(p.get("include_interaction", True))
    if not np.isfinite(L.to_numpy()).all():
        raise ValueError("Matrix contains NaN/inf.")
    Xc = (L - L.mean(axis=0)).to_numpy(float)          # centring is the LAST preprocessing step
    batch = info.loc[L.index, batch_key].astype(str).to_numpy()
    mu = info.loc[L.index, GROWTH].to_numpy(float)
    n = Xc.shape[0]

    # ---- batch, then shared growth slope ------------------------------------
    Qb = _basis(_onehot(batch))
    remove_batch = np.eye(n) - Qb @ Qb.T
    growth_design = remove_batch @ (mu - mu.mean())[:, None]
    Qg = _basis(growth_design)
    if Qg.shape[1] == 0:
        raise ValueError("Growth rate is fully confounded with batch: no estimable growth effect.")
    X_batch, X_growth = Qb @ (Qb.T @ Xc), Qg @ (Qg.T @ Xc)

    # ---- batch x growth interaction: batch-specific slopes beyond the shared one
    # Full model space = one intercept AND one slope per batch. The interaction space is
    # what that adds to (batch + shared slope). Batches with a single condition, or with
    # no spread in growth rate, have no slope of their own and add nothing here.
    Q_reduced = np.column_stack([Qb, Qg])
    if include_interaction:
        per_batch_slopes = _onehot(batch) * (mu - mu.mean())[:, None]
        Qi = _basis(per_batch_slopes - Q_reduced @ (Q_reduced.T @ per_batch_slopes))
    else:
        Qi = np.zeros((n, 0))
    X_inter = Qi @ (Qi.T @ Xc)

    E = Xc - X_batch - X_growth - X_inter
    ss = lambda M: float(np.sum(M * M))
    total = ss(Xc)
    assert np.allclose(Xc, X_batch + X_growth + X_inter + E)
    for A, B in ((X_batch, X_growth), (X_batch, X_inter), (X_growth, X_inter), (X_growth, E), (X_inter, E)):
        assert abs(np.sum(A * B)) < 1e-8 * total, "effect matrices are not orthogonal"

    # ---- test 1, shared growth slope: growth shuffled within batch -----------
    rng = cfg.rng("asca.permutation")
    gram = Xc @ Xc.T
    index_by_batch = [np.where(batch == b)[0] for b in np.unique(batch)]

    def statistic(values: np.ndarray) -> float:
        q = _basis(remove_batch @ (values - values.mean())[:, None])
        return float(np.sum(q * (gram @ q)))

    observed = statistic(mu)
    n_perm = int(p["n_permutations"])
    null = np.empty(n_perm)
    for i in range(n_perm):
        shuffled = mu.copy()
        for idx in index_by_batch:
            shuffled[idx] = mu[rng.permutation(idx)]
        null[i] = statistic(shuffled)
    p_growth = float((1 + np.sum(null >= observed - 1e-12 * total)) / (1 + n_perm))

    # ---- test 2, interaction: Freedman-Lane permutation with a pseudo-F statistic
    # Under H0 (all batches share one slope) the residuals of the batch + growth model are
    # shuffled within batch and added back to that model's fit; the full model is refitted.
    # The statistic is a pseudo-F, F = (SS_interaction / df) / (SS_residual / df_residual),
    # pooled over proteins. A ratio is needed: the raw SS is not exchangeable under residual
    # permutation (shuffling moves part of the growth direction into the interaction space)
    # and gives far too many false positives. Calibration was checked by simulation.
    df_inter = Qi.shape[1]
    df_resid = n - Qb.shape[1] - Qg.shape[1] - df_inter
    if df_inter > 0 and df_resid > 0:
        R = Xc - X_batch - X_growth
        gram_r = R @ R.T
        ss_r = float(np.trace(gram_r))

        def pseudo_f(g: np.ndarray) -> float:
            ss_i = float(np.sum(Qi * (g @ Qi)))                      # ||Qi' R*||^2
            ss_red = float(np.sum(Q_reduced * (g @ Q_reduced)))       # part refitted by batch + growth
            ss_e = ss_r - ss_red - ss_i
            return (ss_i / df_inter) / (ss_e / df_resid) if ss_e > 1e-12 * total else np.inf

        observed_inter = pseudo_f(gram_r)
        rng_i = cfg.rng("asca.permutation.interaction")
        null_inter = np.empty(n_perm)
        for i in range(n_perm):
            order = np.arange(n)
            for idx in index_by_batch:
                order[idx] = idx[rng_i.permutation(idx.size)]
            null_inter[i] = pseudo_f(gram_r[np.ix_(order, order)])
        p_inter = float((1 + np.sum(null_inter >= observed_inter * (1 - 1e-12))) / (1 + n_perm))
    else:
        observed_inter, null_inter, p_inter = np.nan, None, np.nan

    df_batch, df_growth = Qb.shape[1] - 1, Qg.shape[1]
    variance_table = pd.DataFrame({
        "SS": [ss(X_batch), ss(X_growth), ss(X_inter), ss(E), total],
        "SS_fraction": [ss(X_batch) / total, ss(X_growth) / total, ss(X_inter) / total, ss(E) / total, 1.0],
        "df": [df_batch, df_growth, df_inter, n - 1 - df_batch - df_growth - df_inter, n - 1],
        "p_permutation": [np.nan, p_growth, p_inter, np.nan, np.nan],
    }, index=["Dataset", "Growth rate | dataset", "Dataset x growth rate", "Residual", "Total (centred)"])

    # rank-1 shared growth effect: loading = unit vector along the per-protein slopes
    mu_perp = growth_design[:, 0]
    slope = (mu_perp @ Xc) / float(mu_perp @ mu_perp)            # log units per h-1
    loading = slope / np.linalg.norm(slope)
    scores = pd.DataFrame({"effect_score": X_growth @ loading,
                           "effect_plus_residual_score": (X_growth + E) @ loading,
                           GROWTH: mu, batch_key: batch,
                           "Dataset_Source": info.loc[L.index, "Dataset_Source"].to_numpy()}, index=L.index)
    ss_protein = (Xc ** 2).sum(axis=0)
    share = lambda M: np.divide((M ** 2).sum(axis=0), ss_protein, out=np.zeros_like(ss_protein), where=ss_protein > 0)
    protein_table = pd.DataFrame({
        "gene": (gene_names.reindex(L.columns).to_numpy() if gene_names is not None else L.columns),
        "slope_log_per_h-1": slope,
        "loading": loading,
        "growth_share_of_protein_SS": share(X_growth),
        "interaction_share_of_protein_SS": share(X_inter),
        "dataset_share_of_protein_SS": share(X_batch),
    }, index=L.columns)
    protein_table.index.name = "b_number"
    params = {"n_conditions": n, "n_proteins": int(L.shape[1]), "n_permutations": n_perm,
              "batch_key": batch_key, "growth_coding": "continuous, centred, 1 df",
              "interaction_included": include_interaction, "interaction_df": df_inter,
              "p_interaction": p_inter, "interaction_pseudo_F": observed_inter,
              "interaction_test": "Freedman-Lane, pseudo-F, residuals permuted within batch"}
    return ASCAResult(variance_table, protein_table.sort_values("slope_log_per_h-1", ascending=False),
                      scores, p_growth, null, observed, params, p_inter, null_inter, observed_inter)


def cog_summary(protein_table: pd.DataFrame, weights: pd.DataFrame, cfg: Config, min_size: int = 5) -> pd.DataFrame:
    """Are the proteins of a COG category shifted up/down with growth relative to
    all other proteins? Mann-Whitney U on the per-protein slopes, BH-adjusted.
    Exploratory: proteins are not independent observations."""
    slopes = protein_table["slope_log_per_h-1"]
    member = weights.reindex(slopes.index).fillna(0.0) > 0
    rows = []
    for category in member.columns:
        inside, outside = slopes[member[category]], slopes[~member[category]]
        if len(inside) < min_size or len(outside) < min_size:
            continue
        stat, p = mannwhitneyu(inside, outside, alternative="two-sided")
        rows.append({"COG_Category": category, "n_proteins": len(inside),
                     "median_slope": inside.median(), "median_slope_other": outside.median(),
                     "rank_biserial": 2 * stat / (len(inside) * len(outside)) - 1, "p_mannwhitney": p})
    table = pd.DataFrame(rows).set_index("COG_Category")
    table["q_BH"] = false_discovery_control(table["p_mannwhitney"], method="bh")
    return table.sort_values("median_slope", ascending=False)
