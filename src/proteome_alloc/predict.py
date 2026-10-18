"""Step 5 - prediction: can COG allocation predict growth rate in an unseen study?

Leakage control
---------------
* The whole model (scaling + PLS) is one scikit-learn Pipeline, so scaling
  parameters are learned from training conditions only, in every fold.
* The split matches the intended use ("a new laboratory"): leave-one-DATASET-
  out. Conditions of one study never sit on both sides of a split, and
  replicates were averaged before, so related observations stay together.
* Nested design: the inner leave-one-dataset-out loop chooses the number of
  PLS components using outer-training data only; the outer fold evaluates the
  entire model-building procedure. Every trial is recorded.
* Features are within-sample mass fractions computed from each sample alone;
  the COG annotation is external. No step uses the growth rate of a held-out
  condition.
* The score is reported next to a mean-predictor baseline, with a permutation
  null (growth rates shuffled within dataset, whole nested procedure repeated).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .config import Config
from .ingest import GROWTH


def make_model(n_components: int) -> Pipeline:
    """Scaler + PLS regression as one scikit-learn Pipeline.

    Because scaling is a pipeline step, its means and standard deviations are
    learned on the training folds only - this is what prevents preprocessing leakage.
    """
    return Pipeline([("scale", StandardScaler()),
                     ("pls", PLSRegression(n_components=n_components, scale=False))])


def _rmse(y, yhat) -> float:
    """Root-mean-square error between observed and predicted values."""
    return float(np.sqrt(np.mean((np.asarray(y) - np.asarray(yhat)) ** 2)))


def _select_components(X, y, groups, grid):
    """Inner loop: leave-one-dataset-out RMSE for each candidate. Returns (best, trials)."""
    trials = []
    for k in grid:
        errors, feasible = [], True
        for train, test in LeaveOneGroupOut().split(X, y, groups):
            if k > min(len(train) - 1, X.shape[1]):
                feasible = False
                break
            prediction = make_model(k).fit(X[train], y[train]).predict(X[test]).ravel()
            errors.append((y[test] - prediction) ** 2)
        trials.append({"n_components": k, "feasible": feasible,
                       "inner_rmse": float(np.sqrt(np.mean(np.concatenate(errors)))) if feasible else np.nan})
    table = pd.DataFrame(trials)
    if not table["feasible"].any():
        raise ValueError("No feasible number of components for this training set.")
    best = int(table.loc[table["inner_rmse"].idxmin(), "n_components"])   # ties -> fewest components
    return best, table


def nested_cv(X: np.ndarray, y: np.ndarray, groups: np.ndarray, grid) -> dict:
    """Outer leave-one-dataset-out evaluation with inner model selection."""
    prediction = np.full(len(y), np.nan)
    baseline = np.full(len(y), np.nan)
    chosen, trials, folds = np.zeros(len(y), int), [], []
    for train, test in LeaveOneGroupOut().split(X, y, groups):
        held_out = groups[test][0]
        if len(np.unique(groups[train])) >= 2:
            k, table = _select_components(X[train], y[train], groups[train], grid)
        else:                                   # cannot tune with one training dataset
            k, table = int(min(grid)), pd.DataFrame({"n_components": [min(grid)], "feasible": [True], "inner_rmse": [np.nan]})
        trials.append(table.assign(outer_fold_held_out=held_out, selected=lambda d: d["n_components"] == k))
        prediction[test] = make_model(k).fit(X[train], y[train]).predict(X[test]).ravel()
        baseline[test] = y[train].mean()
        chosen[test] = k
        folds.append({"held_out_dataset": held_out, "n_train": len(train), "n_test": len(test),
                      "n_components": k, "rmse_model": _rmse(y[test], prediction[test]),
                      "rmse_baseline": _rmse(y[test], baseline[test])})
    sse, sse_base = np.sum((y - prediction) ** 2), np.sum((y - baseline) ** 2)
    return {"prediction": prediction, "baseline": baseline, "n_components": chosen,
            "folds": pd.DataFrame(folds), "trials": pd.concat(trials, ignore_index=True),
            "rmse_model": _rmse(y, prediction), "rmse_baseline": _rmse(y, baseline),
            "mae_model": float(np.mean(np.abs(y - prediction))),
            "q2_vs_baseline": float(1.0 - sse / sse_base)}


def _permuted_rmse(seed, X, y, groups, grid) -> float:
    """Repeat the full nested cross-validation with growth rates shuffled within each dataset.

    One call = one draw of the permutation null distribution (run in parallel).
    """
    rng = np.random.default_rng(seed)
    shuffled = y.copy()
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        shuffled[idx] = y[rng.permutation(idx)]
    return nested_cv(X, shuffled, groups, grid)["rmse_model"]


def evaluate(Fc: pd.DataFrame, info: pd.DataFrame, cfg: Config) -> dict:
    """Full prediction experiment. Returns a dict of tables and numbers."""
    p = cfg["prediction"]
    features = [c for c in Fc.columns if c not in set(p["exclude_features"]) and Fc[c].std() > 0]
    X, y = Fc[features].to_numpy(float), info.loc[Fc.index, GROWTH].to_numpy(float)
    groups = info.loc[Fc.index, "Dataset_Source"].astype(str).to_numpy()
    grid = [int(k) for k in p["n_components_grid"]]
    if len(np.unique(groups)) < 3:
        raise ValueError("Nested leave-one-dataset-out needs at least 3 datasets in scope.")

    result = nested_cv(X, y, groups, grid)
    seeds = cfg.seed_sequence("prediction.permutation").spawn(int(p["n_permutations"]))
    null = np.array(Parallel(n_jobs=cfg["n_jobs"])(
        delayed(_permuted_rmse)(seed, X, y, groups, grid) for seed in seeds))
    result["null_rmse"] = null
    result["p_permutation"] = float((1 + np.sum(null <= result["rmse_model"])) / (1 + len(null)))

    result["predictions"] = pd.DataFrame({
        "Dataset_Source": groups, "observed": y, "predicted": result["prediction"],
        "baseline_prediction": result["baseline"], "n_components": result["n_components"],
        "outer_fold": groups}, index=Fc.index)

    # final model: procedure re-run on all development data, then jackknife over datasets
    k_final, final_trials = _select_components(X, y, groups, grid)
    final = make_model(k_final).fit(X, y)
    scale = final.named_steps["scale"].scale_
    coefficient = np.ravel(final.named_steps["pls"].coef_) / scale      # per unit of mass fraction
    left_out = {}
    for g in np.unique(groups):
        keep = groups != g
        if k_final > min(keep.sum() - 1, X.shape[1]):
            continue
        model = make_model(k_final).fit(X[keep], y[keep])
        left_out[g] = np.ravel(model.named_steps["pls"].coef_) / model.named_steps["scale"].scale_
    jack = pd.DataFrame(left_out, index=features)
    m = jack.shape[1]
    se = np.sqrt((m - 1) / m * ((jack.sub(coefficient, axis=0)) ** 2).sum(axis=1)) if m > 1 else np.nan
    result["coefficients"] = pd.DataFrame({
        "coefficient_h-1_per_unit_fraction": coefficient, "jackknife_se": se,
        "sign_stable": (np.sign(jack).eq(np.sign(coefficient), axis=0)).all(axis=1)},
        index=pd.Index(features, name="COG_Category")).join(jack.add_prefix("coef_without_"))
    result["final_n_components"], result["final_trials"] = k_final, final_trials
    result["features"], result["final_model"] = features, final
    result["n_configurations_tried"] = int(result["trials"]["feasible"].sum() + final_trials["feasible"].sum())
    return result
