"""Step 2 - analysis scope.

Growth rate changes for different REASONS (nutrient quality, temperature,
osmotic stress, antibiotics, ...), and the proteome responds differently to
each. Pooling them makes "the growth-rate effect" an average of opposing
trends. The primary analysis is therefore restricted to nutrient-quality
conditions: batch growth, reference temperature, growth rate set only by the
medium. Every sample gets an explicit Perturbation_Type and the rule that
assigned it, and the decision table is saved, so the exclusion is auditable.

The classification never uses protein abundances, so it cannot leak
information from the proteome into the model.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from .config import Config
from .ingest import GROWTH

NUTRIENT = "nutrient_quality"


def _text(value) -> str:
    """Missing values become an empty string (never the literal 'nan')."""
    return "" if value is None or (isinstance(value, float) and np.isnan(value)) or value is pd.NA \
        else str(value).strip()


def classify(obs: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """Add Perturbation_Type, Scope_Rule, In_Scope to the sample metadata."""
    s = cfg["scope"]
    rules = {kind: [re.compile(p, re.I) for p in patterns]
             for kind, patterns in s["perturbation_rules"].items()}
    special_medium = [re.compile(p, re.I) for p in s["special_medium_regex"]]
    out = obs.copy()
    kinds, reasons = [], []
    for _, row in out.iterrows():
        text = " ".join(_text(row.get(c)) for c in
                        ("Condition_ID", "Condition_Description", "Special_Procedure"))
        kind, reason = NUTRIENT, "no perturbation rule matched"
        for candidate, patterns in rules.items():
            hit = next((p.pattern for p in patterns if p.search(text)), None)
            if hit:
                kind, reason = candidate, f"description matches /{hit}/"
                break
        temperature = row.get("Temperature_C")
        if kind == NUTRIENT and pd.notna(temperature) and float(temperature) != float(s["reference_temperature_C"]):
            kind, reason = "temperature", f"Temperature_C = {temperature} != {s['reference_temperature_C']}"
        if kind == NUTRIENT and any(p.search(_text(row.get("Base_Medium"))) for p in special_medium):
            kind, reason = "other_perturbation", "non-standard base medium"
        procedure = _text(row.get("Special_Procedure"))
        if kind == NUTRIENT and s["special_procedure_is_perturbation"] and procedure not in ("", "-"):
            kind, reason = "other_perturbation", "special condition/procedure recorded"
        if kind == NUTRIENT and str(row.get("Cultivation_Method", "")) != s["cultivation_method"]:
            kind, reason = "non_batch", f"Cultivation_Method = {row.get('Cultivation_Method')}"
        kinds.append(kind)
        reasons.append(reason)
    out["Perturbation_Type"], out["Scope_Rule"] = kinds, reasons

    overrides = pd.read_csv(cfg.resolve(s["overrides"]), sep="\t", dtype=str, keep_default_na=False)
    unknown = set(overrides["Sample_ID"]) - set(out.index)
    if unknown:
        raise KeyError(f"scope_overrides.tsv refers to unknown Sample_ID(s): {sorted(unknown)}")
    for record in overrides.itertuples(index=False):
        out.loc[record.Sample_ID, "Perturbation_Type"] = record.Perturbation_Type
        out.loc[record.Sample_ID, "Scope_Rule"] = f"manual override: {record.Reason}"

    has_rate = np.isfinite(pd.to_numeric(out[GROWTH], errors="coerce"))
    strain_ok = out["Strain"].astype(str).str.contains(s["strain_contains"], na=False)
    out["In_Scope"] = (out["Perturbation_Type"] == NUTRIENT) & has_rate & strain_ok
    out.loc[~has_rate & (out["Perturbation_Type"] == NUTRIENT), "Scope_Rule"] = "no growth rate available"
    out.loc[~strain_ok & (out["Perturbation_Type"] == NUTRIENT), "Scope_Rule"] = "other strain"
    return out


def design_table(obs: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """Per dataset: samples, distinct conditions and growth-rate range in scope.
    This is the confounding check: a dataset covering a narrow growth-rate range
    cannot separate 'dataset' from 'growth rate'."""
    key = cfg["allocation"]["condition_key"]
    scoped = obs[obs["In_Scope"]]
    table = scoped.groupby("Dataset_Source").agg(
        n_samples=(GROWTH, "size"), n_conditions=(key, "nunique"),
        growth_min=(GROWTH, "min"), growth_median=(GROWTH, "median"), growth_max=(GROWTH, "max"))
    table.insert(0, "n_samples_total", obs.groupby("Dataset_Source").size())
    return table.reindex(obs["Dataset_Source"].unique()).fillna({"n_samples": 0, "n_conditions": 0})


def cross_dataset_offsets(X: pd.DataFrame, obs: pd.DataFrame) -> pd.DataFrame:
    """Unit sanity check. For every pair of datasets: median log10 ratio of the
    mean in-scope profiles over proteins quantified in both, and the number of
    such proteins. Values near 0 mean the datasets are on the same scale; a
    value of +-0.5 would be a ~3-fold unit mismatch."""
    scoped = obs[obs["In_Scope"]]
    profiles = X.loc[scoped.index].groupby(scoped["Dataset_Source"]).mean()
    names = list(profiles.index)
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            both = (profiles.loc[a] > 0) & (profiles.loc[b] > 0)
            ratio = np.log10(profiles.loc[a, both] / profiles.loc[b, both])
            rows.append({"dataset_A": a, "dataset_B": b, "n_shared_proteins": int(both.sum()),
                         "median_log10_ratio_A_over_B": float(ratio.median()),
                         "IQR_log10_ratio": float(ratio.quantile(0.75) - ratio.quantile(0.25)),
                         "spearman_r": float(profiles.loc[a, both].corr(profiles.loc[b, both], method="spearman"))})
    return pd.DataFrame(rows)


def save_scoped(cfg: Config, obs: pd.DataFrame):
    """Write the final metadata table (with scope decisions) used by all later notebooks."""
    path = cfg.path("processed_dir") / "sample_metadata_scoped.tsv"
    obs.to_csv(path, sep="\t")
    return path


def load_scoped(cfg: Config):
    """Matrix, scoped metadata and protein table written by notebooks 01 and 02."""
    from .ingest import load_processed
    X, _, var = load_processed(cfg)
    path = cfg.path("processed_dir") / "sample_metadata_scoped.tsv"
    if not path.exists():
        raise FileNotFoundError("Run notebook 02 first (sample_metadata_scoped.tsv is missing).")
    obs = pd.read_csv(path, sep="\t", index_col=0)
    if not obs.index.equals(X.index):
        raise AssertionError("Scoped metadata and matrix are out of sync; re-run notebooks 01 and 02.")
    return X, obs, var
