"""Step 1 - ingestion and harmonisation.

One loader per publication. Every loader returns the SAME structure
(:class:`Dataset`): a samples x b-number matrix of protein mass fractions, a
sample-metadata table with a fixed schema, and a log of every row/sample that
was dropped and why. Nothing is typed into the code: sample metadata that the
publications only give in text live in ``metadata/*.tsv`` (with a source column)
and all file/sheet/column names live in ``config.yml``.

Value semantics of the merged matrix
------------------------------------
* value > 0 : protein quantified in that sample
* value = 0 : protein listed in that dataset's table but not detected in the sample
* NaN       : protein absent from that dataset's table altogether
For a mass-fraction (allocation) analysis an undetected protein contributes
zero mass, so no imputation is needed at the functional-group level.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .config import Config, sha256

GROWTH = "Specific_Growth_Rate_h-1"

#: fixed schema of the harmonised sample metadata (order = column order on disk)
OBS_COLUMNS = [
    "Original_Sample_ID", "Dataset_Source", "Strain", "Cultivation_Method",
    "Base_Medium", "Carbon_Source", "Supplement", "Temperature_C",
    "Condition_ID", "Condition_Description", "Special_Procedure", "MS_Batch",
    GROWTH, "Growth_Rate_Source", "Original_Unit", "Source_File", "Notes",
]


@dataclass
class Dataset:
    """One loaded study in the common schema: samples x b-numbers mass fractions (X),
    sample metadata (obs), gene names, a log of every row/sample dropped, and file checksums.
    """
    name: str
    X: pd.DataFrame                     # samples x b-numbers, mass fractions
    obs: pd.DataFrame                   # index = Sample_ID, columns = OBS_COLUMNS
    genes: pd.Series                    # b-number -> gene name
    log: dict = field(default_factory=dict)
    files: dict = field(default_factory=dict)   # file name -> sha256


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _raw(cfg: Config, name: str) -> Path:
    """Path of a supplementary file in data/raw/, with a clear error if it is missing."""
    path = cfg.path("raw_dir") / name
    if not path.exists():
        raise FileNotFoundError(
            f"Raw file not found: {path}\nPlace the supplementary file in "
            f"{cfg['paths']['raw_dir']}/ (see data/raw/README.md) or fix the name in config.yml.")
    return path


def _col(df: pd.DataFrame, wanted: str) -> str:
    """Exact column name for `wanted`, tolerant to case and surrounding blanks."""
    lookup = {str(c).strip().lower(): c for c in df.columns}
    key = wanted.strip().lower()
    if key not in lookup:
        raise KeyError(f"Column '{wanted}' not found. Available: {[str(c) for c in df.columns][:40]}")
    return lookup[key]


def _read_metadata(cfg: Config, rel_path: str) -> pd.DataFrame:
    """Read a curated sample table from metadata/ and index it by Original_Sample_ID.

    Stops on duplicated sample IDs, because a duplicate would silently attach
    the wrong growth rate to a proteome.
    """
    meta = pd.read_csv(cfg.resolve(rel_path), sep="\t", dtype=str, keep_default_na=False)
    if meta["Original_Sample_ID"].duplicated().any():
        raise ValueError(f"Duplicated Original_Sample_ID in {rel_path}")
    return meta.set_index("Original_Sample_ID")


def _finalise(name: str, wide: pd.DataFrame, obs: pd.DataFrame, genes: pd.Series,
              cfg: Config, log: dict, files: dict) -> Dataset:
    """Shared clean-up: valid b-numbers, duplicate loci summed, numeric checks,
    unique sample IDs, fixed metadata schema. `wide` is loci (index) x samples."""
    locus_re = re.compile(cfg["ingest"]["locus_regex"])
    wide = wide.copy()
    wide.index = wide.index.astype(str).str.strip()
    valid = np.array([bool(locus_re.match(s)) for s in wide.index], dtype=bool)
    log["rows_invalid_locus_dropped"] = int((~valid).sum())
    log["mass_fraction_in_invalid_locus_rows_max"] = float(
        wide.loc[~valid].apply(pd.to_numeric, errors="coerce").sum().max() if (~valid).any() else 0.0)
    wide = wide.loc[valid].apply(pd.to_numeric, errors="coerce")

    log["cells_non_numeric_set_to_zero"] = int(wide.isna().sum().sum())
    wide = wide.fillna(0.0)             # listed in the table, no value -> not detected
    if (wide < 0).any().any():
        raise ValueError(f"{name}: negative abundances found.")

    log["duplicated_locus_rows_summed"] = int(wide.index.duplicated().sum())
    wide = wide.groupby(level=0).sum()  # isoforms / split entries of one locus are additive in mass

    genes = genes.copy()
    genes.index = genes.index.astype(str).str.strip()
    genes = genes[~genes.index.duplicated()].reindex(wide.index)

    X = wide.T
    X.index = [f"{name}:{s}" for s in X.index]
    X.index.name = "Sample_ID"
    X.columns.name = "b_number"

    obs = obs.copy()
    obs["Original_Sample_ID"] = obs.index.astype(str)
    obs.index = [f"{name}:{s}" for s in obs.index]
    obs.index.name = "Sample_ID"
    obs["Dataset_Source"] = name
    for column in OBS_COLUMNS:
        if column not in obs.columns:
            obs[column] = np.nan
    obs = obs[OBS_COLUMNS].loc[X.index]
    obs[GROWTH] = pd.to_numeric(obs[GROWTH], errors="coerce")
    obs["Temperature_C"] = pd.to_numeric(obs["Temperature_C"], errors="coerce")

    sums = X.sum(axis=1)
    lo, hi = cfg["ingest"]["sample_sum_bounds"]
    bad = sums[(sums < lo) | (sums > hi)]
    log["sample_sum_min"], log["sample_sum_max"] = float(sums.min()), float(sums.max())
    log["samples_outside_sum_bounds"] = bad.round(4).to_dict()
    if len(bad):
        warnings.warn(f"{name}: {len(bad)} sample(s) do not sum to ~1 (bounds {lo}-{hi}): "
                      f"{bad.round(3).to_dict()}. Check units / missing sheets.")
    log["n_samples"], log["n_loci"] = int(X.shape[0]), int(X.shape[1])
    log["detected_per_sample_median"] = float((X > 0).sum(axis=1).median())
    return Dataset(name, X, obs, genes, log, files)


def _missing_samples(cfg: Config, name: str, expected, found, log: dict) -> None:
    """Report samples that are listed in the metadata but have no data column.

    They are logged and warned about (or raise an error when
    `ingest.strict_missing_samples` is true) - never skipped silently.
    """
    missing = [s for s in expected if s not in set(found)]
    log["samples_in_metadata_not_in_files"] = missing
    if missing:
        message = (f"{name}: {len(missing)} sample(s) listed in the metadata table have no data "
                   f"column in the configured files/sheets: {missing}")
        if cfg["ingest"]["strict_missing_samples"]:
            raise KeyError(message)
        warnings.warn(message)


# --------------------------------------------------------------------------- #
# Schmidt et al. 2016 - copies per cell -> mass fraction of the measured proteome
# --------------------------------------------------------------------------- #
def load_schmidt(cfg: Config, annotation: pd.DataFrame) -> Dataset:
    """Copies/cell are converted with  phi_i = N_i * MW_i / sum_j(N_j * MW_j).

    This is the same quantity the other three studies report (fraction of total
    protein mass), needs no cell volume or dry-weight assumption, and its
    denominator includes every quantified protein - also those that cannot be
    mapped to a b-number - so the fractions are not inflated by the mapping.
    UniProt accession -> b-number comes from the cached annotation table."""
    c = cfg["datasets"]["schmidt_2016"]
    path = _raw(cfg, c["file"])
    log, files = {}, {path.name: sha256(path)}
    table = pd.read_excel(path, sheet_name=c["abundance_sheet"], skiprows=c["skiprows"], header=0)
    growth = pd.read_excel(path, sheet_name=c["growth_sheet"], skiprows=c["skiprows"], header=0)

    acc, gene, mw = _col(table, c["accession_col"]), _col(table, c["gene_col"]), _col(table, c["mw_col"])
    prefix = c["copies_prefix"].format(strain=c["strain"])
    sample_cols = [col for col in table.columns if str(col).startswith(prefix)]
    if not sample_cols:
        candidates = [str(col) for col in table.columns if "copies" in str(col).lower()]
        raise KeyError(f"No column starts with '{prefix}'. Copies columns available: {candidates}")
    samples = [str(col)[len(prefix):] for col in sample_cols]

    table = table.dropna(subset=[acc]).copy()
    table[acc] = table[acc].astype(str).str.strip()
    copies = table[sample_cols].apply(pd.to_numeric, errors="coerce")
    copies.columns = samples

    # molecular weight: the value published with the data; UniProt mass only as a documented fallback
    acc_to_mass = annotation.dropna(subset=["UniProt_Accession"]).drop_duplicates("UniProt_Accession") \
        .set_index("UniProt_Accession")["Mass_Da"]
    weight = pd.to_numeric(table[mw], errors="coerce")
    fallback = weight.isna() & table[acc].map(acc_to_mass).notna()
    weight = weight.where(~fallback, table[acc].map(acc_to_mass))
    log["rows_mw_filled_from_uniprot"] = int(fallback.sum())
    usable = weight.notna() & copies.notna().any(axis=1)
    log["rows_dropped_no_mw_or_no_copies"] = int((~usable).sum())

    mass = copies.loc[usable].fillna(0.0).mul(weight.loc[usable], axis=0)
    fractions = mass / mass.sum(axis=0)          # denominator = ALL quantified protein mass

    acc_to_locus = annotation.dropna(subset=["UniProt_Accession"]).drop_duplicates("UniProt_Accession") \
        .set_index("UniProt_Accession")["b_number"]
    locus = table.loc[usable, acc].map(acc_to_locus)
    unmapped = locus.isna()
    log["rows_unmapped_accession_dropped"] = int(unmapped.sum())
    log["mass_fraction_lost_unmapped_per_sample"] = fractions.loc[unmapped].sum().round(5).to_dict()
    fractions = fractions.loc[~unmapped]
    fractions.index = locus.loc[~unmapped].values
    genes = pd.Series(table.loc[usable, gene].loc[~unmapped].values, index=fractions.index)

    meta = _read_metadata(cfg, c["metadata"])
    not_described = [s for s in samples if s not in meta.index]
    if not_described:
        raise KeyError(f"schmidt_2016: add a row to {c['metadata']} for sample(s) {not_described}; "
                       "metadata are never guessed by the code.")
    obs = meta.loc[samples].copy()
    g_strain, g_cond, g_rate = (_col(growth, c[k]) for k in
                                ("growth_strain_col", "growth_condition_col", "growth_rate_col"))
    rates = growth[growth[g_strain].astype(str).str.strip() == c["strain"]] \
        .assign(_cond=lambda d: d[g_cond].astype(str).str.strip()).set_index("_cond")[g_rate]
    obs[GROWTH] = [rates.get(s, np.nan) for s in samples]
    obs["Growth_Rate_Source"] = f"{path.name} / {c['growth_sheet']} / '{c['growth_rate_col']}'"
    obs["Strain"] = c["strain"]
    obs["MS_Batch"] = "schmidt_2016"
    obs["Original_Unit"] = c["original_unit"]
    obs["Source_File"] = f"{path.name} / {c['abundance_sheet']}"
    return _finalise("schmidt_2016", fractions, obs, genes, cfg, log, files)


# --------------------------------------------------------------------------- #
# Mori et al. 2021 - absolute mass fractions, metadata from the samples table
# --------------------------------------------------------------------------- #
def load_mori(cfg: Config) -> Dataset:
    """Load Mori et al. 2021 (EV8 absolute mass fractions + EV2 sample annotation).

    Only NCM3722 samples are kept. Growth rate is computed from the reported
    doubling time as ln(2) / doubling time.
    """
    c = cfg["datasets"]["mori_2021"]
    a_path, s_path = _raw(cfg, c["abundance_file"]), _raw(cfg, c["samples_file"])
    log, files = {}, {a_path.name: sha256(a_path), s_path.name: sha256(s_path)}
    with pd.ExcelFile(a_path) as book:            # recorded so that missing sample series are easy to spot
        log["abundance_file_sheets_available"] = list(book.sheet_names)
    with pd.ExcelFile(s_path) as book:
        log["samples_file_sheets_available"] = list(book.sheet_names)

    samples = pd.concat([pd.read_excel(s_path, sheet_name=s) for s in c["samples_sheets"]],
                        ignore_index=True)
    sid, strain = _col(samples, c["sample_id_col"]), _col(samples, c["strain_col"])
    samples = samples.dropna(subset=[sid]).copy()
    samples[sid] = samples[sid].astype(str).str.strip()
    keep = samples[strain].astype(str).str.contains(c["strain_contains"], na=False)
    log["samples_other_strain_dropped"] = int((~keep).sum())
    samples = samples.loc[keep].drop_duplicates(sid).set_index(sid)

    blocks, genes = [], pd.Series(dtype=object)
    for sheet in c["abundance_sheets"]:
        df = pd.read_excel(a_path, sheet_name=sheet)
        locus, gene = _col(df, c["locus_col"]), _col(df, c["gene_col"])
        df = df.dropna(subset=[locus]).copy()
        df[locus] = df[locus].astype(str).str.strip()
        genes = pd.concat([genes, df.drop_duplicates(locus).set_index(locus)[gene]])
        cols = [col for col in df.columns if str(col).strip() in samples.index]
        blocks.append(df.set_index(locus)[cols].rename(columns=lambda s: str(s).strip())
                      .apply(pd.to_numeric, errors="coerce").groupby(level=0).sum())
    wide = pd.concat(blocks, axis=1)
    if wide.columns.duplicated().any():
        raise ValueError("mori_2021: the same sample occurs in more than one abundance sheet.")
    _missing_samples(cfg, "mori_2021", samples.index, wide.columns, log)

    def text(column: str) -> pd.Series:
        """Text column of the sample sheet, stripped; an empty column if the sheet does not have it."""
        try:
            return samples[_col(samples, column)].astype("string").fillna("").str.strip()
        except KeyError:
            return pd.Series("", index=samples.index, dtype="string")

    doubling = pd.to_numeric(samples[_col(samples, c["doubling_time_col"])], errors="coerce")
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.log(2) / (doubling / 60.0)       # mu = ln2 / t_d, t_d in hours
    rate = rate.where(np.isfinite(rate) & (doubling > 0))
    log["samples_without_growth_rate"] = rate.index[rate.isna()].tolist()

    obs = pd.DataFrame(index=samples.index)
    obs["Strain"] = samples[strain].astype(str).str.strip()
    obs["Cultivation_Method"] = "Batch"
    obs["Base_Medium"] = text("Base medium")
    obs["Carbon_Source"] = text("Carbon Source")
    obs["Supplement"] = text("Supplement")
    # Temperature is not a column of the EV2 samples table. The growth temperatures stated in the
    # Mori 2021 Methods are set in config.yml: default_temperature_C for all samples, and
    # temperature_shift_C for samples whose description names a temperature shift ("Hi Temp").
    default_temperature = c.get("default_temperature_C")
    shifted_temperature = c.get("temperature_shift_C")
    temperature_rule = re.compile("|".join(cfg["scope"]["perturbation_rules"]["temperature"]), re.I)
    description = (text("Short Description") + " " + text("Description")).astype(str)
    obs["Temperature_C"] = [
        (float(shifted_temperature) if shifted_temperature is not None else np.nan)
        if temperature_rule.search(d)
        else (float(default_temperature) if default_temperature is not None else np.nan)
        for d in description]
    obs["Condition_ID"] = text("Short Description")
    obs["Condition_Description"] = (text("Short Description") + " | " + text("Description")).str.strip(" |")
    obs["Special_Procedure"] = text("Special condition/procedure")
    obs["MS_Batch"] = "mori_2021"
    obs[GROWTH] = rate
    obs["Growth_Rate_Source"] = f"{s_path.name} / '{c['doubling_time_col']}' (mu = ln2 / t_d)"
    obs["Original_Unit"] = c["original_unit"]
    obs["Source_File"] = f"{a_path.name} / {', '.join(c['abundance_sheets'])}"
    obs["Notes"] = "Nitrogen source: " + text("Nitrogen Source")
    return _finalise("mori_2021", wide, obs.loc[wide.columns], genes, cfg, log, files)


# --------------------------------------------------------------------------- #
# Wu et al. 2023 - three measurement groups (= MS batches), metadata from TSV
# --------------------------------------------------------------------------- #
def load_wu(cfg: Config) -> Dataset:
    """Load Wu et al. 2023: one sheet per MS batch (Group #1-#3) of protein mass fractions.

    Sample metadata and growth rates come from metadata/wu_2023_samples.tsv;
    the MS batch is kept so batch effects can be modelled later.
    """
    c = cfg["datasets"]["wu_2023"]
    path = _raw(cfg, c["file"])
    log, files = {}, {path.name: sha256(path)}
    meta = _read_metadata(cfg, c["metadata"])

    blocks, genes = [], pd.Series(dtype=object)
    for batch, sheet in c["sheets"].items():
        df = pd.read_excel(path, sheet_name=sheet)
        locus = _col(df, c["locus_col"])
        df = df.dropna(subset=[locus]).copy()
        df[locus] = df[locus].astype(str).str.strip()
        try:
            genes = pd.concat([genes, df.drop_duplicates(locus).set_index(locus)[_col(df, c["gene_col"])]])
        except KeyError:
            pass                                    # gene names only needed once
        wanted = meta.index[meta["MS_Batch"] == batch]
        cols = [col for col in df.columns if str(col).strip() in set(wanted)]
        blocks.append(df.set_index(locus)[cols].rename(columns=lambda s: str(s).strip())
                      .apply(pd.to_numeric, errors="coerce").groupby(level=0).sum())
    wide = pd.concat(blocks, axis=1)
    _missing_samples(cfg, "wu_2023", meta.index, wide.columns, log)

    obs = meta.loc[wide.columns].copy()
    obs["Strain"] = cfg["scope"]["strain_contains"] if "Strain" not in obs else obs["Strain"]
    obs["Original_Unit"] = c["original_unit"]
    obs["Source_File"] = path.name
    return _finalise("wu_2023", wide, obs, genes[~genes.index.duplicated()], cfg, log, files)


# --------------------------------------------------------------------------- #
# Zhu et al. 2025 - every sheet of every configured file is searched
# --------------------------------------------------------------------------- #
def load_zhu(cfg: Config) -> Dataset:
    """Load the E. coli samples of Zhu et al. 2025 (Supplementary Data 1/2).

    Every sheet of every configured file is searched for the sample columns listed in
    metadata/zhu_2025_samples.tsv; the sheet each sample was found in is logged.
    """
    c = cfg["datasets"]["zhu_2025"]
    meta = _read_metadata(cfg, c["metadata"])
    log, files, blocks, found_in = {}, {}, [], {}
    genes = pd.Series(dtype=object)
    for name in c["files"]:
        path = _raw(cfg, name)
        files[path.name] = sha256(path)
        for sheet, df in pd.read_excel(path, sheet_name=None).items():
            try:
                locus = _col(df, c["locus_col"])
            except KeyError:
                continue                            # sheet without a locus column holds no abundances
            cols = [col for col in df.columns
                    if str(col).strip() in meta.index and str(col).strip() not in found_in]
            if not cols:
                continue
            df = df.dropna(subset=[locus]).copy()
            df[locus] = df[locus].astype(str).str.strip()
            try:
                genes = pd.concat([genes, df.drop_duplicates(locus).set_index(locus)[_col(df, c["gene_col"])]])
            except KeyError:
                pass
            for col in cols:
                found_in[str(col).strip()] = f"{path.name} / {sheet}"
            blocks.append(df.set_index(locus)[cols].rename(columns=lambda s: str(s).strip())
                          .apply(pd.to_numeric, errors="coerce").groupby(level=0).sum())
    if not blocks:
        raise KeyError("zhu_2025: no sheet contains both the locus column and a known sample column.")
    wide = pd.concat(blocks, axis=1)
    log["sample_found_in"] = found_in
    _missing_samples(cfg, "zhu_2025", meta.index, wide.columns, log)

    obs = meta.loc[wide.columns].copy()
    obs["Strain"] = cfg["scope"]["strain_contains"]
    obs["Original_Unit"] = c["original_unit"]
    obs["Source_File"] = [found_in[s] for s in wide.columns]
    return _finalise("zhu_2025", wide, obs, genes[~genes.index.duplicated()], cfg, log, files)


# --------------------------------------------------------------------------- #
# merge + disk I/O
# --------------------------------------------------------------------------- #
def combine(datasets: list[Dataset]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Outer join on b-number. NaN marks 'protein not in this dataset's table'."""
    X = pd.concat([d.X for d in datasets], axis=0, join="outer").sort_index(axis=1)
    obs = pd.concat([d.obs for d in datasets], axis=0)
    if X.index.has_duplicates or not X.index.equals(obs.index):
        raise AssertionError("Sample IDs are not unique / matrix and metadata rows do not match.")
    var = pd.DataFrame(index=X.columns)
    var.index.name = "b_number"
    names = pd.concat([d.genes.rename(d.name) for d in datasets], axis=1).reindex(X.columns)
    var["Gene_Name"] = names.bfill(axis=1).iloc[:, 0]
    for d in datasets:
        present = pd.Series(False, index=X.columns)
        present.loc[d.X.columns] = True
        var[f"in_table_{d.name}"] = present
        var[f"detected_fraction_{d.name}"] = (d.X > 0).mean(axis=0).reindex(X.columns)
    var["n_datasets_in_table"] = var.filter(like="in_table_").sum(axis=1)
    return X, obs, var


def ingest_report(datasets: list[Dataset]) -> pd.DataFrame:
    """One-row-per-dataset summary of the ingest log (sizes, sum checks, drops)."""
    keys = ["n_samples", "n_loci", "detected_per_sample_median", "sample_sum_min", "sample_sum_max",
            "rows_invalid_locus_dropped", "duplicated_locus_rows_summed",
            "cells_non_numeric_set_to_zero", "samples_in_metadata_not_in_files"]
    return pd.DataFrame({d.name: {k: d.log.get(k) for k in keys} for d in datasets}).T


def save_processed(cfg: Config, X: pd.DataFrame, obs: pd.DataFrame, var: pd.DataFrame) -> list[Path]:
    """Write the harmonised matrix, sample metadata and protein table to data/processed/."""
    out = cfg.path("processed_dir")
    paths = [out / "mass_fractions.tsv.gz", out / "sample_metadata.tsv", out / "protein_table.tsv"]
    X.to_csv(paths[0], sep="\t", float_format="%.10g")
    obs.to_csv(paths[1], sep="\t")
    var.to_csv(paths[2], sep="\t")
    return paths


def load_processed(cfg: Config) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read the files written by `save_processed` and check that rows and columns still match."""
    out = cfg.path("processed_dir")
    X = pd.read_csv(out / "mass_fractions.tsv.gz", sep="\t", index_col=0)
    obs = pd.read_csv(out / "sample_metadata.tsv", sep="\t", index_col=0)
    var = pd.read_csv(out / "protein_table.tsv", sep="\t", index_col=0)
    if not X.index.equals(obs.index) or not X.columns.equals(var.index):
        raise AssertionError("Processed files are out of sync; re-run notebook 01.")
    return X, obs, var


def to_anndata(X: pd.DataFrame, obs: pd.DataFrame, var: pd.DataFrame):
    """Optional AnnData export (X = mass fractions, NaN = not in dataset table)."""
    import anndata as ad
    return ad.AnnData(X=X.to_numpy(dtype=float), obs=obs.copy(), var=var.copy())
