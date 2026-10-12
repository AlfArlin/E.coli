"""Create SYNTHETIC stand-ins for the raw supplementary files.

Purpose: a smoke test. The files have the real file names, sheet names and
column names, but RANDOM numbers with a known planted signal (COG category J
rises with growth rate; C, E and G fall). They let anyone check, without the
real data, that the pipeline installs and runs from beginning to end and
recovers the planted signal. They are written to tests/_synthetic/ only and
must never be copied into data/raw/.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tests" / "_synthetic"
LETTERS = ["J", "C", "E", "G", "K", "M", "O", "P", "F", "H"]
SLOPE = {"J": 0.9, "C": -0.35, "E": -0.45, "G": -0.5}        # log-fold change per h-1


def main() -> Path:
    rng = np.random.default_rng(1)
    if OUT.exists():
        shutil.rmtree(OUT)
    raw, ann = OUT / "data" / "raw", OUT / "data" / "annotation"
    raw.mkdir(parents=True), ann.mkdir(parents=True)
    shutil.copytree(ROOT / "metadata", OUT / "metadata")

    n = 400
    loci = [f"b{i:04d}" for i in range(1, n + 1)]
    genes = [f"gen{i}" for i in range(1, n + 1)]
    acc = [f"P{i:05d}" for i in range(1, n + 1)]
    mw = rng.uniform(1.5e4, 1.2e5, n)
    letter = rng.choice(LETTERS, n)
    cog = {L: f"COG{1000 + i:04d}" for i, L in enumerate(LETTERS)}
    has_cog = rng.random(n) > 0.12
    two = rng.random(n) > 0.9
    base = rng.lognormal(0, 1.6, n)

    def profile(mu: float, bias: np.ndarray, noise: float = 0.08, drop: float = 0.05) -> np.ndarray:
        effect = np.array([SLOPE.get(L, 0.0) for L in letter])
        values = base * np.exp(effect * mu) * bias * rng.lognormal(0, noise, n)
        values[rng.random(n) < drop * (base < np.median(base))] = 0.0
        return values / values.sum()

    # ---- annotation sources -------------------------------------------------
    eggnog = [(f"{cog[L]};" + (f"{cog['K']};" if t else "")) if h else "ENOG502ZZZZ;"
              for L, h, t in zip(letter, has_cog, two)]
    pd.DataFrame({"Entry": acc, "Gene Names (ordered locus)": [f"{b} JW{i:04d}" for i, b in enumerate(loci)],
                  "Gene Names (primary)": genes, "Mass": mw.round().astype(int), "eggNOG": eggnog}
                 ).to_csv(ann / "uniprot_ecoli_k12.tsv", sep="\t", index=False)
    with open(ann / "cog-20.def.tab", "w", encoding="latin-1") as fh:
        for L in LETTERS:
            fh.write(f"{cog[L]}\t{L}\tsynthetic COG {L}\tgene\tpathway\t\t\n")

    # ---- Schmidt: copies per cell, UniProt accessions, 60% of the proteome ---
    bias = rng.lognormal(0, 0.25, n)
    seen = rng.random(n) < 0.6
    table = pd.DataFrame({"Uniprot Accession": acc, "Description": "synthetic", "Gene": genes, "proteinMW": mw})
    for name, mu in (("LB", 2.30), ("Glucose", 1.03)):
        table[f"Copies/Cell_NCM3722.{name}"] = profile(mu, bias) / mw * 6e10
        table[f"Copies/Cell_BW25113.{name}"] = profile(mu * 0.8, bias) / mw * 6e10
    table = table[seen]
    table.loc[len(acc) + 1] = ["Q99999", "synthetic", "unmapped", 5e4, 10.0, 10.0, 10.0, 10.0]
    growth = pd.DataFrame({"Growth condition": ["LB", "Glucose", "LB", "Glucose"],
                           "Strain": ["NCM3722", "NCM3722", "BW25113", "BW25113"],
                           "Growth rate (h-1)": [2.30, 1.03, 1.90, 0.58], "Stdev": 0.05})
    with pd.ExcelWriter(raw / "1.1The quantitative and condition-dependent Escherichia coli proteome.xlsx") as xl:
        table.to_excel(xl, sheet_name="Table S9", index=False, startrow=2)
        growth.to_excel(xl, sheet_name="Table S23", index=False, startrow=2)

    # ---- Mori: 19 library samples, metadata in a second file ----------------
    bias = rng.lognormal(0, 0.25, n)
    rows = [("Hi Temp", "MOPS (Record's)", "0.2% glucose", "", 47), ("Hi Osm 1", "MOPS (Record's)", "0.2% glucose", "", 170),
            ("Hi Osm 2", "MOPS (Record's)", "0.2% fructose", "", 110), ("Hi Osm 3", "MOPS (Record's)", "0.2% glycerol", "", 150),
            ("Anaerobic", "Special N-C- (see Methods)", "0.2% glucose", "", 80),
            ("Glucose", "MOPS (Record's)", "0.2% glucose", "", 43), ("Glycerol", "MOPS (Record's)", "0.2% glycerol", "", 62),
            ("Acetate", "MOPS (Record's)", "0.2% acetate", "", 95), ("Mannose", "MOPS (Record's)", "0.2% mannose", "", 110),
            ("Glucose+CAA", "MOPS (Record's)", "0.2% glucose", "", 30), ("Rich", "MOPS (Record's)", "0.2% glucose", "", 24),
            ("Stationary", "MOPS (Record's)", "0.2% glucose", "", 0), ("Low pH", "MOPS (Record's)", "0.2% glucose", "", 70),
            ("Cm 2uM", "MOPS (Record's)", "0.2% glucose", "", 75), ("H2O2", "MOPS (Record's)", "0.2% glucose", "", 60),
            ("Xylose", "MOPS (Record's)", "0.2% xylose", "", 58), ("Lactose", "MOPS (Record's)", "0.2% lactose", "", 46),
            ("Glucose shift", "MOPS (Record's)", "0.2% glucose", "sampled after upshift", 40),
            ("Pyruvate", "MOPS (Record's)", "0.2% pyruvate", "", np.nan)]
    ids = [f"Lib-{i:02d}" for i in range(1, 20)]
    meta = pd.DataFrame(rows, columns=["Short Description", "Base medium", "Carbon Source",
                                       "Special condition/procedure", "Doubling time (min)"])
    meta.insert(0, "Sample ID", ids)
    meta["Description"] = "synthetic"
    meta["Strain"] = "NCM3722"
    meta["Nitrogen Source"], meta["Supplement"] = "NH4Cl", ""
    meta.loc[len(meta)] = ["Lib-99", "Other strain", "MOPS (Record's)", "0.2% glucose", "", 40, "synthetic", "MG1655", "NH4Cl", ""]
    meta.to_excel(raw / "2.6MOESM2_ESM.xlsx", sheet_name="EV2-Samples-1", index=False)
    mori = pd.DataFrame({"Gene locus": loci, "Gene name": genes})
    for sid, row in zip(ids, rows):
        mu = np.log(2) / (row[4] / 60) if row[4] and row[4] == row[4] and row[4] > 0 else 0.3
        stress = rng.lognormal(0, 0.6, n) if row[0] in {"Hi Temp", "Hi Osm 1", "Anaerobic", "Cm 2uM"} else 1.0
        mori[sid] = profile(mu, bias * stress)
    mori.loc[len(mori)] = [np.nan, "no locus"] + [0.0] * 19
    mori.to_excel(raw / "2.3MOESM8_ESM.xlsx", sheet_name="EV8-AbsoluteMassFractions-1", index=False)

    # ---- Wu: three measurement groups ---------------------------------------
    wu = pd.read_csv(ROOT / "metadata" / "wu_2023_samples.tsv", sep="\t")
    with pd.ExcelWriter(raw / "3.1CHEN41564_2022_1310_MOESM3_ESM.xlsx") as xl:
        for group, sheet in (("Group1", "Group #1 protein mass fractions"), ("Group2", "Group #2 protein mass fractions"),
                             ("Group3", "Group #3 protein mass fractions")):
            bias = rng.lognormal(0, 0.25, n) * rng.lognormal(0, 0.05, n)
            df = pd.DataFrame({"Gene name": genes, "Gene locus": loci})
            for r in wu[(wu["MS_Batch"] == group) & (wu["Original_Sample_ID"] != "A3")].itertuples():
                df[r.Original_Sample_ID] = profile(r._10, bias)
            df.loc[len(df)] = ["zero locus", 0] + [0.0] * (df.shape[1] - 2)
            df.to_excel(xl, sheet_name=sheet, index=False)

    # ---- Zhu: 37 C samples in SD01, 30 C samples in a sheet of SD02 ----------
    zhu = pd.read_csv(ROOT / "metadata" / "zhu_2025_samples.tsv", sep="\t")
    bias = rng.lognormal(0, 0.25, n)
    cold = rng.lognormal(0, 0.3, n)
    frames = {}
    for flag, part in zhu.groupby(zhu["Original_Sample_ID"].str.startswith("K")):
        df = pd.DataFrame({"gene name": genes, "gene locus": loci})
        for r in part.itertuples():
            df[r.Original_Sample_ID] = profile(r._10, bias * (cold if r.Temperature_C == 30 else 1.0))
        frames[flag] = df
    with pd.ExcelWriter(raw / "4.1. Zhu_SD01.xlsx") as xl:
        frames[False].to_excel(xl, sheet_name="Sheet1", index=False)
    with pd.ExcelWriter(raw / "4.2Zhu_SD02.xlsx") as xl:
        pd.DataFrame({"note": ["synthetic sheet without abundances"]}).to_excel(xl, sheet_name="readme", index=False)
        frames[True].to_excel(xl, sheet_name="abundances", index=False)

    # ---- test configuration: same parameters, smaller resampling, own root ---
    cfg = yaml.safe_load((ROOT / "config.yml").read_text())
    cfg["root_dir"] = str(OUT)
    cfg["allocation"].update(n_permutations=999, n_bootstrap=400)
    cfg["asca"]["n_permutations"] = 499
    cfg["prediction"]["n_permutations"] = 30
    target = OUT / "config_test.yml"
    target.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return target


if __name__ == "__main__":
    print(main())
    sys.exit(0)
