"""COG functional annotation - ONE table, keyed by b-number, used for all datasets.

Source chain (both files are downloaded once, cached in data/annotation/ with
their SHA-256, download date and UniProt release, then reused offline):

1. UniProtKB/Swiss-Prot, E. coli K-12 (taxon 83333):
   accession, ordered locus name (b-number), gene, mass, eggNOG cross-reference.
   The eggNOG cross-reference carries the COG identifier(s) of the protein.
2. NCBI COG 2020 ``cog-20.def.tab``: COG identifier -> functional category letter(s).

Using one external annotation for every dataset (instead of the category
columns each paper ships) guarantees that a b-number belongs to the same
functional group in all four studies. The same cached UniProt table also
provides the accession -> b-number map for Schmidt 2016, replacing a live,
non-reproducible API call.

Run once:  python -m proteome_alloc.annotation        (add --refresh to re-download)
"""
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from .config import Config, load_config, sha256

B_NUMBER = re.compile(r"^b\d{4}$")
COG_ID = re.compile(r"^COG\d{4}$")


def _download(url: str, target: Path, params: dict | None = None) -> dict:
    """Download one annotation file and return its provenance record
    (final URL, UTC date, UniProt release header, SHA-256 and size).
    """
    import requests
    response = requests.get(url, params=params, timeout=300)
    response.raise_for_status()
    target.write_bytes(response.content)
    return {"url": response.url, "downloaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "uniprot_release": response.headers.get("X-UniProt-Release"),
            "sha256": sha256(target), "bytes": target.stat().st_size}


def fetch_sources(cfg: Config, refresh: bool = False) -> dict:
    """Download the two source files unless they are already cached."""
    a, folder = cfg["annotation"], cfg.path("annotation_dir")
    provenance_file = folder / "provenance.json"
    provenance = json.loads(provenance_file.read_text()) if provenance_file.exists() else {}
    uniprot, cog_def = folder / a["uniprot_file"], folder / a["cog_def_file"]
    if refresh or not uniprot.exists():
        provenance[uniprot.name] = _download(
            a["uniprot_url"], uniprot,
            {"query": a["uniprot_query"], "fields": a["uniprot_fields"], "format": "tsv"})
    if refresh or not cog_def.exists():
        provenance[cog_def.name] = _download(a["cog_def_url"], cog_def)
    for file in (uniprot, cog_def):                      # manual downloads are fingerprinted too
        provenance.setdefault(file.name, {"url": "placed manually", "sha256": sha256(file)})
        if provenance[file.name]["sha256"] != sha256(file):
            raise ValueError(f"{file.name} changed since it was recorded in provenance.json; "
                             "delete provenance.json or run with --refresh.")
    provenance_file.write_text(json.dumps(provenance, indent=2))
    return provenance


def read_cog_definitions(path: Path) -> pd.Series:
    """COG id -> category letters. Columns of cog-20.def.tab: COG ID, functional
    category (one or more letters, most relevant first), COG name, gene, pathway, PubMed, PDB."""
    table = pd.read_csv(path, sep="\t", header=None, dtype=str, encoding="latin-1",
                        usecols=[0, 1], names=["COG", "letters"], quoting=3)
    table = table[table["COG"].str.match(COG_ID, na=False)]
    return table.drop_duplicates("COG").set_index("COG")["letters"].fillna("")


def build_annotation(cfg: Config, refresh: bool = False) -> pd.DataFrame:
    """Build data/annotation/protein_cog_annotation.tsv (one row per b-number)."""
    a, folder = cfg["annotation"], cfg.path("annotation_dir")
    fetch_sources(cfg, refresh)
    cog_letters = read_cog_definitions(folder / a["cog_def_file"])
    known = set(pd.read_csv(cfg.resolve(a["cog_categories"]), sep="\t", dtype=str,
                            keep_default_na=False)["COG_Category"])
    uni = pd.read_csv(folder / a["uniprot_file"], sep="\t", dtype=str).fillna("")
    columns = {"Entry": "UniProt_Accession", "Gene Names (ordered locus)": "loci",
               "Gene Names (primary)": "Gene_Name", "Mass": "Mass_Da", "eggNOG": "eggnog"}
    missing = set(columns) - set(uni.columns)
    if missing:
        raise KeyError(f"UniProt table lacks columns {missing}; check annotation.uniprot_fields.")
    uni = uni.rename(columns=columns)

    rows = []
    for record in uni.itertuples(index=False):
        loci = [t for t in re.split(r"[\s;,]+", record.loci) if B_NUMBER.match(t)]
        cogs = [t for t in re.split(r"[\s;,]+", record.eggnog) if COG_ID.match(t)]
        letters: list[str] = []
        for cog in cogs:                                  # union of letters, order preserved
            for letter in cog_letters.get(cog, ""):
                if letter in known and letter not in letters:
                    letters.append(letter)
        for locus in loci:
            rows.append({"b_number": locus, "UniProt_Accession": record.UniProt_Accession,
                         "Gene_Name": record.Gene_Name,
                         "Mass_Da": pd.to_numeric(str(record.Mass_Da).replace(",", ""), errors="coerce"),
                         "COG_IDs": ";".join(cogs),
                         "COG_IDs_not_in_COG2020": ";".join(c for c in cogs if c not in cog_letters.index),
                         "COG_Categories": "".join(letters)})
    table = pd.DataFrame(rows)
    table["n_entries_for_locus"] = table.groupby("b_number")["b_number"].transform("size")
    table = table.sort_values(["b_number", "UniProt_Accession"]).drop_duplicates("b_number")
    table.to_csv(folder / a["output_file"], sep="\t", index=False)
    return table


def load_annotation(cfg: Config) -> pd.DataFrame:
    """Read the cached b-number -> COG table built by `build_annotation`.

    Empty COG fields are kept as empty strings (not NaN) so that proteins without
    a COG stay in the analysis as their own group instead of being dropped.
    """
    path = cfg.path("annotation_dir") / cfg["annotation"]["output_file"]
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Run:  python -m proteome_alloc.annotation")
    return pd.read_csv(path, sep="\t", dtype={"COG_Categories": str, "COG_IDs": str},
                       keep_default_na=False, na_values=[""]).fillna(
        {"COG_Categories": "", "COG_IDs": "", "COG_IDs_not_in_COG2020": ""})


def category_weights(loci: pd.Index, annotation: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """Proteins x categories weight matrix whose rows sum to exactly 1.

    multi_category = 'split': a protein in k categories gives 1/k of its mass to
    each (no double counting, category fractions still sum to the proteome);
    'first': all mass to the first (most relevant) letter.
    Proteins without any COG go to the unassigned label so that mass is never lost."""
    a = cfg["annotation"]
    categories = pd.read_csv(cfg.resolve(a["cog_categories"]), sep="\t", dtype=str,
                             keep_default_na=False)["COG_Category"].tolist()
    unassigned = a["unassigned_label"]
    letters = annotation.set_index("b_number")["COG_Categories"].reindex(loci).fillna("")
    W = pd.DataFrame(0.0, index=loci, columns=categories)
    for locus, value in letters.items():
        chosen = list(value) or [unassigned]
        if a["multi_category"] == "first":
            chosen = chosen[:1]
        elif a["multi_category"] != "split":
            raise ValueError("annotation.multi_category must be 'split' or 'first'")
        W.loc[locus, chosen] = 1.0 / len(chosen)
    assert np.allclose(W.sum(axis=1), 1.0)
    return W


def category_labels(cfg: Config) -> pd.Series:
    """COG category letter -> human-readable description (metadata/cog_categories.tsv)."""
    table = pd.read_csv(cfg.resolve(cfg["annotation"]["cog_categories"]), sep="\t", dtype=str,
                        keep_default_na=False)
    return table.set_index("COG_Category")["Description"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build the b-number -> COG annotation table.")
    parser.add_argument("--refresh", action="store_true", help="re-download the source files")
    arguments = parser.parse_args()
    config = load_config()
    result = build_annotation(config, refresh=arguments.refresh)
    share = (result["COG_Categories"] != "").mean()
    print(f"{len(result)} b-numbers annotated; {share:.1%} have at least one COG category.")
    print("Written:", config.path("annotation_dir") / config["annotation"]["output_file"])
