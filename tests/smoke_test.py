"""End-to-end smoke test on SYNTHETIC data (no real data, no internet, no Jupyter).

    python tests/smoke_test.py

1. builds synthetic raw files with the real file / sheet / column layout;
2. executes the code cells of the five notebooks in order against them;
3. checks that the planted signal is recovered and that a second run gives
   byte-identical result tables (reproducibility check).
Outputs go to tests/_synthetic/ only.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
os.environ.setdefault("MPLBACKEND", "Agg")

import make_synthetic_data  # noqa: E402


def run_notebooks() -> dict:
    import matplotlib.pyplot as plt
    hashes = {}
    for notebook in sorted((ROOT / "notebooks").glob("0*.ipynb")):
        cells = json.loads(notebook.read_text(encoding="utf-8"))["cells"]
        namespace = {"display": lambda *objects: [print(o) for o in objects], "__name__": "__main__"}
        os.chdir(ROOT / "notebooks")
        for number, cell in enumerate(cells):
            if cell["cell_type"] != "code":
                continue
            try:
                exec(compile("".join(cell["source"]), f"{notebook.name}[cell {number}]", "exec"), namespace)
            except Exception:
                print(f"\nFAILED in {notebook.name}, cell {number}", file=sys.stderr)
                raise
        plt.close("all")
        print(f"--- {notebook.name}: ok")
    tables = Path(namespace["cfg"].tables)
    for table in sorted(tables.glob("*.tsv")):
        hashes[table.name] = hashlib.sha256(table.read_bytes()).hexdigest()
    return hashes


def main() -> None:
    import pandas as pd
    os.environ["PALLOC_CONFIG"] = str(make_synthetic_data.main())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        first = run_notebooks()
        second = run_notebooks()

    tables = make_synthetic_data.OUT / "results" / "tables"
    slopes = pd.read_csv(tables / "03_cog_growth_slopes.tsv", sep="\t", index_col=0)
    assert slopes.loc["J", "slope_per_h-1"] > 0 and bool(slopes.loc["J", "significant"]), "planted rise of J not found"
    for category in ("C", "E", "G"):
        assert slopes.loc[category, "slope_per_h-1"] < 0, f"planted fall of {category} not found"
    assert abs(slopes["slope_per_h-1"].sum()) < 1e-9, "category slopes must sum to zero"
    scope = pd.read_csv(tables / "02_scope_decisions.tsv", sep="\t", index_col=0)
    assert not scope.loc[scope["Temperature_C"] == 30, "In_Scope"].any(), "30 C samples must be out of scope"
    changed = [name for name in first if first[name] != second[name]]
    assert not changed, f"tables differ between two identical runs: {changed}"
    print(f"\nSMOKE TEST PASSED - {len(first)} result tables identical across two runs.")


if __name__ == "__main__":
    main()
