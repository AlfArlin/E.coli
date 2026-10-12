"""Single entry point: build the COG annotation, then execute notebooks 01-05 in order.

    python run_all.py            # full run
    python run_all.py --from 3   # restart from notebook 03

Each notebook is executed top to bottom in a fresh kernel and saved in place with
its outputs, so the committed notebooks always show a clean, complete run.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--from", dest="start", type=int, default=1, help="first notebook number to run")
    parser.add_argument("--refresh-annotation", action="store_true", help="re-download UniProt / COG files")
    arguments = parser.parse_args()

    from proteome_alloc import annotation, config
    cfg = config.load_config(ROOT / "config.yml")
    table = annotation.build_annotation(cfg, refresh=arguments.refresh_annotation)
    print(f"[annotation] {len(table)} b-numbers, {(table['COG_Categories'] != '').mean():.1%} with a COG category")

    for notebook in sorted((ROOT / "notebooks").glob("0*.ipynb")):
        if int(notebook.name[:2]) < arguments.start:
            continue
        print(f"[run] {notebook.name}")
        subprocess.run([sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute",
                        "--inplace", "--ExecutePreprocessor.timeout=3600", str(notebook)],
                       check=True, cwd=notebook.parent)
    print("Done. Figures and tables are in results/; the run record is results/run_log.json")


if __name__ == "__main__":
    main()
