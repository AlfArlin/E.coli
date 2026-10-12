"""Configuration, paths, deterministic random streams and the run log.

Reproducibility rules implemented here
--------------------------------------
* one YAML file holds every parameter (no values hidden in notebooks);
* every path is relative to the repository root (no machine-specific paths);
* every stochastic step draws from a stream spawned from ONE master seed, so
  results are identical whatever the number of parallel workers;
* every input file is fingerprinted (SHA-256) and every notebook appends its
  parameters, package versions and outputs to results/run_log.json.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path

import numpy as np
import yaml

TRACKED_PACKAGES = ("numpy", "pandas", "scipy", "scikit-learn", "matplotlib",
                    "openpyxl", "PyYAML", "joblib", "requests", "anndata")


def find_root(start: Path | None = None) -> Path:
    """Walk upwards until the folder containing config.yml is found."""
    here = Path(start or Path.cwd()).resolve()
    for folder in (here, *here.parents):
        if (folder / "config.yml").exists() and (folder / "src").exists():
            return folder
    raise FileNotFoundError("config.yml not found: run from inside the repository.")


class Config(dict):
    """Dictionary with the repository root attached and path helpers."""

    root: Path

    def path(self, key: str) -> Path:
        p = self.root / self["paths"][key]
        return p

    def resolve(self, relative: str | os.PathLike) -> Path:
        p = Path(relative)
        return p if p.is_absolute() else self.root / p

    @property
    def figures(self) -> Path:
        p = self.path("results_dir") / "figures"
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def tables(self) -> Path:
        p = self.path("results_dir") / "tables"
        p.mkdir(parents=True, exist_ok=True)
        return p

    def rng(self, stream: str) -> np.random.Generator:
        """Independent, named random stream derived from the master seed."""
        return np.random.default_rng(self.seed_sequence(stream))

    def seed_sequence(self, stream: str) -> np.random.SeedSequence:
        tag = int.from_bytes(hashlib.sha256(stream.encode()).digest()[:4], "big")
        return np.random.SeedSequence([int(self["seed"]), tag])


def load_config(path: str | os.PathLike | None = None) -> Config:
    """Load config.yml. The environment variable PALLOC_CONFIG can point to an
    alternative file (used by the smoke test so real data are never touched)."""
    env = os.environ.get("PALLOC_CONFIG")
    if path is None and env:
        path = env
    if path is None:
        root = find_root()
        path = root / "config.yml"
    path = Path(path).resolve()
    with open(path, encoding="utf-8") as fh:
        cfg = Config(yaml.safe_load(fh))
    cfg.root = Path(cfg.get("root_dir") or path.parent).resolve()
    cfg.config_file = path
    for key in ("processed_dir", "annotation_dir", "results_dir"):
        cfg.path(key).mkdir(parents=True, exist_ok=True)
    return cfg


def sha256(path: str | os.PathLike, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def environment_record() -> dict:
    versions = {}
    for pkg in TRACKED_PACKAGES:
        try:
            versions[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            versions[pkg] = None
    return {"python": sys.version.split()[0], "platform": platform.platform(),
            "packages": versions}


def log_run(cfg: Config, step: str, params: dict | None = None,
            inputs: dict | None = None, outputs: list | None = None,
            notes: dict | None = None) -> dict:
    """Append one entry per pipeline step to results/run_log.json."""
    log_file = cfg.path("results_dir") / "run_log.json"
    log = json.loads(log_file.read_text()) if log_file.exists() else {}
    entry = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "config_file": str(cfg.config_file.name),
        "config_sha256": sha256(cfg.config_file),
        "seed": cfg["seed"],
        "environment": environment_record(),
        "parameters": params or {},
        "inputs_sha256": inputs or {},
        "outputs": [str(o) for o in (outputs or [])],
        "notes": notes or {},
    }
    log[step] = entry          # re-running a step replaces its entry
    log_file.write_text(json.dumps(log, indent=2, default=str))
    return entry
