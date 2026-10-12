"""Combined view of all datasets: removing the dataset offset for DISPLAY.

Why a separate step
-------------------
The statistics (notebooks 03-05) never alter the data: each batch keeps its own
intercept inside the model. To SHOW the combined proteome - all conditions of
all studies on one axis, PCA, clustering - the batch offsets have to be taken
out of the values. This module does that with the same model:

    y = a_batch + b * growth_rate + residual          (fitted jointly, per variable)
    y_combined = y - (a_batch - mean(a))

Only the batch intercepts are removed; the growth term and the residuals stay.
Fitting batch and growth JOINTLY matters: a study that only measured fast
conditions keeps its fast position. Simply subtracting each batch's mean (what
the sequential ASCA batch matrix does) would pull such a study to the centre
and erase real growth differences between studies.

The adjusted values are for visualisation and exploration only; no test or
model in notebooks 03-05 is fitted on them, so nothing is used twice.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .ingest import GROWTH


def remove_batch_effect(values: pd.DataFrame, info: pd.DataFrame, batch_key: str) -> pd.DataFrame:
    """Rows = conditions (or samples) aligned with `info`; columns = variables."""
    batch = info.loc[values.index, batch_key].astype(str)
    mu = info.loc[values.index, GROWTH].to_numpy(float)
    onehot = pd.get_dummies(batch, dtype=float).to_numpy()
    design = np.column_stack([onehot, mu - mu.mean()])
    Y = values.to_numpy(float)
    beta, *_ = np.linalg.lstsq(design, Y, rcond=None)
    offsets = onehot @ beta[: onehot.shape[1]]                 # a_batch for every row
    adjusted = Y - (offsets - offsets.mean(axis=0, keepdims=True))
    return pd.DataFrame(adjusted, index=values.index, columns=values.columns)


def pca_scores(values: pd.DataFrame, n_components: int = 2, orient: pd.Series | None = None):
    """PCA by SVD of the column-centred matrix (no scaling: log units are comparable).
    PC signs are arbitrary; with `orient` (e.g. growth rate) each PC is flipped so that
    its scores correlate positively with it. Returns (scores DataFrame, explained-variance ratios)."""
    Xc = values.to_numpy(float) - values.to_numpy(float).mean(axis=0)
    U, s, _ = np.linalg.svd(Xc, full_matrices=False)
    k = min(n_components, len(s))
    scores = U[:, :k] * s[:k]
    if orient is not None:
        reference = orient.loc[values.index].to_numpy(float)
        for j in range(k):
            if np.corrcoef(scores[:, j], reference)[0, 1] < 0:
                scores[:, j] *= -1
    ratio = (s ** 2) / np.sum(s ** 2)
    return pd.DataFrame(scores, index=values.index, columns=[f"PC{i + 1}" for i in range(k)]), ratio[:k]


def condition_labels(info: pd.DataFrame) -> pd.Series:
    """Readable condition names: '<condition> (<growth rate>)'."""
    batches = [label.split("|", 1)[0] for label in info.index]
    names = [label.split("|", 1)[-1] for label in info.index]
    text = pd.Series([f"{n} ({g:.2f})" for n, g in zip(names, info[GROWTH])], index=info.index)
    repeated = text.duplicated(keep=False)              # same medium measured in several batches
    text[repeated] = [f"{t} [{b}]" for t, b in zip(text[repeated], pd.Series(batches, index=info.index)[repeated])]
    return text


def proteome_shift_table(adjusted_fractions: pd.DataFrame, info: pd.DataFrame, labels: pd.Series,
                         reference_quantile: float = 0.25) -> pd.DataFrame:
    """For every condition: change of each COG category (percentage points) relative to
    the mean of the slowest conditions (bottom quartile of growth rate), on batch-adjusted
    values. A compact 'how does the proteome change with the condition' table."""
    slow = info[GROWTH] <= info[GROWTH].quantile(reference_quantile)
    reference = adjusted_fractions.loc[slow].mean()
    shift = (adjusted_fractions - reference) * 100
    shift.index = [labels.loc[i] for i in shift.index]
    return shift.round(2)
