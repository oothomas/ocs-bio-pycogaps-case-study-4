#!/usr/bin/env python3
"""
Shared helpers for the local Python CoGAPS experiment scripts.

These utilities keep the single-process and distributed experiment code aligned
on:
- path handling
- runtime logging
- condition alignment
- top-gene extraction
- lightweight diagnostic summaries
"""

from __future__ import annotations

import json
import logging
import os
import sys
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

warnings.filterwarnings("ignore", category=FutureWarning, module=r"anndata(\.|$)")

import anndata as ad
import numpy as np
import pandas as pd


def now_stamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def setup_logger(log_path: Path, name: Optional[str] = None) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(name or f"cogaps_{now_stamp()}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def set_blas_threads(n: int) -> None:
    n = str(int(n))
    os.environ["OMP_NUM_THREADS"] = n
    os.environ["OPENBLAS_NUM_THREADS"] = n
    os.environ["MKL_NUM_THREADS"] = n
    os.environ["NUMEXPR_NUM_THREADS"] = n
    os.environ["VECLIB_MAXIMUM_THREADS"] = n


def atomic_write_h5ad(adata: ad.AnnData, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()
    adata.write_h5ad(tmp_path)
    if tmp_path.stat().st_size == 0:
        raise RuntimeError(f"Atomic write failed: temp file is 0 bytes: {tmp_path}")
    tmp_path.replace(out_path)


def atomic_write_json(payload: Dict[str, Any], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    tmp_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp_path.replace(out_path)


def find_preprocessed_cells(outdir: Path) -> Optional[Path]:
    cache = outdir / "cache"
    if not cache.exists():
        return None
    candidates = sorted(cache.glob("preprocessed_cells_*.h5ad"))
    if not candidates:
        return None
    for path in candidates:
        if "hvg3000" in path.name:
            return path
    return candidates[0]


def ensure_condition_column(adata_cells: ad.AnnData) -> None:
    if "condition" not in adata_cells.obs:
        if "label" in adata_cells.obs:
            adata_cells.obs["condition"] = adata_cells.obs["label"]
        elif "stim" in adata_cells.obs:
            adata_cells.obs["condition"] = adata_cells.obs["stim"]
    if "condition" not in adata_cells.obs:
        raise ValueError("Could not find/create 'condition' in cell metadata.")


def pattern_columns(columns: Iterable[str]) -> List[str]:
    pats = [str(c) for c in columns if str(c).lower().startswith("pattern")]
    return sorted(
        pats,
        key=lambda s: int(str(s).replace("Pattern", "")) if str(s).replace("Pattern", "").isdigit() else 10**9,
    )


def canonicalize_pattern_name(name: str) -> str:
    text = str(name)
    if text.startswith("Pattern_"):
        return "Pattern" + text.split("_", 1)[1]
    return text


def canonicalize_pattern_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [canonicalize_pattern_name(col) for col in out.columns]
    return out


def load_condition_series(
    *,
    result_cell_ids: Sequence[str],
    cogaps_input: Optional[ad.AnnData] = None,
    preprocessed_cells_h5ad: Optional[Path] = None,
    outdir: Optional[Path] = None,
) -> pd.Series:
    cells = pd.Index([str(x) for x in result_cell_ids])

    if cogaps_input is not None and "condition" in cogaps_input.var.columns:
        cond = cogaps_input.var.reindex(cells)["condition"]
        if cond.notna().all():
            return cond

    pre_path = preprocessed_cells_h5ad
    if pre_path is None and outdir is not None:
        pre_path = find_preprocessed_cells(outdir)
    if pre_path is None:
        raise ValueError("No usable condition source found for IFN pattern detection.")

    ad_cells = ad.read_h5ad(str(pre_path))
    ensure_condition_column(ad_cells)
    cond = ad_cells.obs.reindex(cells)["condition"]
    if cond.isna().any():
        raise ValueError("Could not align condition labels to CoGAPS result cells.")
    return cond


def compute_ifn_pattern(
    *,
    result: ad.AnnData,
    cogaps_input: Optional[ad.AnnData],
    preprocessed_cells_h5ad: Optional[Path],
    outdir: Optional[Path],
    top_n: int,
    stim_label: str = "stim",
) -> Tuple[str, float, List[str]]:
    pats = pattern_columns(result.var.columns)
    if not pats:
        raise ValueError("No Pattern* columns found in result.var")

    cond_series = load_condition_series(
        result_cell_ids=result.var_names.astype(str),
        cogaps_input=cogaps_input,
        preprocessed_cells_h5ad=preprocessed_cells_h5ad,
        outdir=outdir,
    )
    stim_indicator = (cond_series.astype(str) == stim_label).astype(int)
    corrs = {pat: float(result.var[pat].corr(stim_indicator)) for pat in pats}
    ifn_pattern = max(corrs, key=lambda key: corrs[key])
    ifn_corr = corrs[ifn_pattern]

    obs_pats = pattern_columns(result.obs.columns)
    if ifn_pattern not in obs_pats:
        raise ValueError(f"Expected {ifn_pattern} in result.obs pattern columns.")
    top_genes = (
        result.obs[ifn_pattern]
        .sort_values(ascending=False)
        .head(int(top_n))
        .index.astype(str)
        .tolist()
    )
    return ifn_pattern, ifn_corr, top_genes


def top_genes_table(result: ad.AnnData, n_top: int) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for pattern in pattern_columns(result.obs.columns):
        ranked = result.obs[pattern].sort_values(ascending=False).head(int(n_top))
        for rank, (gene, weight) in enumerate(ranked.items(), start=1):
            rows.append(
                {
                    "pattern": pattern,
                    "rank": rank,
                    "gene": str(gene),
                    "weight": float(weight),
                }
            )
    return pd.DataFrame(rows)


def safe_shape(value: Any) -> Optional[List[int]]:
    if value is None:
        return None
    arr = np.asarray(value)
    return list(arr.shape)


def summarize_result_uns(result: ad.AnnData) -> Dict[str, Any]:
    uns = result.uns
    summary: Dict[str, Any] = {}
    for key in ("meanChiSq", "totalRunningTime", "totalUpdates", "seed", "averageQueueLengthA", "averageQueueLengthP"):
        if key in uns:
            value = uns[key]
            if isinstance(value, np.generic):
                value = value.item()
            summary[key] = value

    for key in ("chisqHistory", "atomhistoryA", "atomhistoryP"):
        if key in uns:
            try:
                summary[f"{key}_length"] = int(len(uns[key]))
            except TypeError:
                summary[f"{key}_length"] = None

    for key in (
        "equilibrationSnapshotsA",
        "equilibrationSnapshotsP",
        "samplingSnapshotsA",
        "samplingSnapshotsP",
        "meanPatternAssignment",
        "pumpMatrix",
    ):
        if key in uns:
            summary[f"{key}_shape"] = safe_shape(uns[key])

    return summary


def trace_dataframe(result: ad.AnnData) -> pd.DataFrame:
    uns = result.uns
    chisq = list(uns.get("chisqHistory", []))
    atoms_a = list(uns.get("atomhistoryA", []))
    atoms_p = list(uns.get("atomhistoryP", []))
    n = max(len(chisq), len(atoms_a), len(atoms_p))
    if n == 0:
        return pd.DataFrame(columns=["trace_index", "chisq", "atomsA", "atomsP"])
    rows = []
    for i in range(n):
        rows.append(
            {
                "trace_index": i + 1,
                "chisq": chisq[i] if i < len(chisq) else np.nan,
                "atomsA": atoms_a[i] if i < len(atoms_a) else np.nan,
                "atomsP": atoms_p[i] if i < len(atoms_p) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def load_explicit_sets_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    sets = []
    for entry in payload["sets"]:
        sets.append(
            {
                "set_id": entry["set_id"],
                "zero_based_indices": np.asarray(entry["zero_based_indices"], dtype=int),
                "cell_ids": [str(x) for x in entry.get("cell_ids", [])],
            }
        )
    payload["sets"] = sets
    return payload


def explicit_sets_summary_df(payload: Dict[str, Any]) -> pd.DataFrame:
    rows = []
    for entry in payload["sets"]:
        rows.append(
            {
                "set_id": entry["set_id"],
                "n_cells": int(len(entry["zero_based_indices"])),
                "first_index": int(entry["zero_based_indices"][0]) if len(entry["zero_based_indices"]) else np.nan,
                "last_index": int(entry["zero_based_indices"][-1]) if len(entry["zero_based_indices"]) else np.nan,
            }
        )
    return pd.DataFrame(rows)
