#!/usr/bin/env python3
"""
cogaps_run_one_singleprocess.py

Run ONE CoGAPS job (single-process) for a given (K, seed, n_iter) using a cached
genes×cells input matrix.

Key fixes vs earlier version
----------------------------
1) The cache-skip logic now only skips when metrics.status == "ok".
   If a previous run failed, re-running the same job will try again automatically.

2) IFN/condition correlation no longer *requires* 'condition' in cogaps_input.var.
   If missing, it falls back to loading the preprocessed cells file from:
       outdir/cache/preprocessed_cells_*.h5ad
   and aligns by cell IDs.

Usage (typical)
---------------
python cogaps_run_one_singleprocess.py \
  --cogaps-input-h5ad results_cogaps_singleprocess_hpc/cache/cogaps_input_genesxcells_hvg3000_float64.h5ad \
  --outdir results_cogaps_singleprocess_hpc \
  --k 7 --seed 1 --n-iter 2000 \
  --use-sparse-opt \
  --cogaps-threads 1

Notes
-----
- Keep this "single-process" by NOT setting the 'distributed' parameter at all.
- Control BLAS threads via --blas-threads (important on HPC).
"""

from __future__ import annotations

import sys
import json
import time
import argparse
import traceback
import warnings
import os
from pathlib import Path
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

warnings.filterwarnings("ignore", category=FutureWarning, module=r"anndata(\.|$)")
warnings.filterwarnings("ignore", message=r"Importing .* from `anndata` is deprecated.*", category=FutureWarning)

import anndata as ad

from PyCoGAPS.parameters import CoParams, setParams
from PyCoGAPS.pycogaps_main import CoGAPS
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from cogaps_python_utils import (  # noqa: E402
    atomic_write_h5ad,
    atomic_write_json,
    compute_ifn_pattern,
    setup_logger,
    set_blas_threads,
    summarize_result_uns,
    top_genes_table,
    trace_dataframe,
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cogaps-input-h5ad", required=True, help="Cached genes×cells AnnData (dense float64)")
    ap.add_argument("--outdir", default="results_cogaps_singleprocess_hpc", help="Output dir containing runs/ logs/ cache/")
    ap.add_argument("--k", type=int, required=True, help="Number of patterns (K)")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--n-iter", type=int, required=True)
    ap.add_argument("--top-genes", type=int, default=50)
    ap.add_argument("--stim-label", type=str, default="stim", help="Which condition value counts as 'stim' (default: stim)")
    ap.add_argument("--blas-threads", type=int, default=1, help="BLAS/OpenMP threads (important on HPC)")
    ap.add_argument("--cogaps-threads", type=int, default=1, help="Threads passed to CoGAPS(nThreads=...)")
    ap.add_argument("--use-sparse-opt", action="store_true", help="Use sparseOptimization / useSparseOptimization")
    ap.add_argument("--no-sparse-opt", action="store_true", help="Disable sparse optimization")
    ap.add_argument("--output-frequency", type=int, default=1000, help="How often CoGAPS prints trace updates")
    ap.add_argument("--checkpoint-interval", type=int, default=0, help="Checkpoint interval (0 disables checkpoints)")
    ap.add_argument("--checkpoint-out-file", type=str, default="", help="Optional checkpoint output file")
    ap.add_argument("--checkpoint-in-file", type=str, default="", help="Optional checkpoint input file")
    ap.add_argument("--n-snapshots", type=int, default=0, help="Snapshots per phase (0 disables snapshots)")
    ap.add_argument("--snapshot-phase", type=str, default="sampling", choices=["sampling", "equilibration", "all"], help="Which phase(s) to snapshot")
    ap.add_argument("--take-pump-samples", action="store_true", help="Enable pump statistics collection")
    ap.add_argument("--asynchronous-updates", action="store_true", help="Request asynchronous updates when CoGAPS threads > 1")
    ap.add_argument("--preprocessed-h5ad", type=str, default="", help="Optional processed cells x genes H5AD for condition alignment fallback")
    ap.add_argument("--force-rerun", action="store_true", help="Re-run even if metrics exist and status==ok")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    runs_dir = outdir / "runs"
    logs_dir = outdir / "logs"
    runs_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    tag = f"K{args.k}_seed{args.seed}_iter{args.n_iter}"
    run_log = logs_dir / f"run_{tag}.log"
    logger = setup_logger(run_log)

    set_blas_threads(args.blas_threads)

    result_path = runs_dir / f"cogaps_{tag}.h5ad"
    metrics_path = runs_dir / f"cogaps_{tag}.metrics.json"
    diagnostics_path = runs_dir / f"cogaps_{tag}.diagnostics.json"
    trace_path = runs_dir / f"cogaps_{tag}.trace.csv"
    top_genes_path = runs_dir / f"cogaps_{tag}.pattern_top_genes.csv"

    # Cache skip ONLY if status==ok
    if result_path.exists() and metrics_path.exists() and (not args.force_rerun):
        try:
            m = json.loads(metrics_path.read_text(encoding="utf-8"))
            if m.get("status") == "ok":
                logger.info(f"[CACHE] status=ok; skipping {tag}")
                logger.info(f"[CACHE] result:  {result_path}")
                logger.info(f"[CACHE] metrics: {metrics_path}")
                return
            else:
                logger.info(f"[CACHE] Found previous status={m.get('status')}; re-running {tag}")
        except Exception:
            logger.info("[CACHE] Could not parse metrics; re-running.")

    # If previous partial result exists, keep it for debugging but avoid mixing with new writes
    # (Atomic write will overwrite cleanly at the end.)
    t0 = time.time()
    status = "ok"
    err_txt = None
    ifn_pattern = None
    ifn_corr = None
    top_genes = None
    preprocessed_h5ad = Path(args.preprocessed_h5ad) if args.preprocessed_h5ad else None

    try:
        cogaps_input = ad.read_h5ad(args.cogaps_input_h5ad)

        params = CoParams(adata=cogaps_input)

        # IMPORTANT: do NOT set 'distributed' at all (keep single-process).
        sparse_opt = True
        if args.no_sparse_opt:
            sparse_opt = False
        elif args.use_sparse_opt:
            sparse_opt = True

        setParams(
            params,
            {
                "nPatterns": int(args.k),
                "nIterations": int(args.n_iter),
                "seed": int(args.seed),
                "useSparseOptimization": bool(sparse_opt),
                "takePumpSamples": bool(args.take_pump_samples),
            },
        )

        logger.info(
            "[RUN] Starting %s | sparse_opt=%s | blas_threads=%s | cogaps_threads=%s | "
            "output_frequency=%s | checkpoint_interval=%s | n_snapshots=%s | snapshot_phase=%s | pump=%s",
            tag,
            sparse_opt,
            args.blas_threads,
            args.cogaps_threads,
            args.output_frequency,
            args.checkpoint_interval,
            args.n_snapshots,
            args.snapshot_phase,
            args.take_pump_samples,
        )
        result = CoGAPS(
            cogaps_input,
            params,
            nThreads=int(args.cogaps_threads),
            outputFrequency=int(args.output_frequency),
            checkpointOutFile=args.checkpoint_out_file,
            checkpointInterval=int(args.checkpoint_interval),
            checkpointInFile=args.checkpoint_in_file,
            asynchronousUpdates=(True if args.asynchronous_updates else None),
            nSnapshots=int(args.n_snapshots),
            snapshotPhase=args.snapshot_phase,
        )

        logger.info(f"[INFO] Writing result: {result_path}")
        atomic_write_h5ad(result, result_path)

        ifn_pattern, ifn_corr, top_genes = compute_ifn_pattern(
            result=result,
            cogaps_input=cogaps_input,
            preprocessed_cells_h5ad=preprocessed_h5ad,
            outdir=outdir,
            top_n=int(args.top_genes),
            stim_label=args.stim_label,
        )
        logger.info(f"[METRICS] IFN pattern: {ifn_pattern} corr={ifn_corr:.3f}")
        logger.info(f"[METRICS] Top {args.top_genes} genes (first 10): {top_genes[:10]}")

    except Exception as e:
        status = "failed"
        err_txt = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        logger.error("[ERROR] Run failed:\n" + err_txt)

    runtime = time.time() - t0
    diagnostics: Dict[str, Any] = {}
    payload = {
        "K": int(args.k),
        "seed": int(args.seed),
        "n_iter": int(args.n_iter),
        "status": status,
        "language": "Python",
        "package": "PyCoGAPS",
        "sparseOptimization": bool(False if args.no_sparse_opt else True),
        "outputFrequency": int(args.output_frequency),
        "checkpointInterval": int(args.checkpoint_interval),
        "nSnapshots": int(args.n_snapshots),
        "snapshotPhase": args.snapshot_phase,
        "takePumpSamples": bool(args.take_pump_samples),
        "asynchronousUpdates": bool(args.asynchronous_updates),
        "runtime_sec": float(runtime),
        "result_path": str(result_path),
        "metrics_path": str(metrics_path),
        "ifn_pattern": ifn_pattern,
        "ifn_corr": ifn_corr,
        "top_genes": top_genes,
        "error": err_txt,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    if status == "ok" and result_path.exists():
        result_saved = ad.read_h5ad(result_path)
        diagnostics = summarize_result_uns(result_saved)
        atomic_write_json(diagnostics, diagnostics_path)
        trace_dataframe(result_saved).to_csv(trace_path, index=False)
        top_genes_table(result_saved, n_top=int(args.top_genes)).to_csv(top_genes_path, index=False)
        payload.update(diagnostics)

    atomic_write_json(payload, metrics_path)

    logger.info(f"[DONE] {tag} | {status.upper()} | {runtime/60:.2f} min | metrics: {metrics_path}")
    if status != "ok":
        sys.exit(1)


if __name__ == "__main__":
    main()
