#!/usr/bin/env python3
"""
Run a controlled distributed PyCoGAPS experiment with fixed explicit subsets.

This runner is intentionally separate from the older notebook-oriented Python
distributed script. It is designed to mirror the current R distributed
comparison as closely as PyCoGAPS allows:

- same genes x cells input matrix
- same explicit single-cell subsets
- same K / seed / iterations
- same sparse on/off toggle
- same output-frequency / snapshot settings

Important caveat
----------------
PyCoGAPS reimplements the distributed consensus step in Python, so we should
expect biological concordance, not bitwise identity, when comparing to the R
distributed runs.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing
import time
import traceback
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

warnings.filterwarnings("ignore", category=FutureWarning, module=r"anndata(\.|$)")

import anndata as ad
import numpy as np
import pandas as pd
import pycogaps
from scipy import sparse as sp

from PyCoGAPS.distributed_functions import corrToMeanPattern, findConsensusMatrix, stitchTogether
from PyCoGAPS.parameters import CoParams, setParams
from PyCoGAPS.pycogaps_main import standardCoGAPS

SCRIPT_DIR = Path(__file__).resolve().parent
import sys

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from cogaps_python_utils import (  # noqa: E402
    atomic_write_h5ad,
    atomic_write_json,
    compute_ifn_pattern,
    explicit_sets_summary_df,
    load_explicit_sets_json,
    pattern_columns,
    set_blas_threads,
    setup_logger,
    summarize_result_uns,
    top_genes_table,
    trace_dataframe,
)

def build_base_params(
    adata: ad.AnnData,
    *,
    n_patterns: int,
    n_iterations: int,
    seed: int,
    use_sparse_opt: bool,
    distributed_mode: str,
    n_sets: int,
    cut: int,
    min_ns: int,
    max_ns: int,
) -> CoParams:
    params = CoParams(adata=adata)
    setParams(
        params,
        {
            "nPatterns": int(n_patterns),
            "nIterations": int(n_iterations),
            "seed": int(seed),
            "useSparseOptimization": bool(use_sparse_opt),
        },
    )
    params = setParamSafe(params, "distributed", distributed_mode)
    params.setDistributedParams(nSets=int(n_sets), cut=int(cut), minNS=int(min_ns), maxNS=int(max_ns))
    return params


def setParamSafe(params: CoParams, key: str, value: Any) -> CoParams:
    setParams(params, {key: value})
    return params


def run_subset_worker(task: Dict[str, Any]) -> Dict[str, Any]:
    data_path = task["data_path"]
    subset_indices = np.asarray(task["subset_indices"], dtype=int)
    worker_id = int(task["worker_id"])
    run_cfg = task["run_cfg"]
    fixed_consensus = task.get("fixed_consensus")

    adata = ad.read_h5ad(data_path)
    if run_cfg["worker_source_layout"] == "cellsxgenes":
        if run_cfg["distributed_mode"] == "single-cell":
            adata = adata[subset_indices, :].copy()
        else:
            adata = adata[:, subset_indices].copy()
        adata = adata.T.copy()
    else:
        if run_cfg["distributed_mode"] == "single-cell":
            adata = adata[:, subset_indices].copy()
        else:
            adata = adata[subset_indices, :].copy()

    if sp.issparse(adata.X):
        adata.X = np.asarray(adata.X.toarray(), dtype=np.float64)
    else:
        adata.X = np.asarray(adata.X, dtype=np.float64)

    params = CoParams(adata=adata)
    setParams(
        params,
        {
            "nPatterns": int(run_cfg["n_patterns"]),
            "nIterations": int(run_cfg["n_iterations"]),
            "seed": int(run_cfg["seed"]),
            "useSparseOptimization": bool(run_cfg["use_sparse_opt"]),
            "takePumpSamples": bool(run_cfg["take_pump_samples"]),
        },
    )

    if fixed_consensus is not None:
        fixed_arr = np.asarray(fixed_consensus, dtype=float)
        params.gaps.nPatterns = fixed_arr.shape[1]
        params.gaps.useFixedPatterns = True
        params.gaps.fixedPatterns = pycogaps.Matrix(fixed_arr)
        params.gaps.whichMatrixFixed = "A" if run_cfg["distributed_mode"] == "single-cell" else "P"

    result = standardCoGAPS(
        adata,
        params,
        nThreads=1,
        messages=bool(run_cfg["messages"]),
        outputFrequency=int(run_cfg["output_frequency"]),
        checkpointOutFile="",
        checkpointInterval=0,
        checkpointInFile="",
        transposeData=False,
        workerID=worker_id,
        asynchronousUpdates=False,
        nSnapshots=int(run_cfg["n_snapshots"]),
        snapshotPhase=run_cfg["snapshot_phase"],
    )

    return {
        "worker_id": worker_id,
        "subset_size": int(len(subset_indices)),
        "result": result,
    }


def summarize_worker(entry: Dict[str, Any]) -> Dict[str, Any]:
    result = entry["result"]
    uns_summary = summarize_result_uns(result)
    return {
        "worker_id": int(entry["worker_id"]),
        "subset_size": int(entry["subset_size"]),
        "n_patterns": int(len(pattern_columns(result.obs.columns if result.obs.shape[1] else result.var.columns))),
        **uns_summary,
    }


def summarize_consensus(matched: Dict[str, Any]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for idx, cluster in enumerate(matched["clustered"], start=1):
        if cluster is None:
            continue
        corr = np.asarray(corrToMeanPattern(cluster), dtype=float)
        rows.append(
            {
                "cluster": idx,
                "n_members": int(cluster.shape[1]),
                "min_corr_to_consensus": float(np.min(corr)),
                "mean_corr_to_consensus": float(np.mean(corr)),
                "max_corr_to_consensus": float(np.max(corr)),
            }
        )
    return pd.DataFrame(rows)


def build_result_adata(
    full_data: ad.AnnData,
    stitched: Dict[str, pd.DataFrame],
    final_results: List[Dict[str, Any]],
) -> ad.AnnData:
    adata = full_data.copy()
    first_final = final_results[0]["result"]

    adata.obs = stitched["Amean"]
    adata.var = stitched["Pmean"]
    adata.uns["asd"] = stitched["Asd"]
    adata.uns["psd"] = stitched["Psd"]

    for key in (
        "atomhistoryA",
        "atomhistoryP",
        "averageQueueLengthA",
        "averageQueueLengthP",
        "chisqHistory",
        "equilibrationSnapshotsA",
        "equilibrationSnapshotsP",
        "meanChiSq",
        "meanPatternAssignment",
        "pumpMatrix",
        "samplingSnapshotsA",
        "samplingSnapshotsP",
        "seed",
        "totalRunningTime",
        "totalUpdates",
    ):
        if key in first_final.uns:
            adata.uns[key] = first_final.uns[key]
    return adata


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cogaps-input-h5ad", required=True, help="Genes x cells cached CoGAPS AnnData (.h5ad)")
    ap.add_argument("--preprocessed-h5ad", default="", help="Optional processed cells x genes H5AD for condition alignment")
    ap.add_argument("--worker-source-h5ad", default="", help="Optional H5AD loaded by distributed workers")
    ap.add_argument("--worker-source-layout", default="genesxcells", choices=["genesxcells", "cellsxgenes"], help="Orientation of the worker source H5AD")
    ap.add_argument("--explicit-sets-json", required=True, help="Python-friendly explicit set JSON exported from the R subsets")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--n-iter", type=int, required=True)
    ap.add_argument("--top-genes", type=int, default=50)
    ap.add_argument("--stim-label", type=str, default="stim")
    ap.add_argument("--distributed-mode", type=str, default="single-cell", choices=["single-cell", "genome-wide"])
    ap.add_argument("--n-sets", type=int, default=4)
    ap.add_argument("--cut", type=int, default=None)
    ap.add_argument("--min-ns", type=int, default=None)
    ap.add_argument("--max-ns", type=int, default=None)
    ap.add_argument("--pool-workers", type=int, default=None, help="Worker processes for the subset pools [default nSets]")
    ap.add_argument("--output-frequency", type=int, default=100)
    ap.add_argument("--n-snapshots", type=int, default=10)
    ap.add_argument("--snapshot-phase", type=str, default="all", choices=["sampling", "equilibration", "all"])
    ap.add_argument("--take-pump-samples", action="store_true")
    ap.add_argument("--messages", action="store_true", help="Print CoGAPS worker status messages")
    ap.add_argument("--use-sparse-opt", action="store_true", help="Enable sparse optimization")
    ap.add_argument("--no-sparse-opt", action="store_true", help="Disable sparse optimization")
    ap.add_argument("--force-rerun", action="store_true")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    logs_dir = outdir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    sparse_opt = False if args.no_sparse_opt else True
    tag = f"K{args.k}_seed{args.seed}_iter{args.n_iter}_{args.distributed_mode.replace('-', '_')}"
    logger = setup_logger(logs_dir / f"run_{tag}.log", name=f"py_dist_{tag}")

    result_path = outdir / f"cogaps_{tag}.h5ad"
    metrics_path = outdir / f"cogaps_{tag}.metrics.json"
    diagnostics_path = outdir / f"cogaps_{tag}.diagnostics.json"
    subset_summary_path = outdir / f"cogaps_{tag}.subset_summary.csv"
    first_pass_summary_path = outdir / f"cogaps_{tag}.first_pass_worker_summary.csv"
    final_stage_summary_path = outdir / f"cogaps_{tag}.final_stage_worker_summary.csv"
    consensus_summary_path = outdir / f"cogaps_{tag}.consensus_summary.csv"
    top_genes_path = outdir / f"cogaps_{tag}.pattern_top_genes.csv"
    trace_path = outdir / f"cogaps_{tag}.trace.csv"
    consensus_path = outdir / f"cogaps_{tag}.consensus_patterns.csv"

    if result_path.exists() and metrics_path.exists() and not args.force_rerun:
        try:
            prior = json.loads(metrics_path.read_text(encoding="utf-8"))
            if prior.get("status") == "ok":
                logger.info("[CACHE] status=ok; skipping %s", tag)
                return
        except Exception:
            logger.info("[CACHE] Could not parse cached metrics; re-running.")

    set_blas_threads(1)

    status = "ok"
    err_txt = None
    t0 = time.time()

    try:
        full_data = ad.read_h5ad(args.cogaps_input_h5ad)
        explicit_payload = load_explicit_sets_json(Path(args.explicit_sets_json))
        explicit_sets = [entry["zero_based_indices"] for entry in explicit_payload["sets"]]

        if len(explicit_sets) != int(args.n_sets):
            raise ValueError(f"JSON contains {len(explicit_sets)} sets, expected {args.n_sets}.")

        if args.distributed_mode == "single-cell":
            if full_data.n_vars == 0:
                raise ValueError("Expected genes x cells input with cells in .var")
            for entry in explicit_payload["sets"]:
                if entry["cell_ids"]:
                    observed = full_data.var_names[entry["zero_based_indices"]].astype(str).tolist()
                    if observed != entry["cell_ids"]:
                        raise ValueError(
                            f"Cell ordering mismatch for {entry['set_id']} between the exported R sets and Python input."
                        )

        subset_summary = explicit_sets_summary_df(explicit_payload)
        subset_summary.to_csv(subset_summary_path, index=False)

        cut_value = int(args.cut if args.cut is not None else args.k)
        min_ns_value = int(args.min_ns if args.min_ns is not None else np.ceil(args.n_sets / 2))
        max_ns_value = int(args.max_ns if args.max_ns is not None else (min_ns_value + args.n_sets))
        pool_workers = int(args.pool_workers if args.pool_workers is not None else args.n_sets)
        worker_source_h5ad = args.worker_source_h5ad or args.cogaps_input_h5ad

        base_params = build_base_params(
            full_data,
            n_patterns=int(args.k),
            n_iterations=int(args.n_iter),
            seed=int(args.seed),
            use_sparse_opt=bool(sparse_opt),
            distributed_mode=args.distributed_mode,
            n_sets=int(args.n_sets),
            cut=cut_value,
            min_ns=min_ns_value,
            max_ns=max_ns_value,
        )

        run_cfg = {
            "distributed_mode": args.distributed_mode,
            "n_patterns": int(args.k),
            "n_iterations": int(args.n_iter),
            "seed": int(args.seed),
            "use_sparse_opt": bool(sparse_opt),
            "take_pump_samples": bool(args.take_pump_samples),
            "messages": bool(args.messages),
            "output_frequency": int(args.output_frequency),
            "n_snapshots": int(args.n_snapshots),
            "snapshot_phase": args.snapshot_phase,
            "worker_source_layout": args.worker_source_layout,
        }

        worker_tasks = [
            {
                "data_path": worker_source_h5ad,
                "subset_indices": subset.tolist(),
                "worker_id": i + 1,
                "run_cfg": run_cfg,
                "fixed_consensus": None,
            }
            for i, subset in enumerate(explicit_sets)
        ]

        logger.info("[RUN] First pass | sparse_opt=%s | distributed_mode=%s | nSets=%s", sparse_opt, args.distributed_mode, args.n_sets)
        with multiprocessing.get_context("spawn").Pool(processes=pool_workers) as pool:
            first_pass_results = pool.map(run_subset_worker, worker_tasks)

        first_pass_summary = pd.DataFrame([summarize_worker(entry) for entry in first_pass_results])
        first_pass_summary.to_csv(first_pass_summary_path, index=False)

        if args.distributed_mode == "genome-wide":
            unmatched = np.array(first_pass_results[0]["result"].var)
            for entry in first_pass_results[1:]:
                unmatched = np.append(unmatched, np.array(entry["result"].var), axis=1)
        else:
            unmatched = np.array(first_pass_results[0]["result"].obs)
            for entry in first_pass_results[1:]:
                unmatched = np.append(unmatched, np.array(entry["result"].obs), axis=1)

        logger.info("[RUN] Consensus matching")
        matched = findConsensusMatrix(unmatched, base_params)
        matched["consensus"].to_csv(consensus_path, index=False)
        consensus_summary = summarize_consensus(matched)
        consensus_summary.to_csv(consensus_summary_path, index=False)

        fixed_consensus = matched["consensus"].to_numpy(dtype=float)
        final_tasks = []
        for i, subset in enumerate(explicit_sets):
            final_tasks.append(
                {
                    "data_path": worker_source_h5ad,
                    "subset_indices": subset.tolist(),
                    "worker_id": i + 1,
                    "run_cfg": {
                        **run_cfg,
                        "n_patterns": int(fixed_consensus.shape[1]),
                    },
                    "fixed_consensus": fixed_consensus,
                }
            )

        base_params.gaps.nPatterns = fixed_consensus.shape[1]
        base_params.gaps.useFixedPatterns = True
        base_params.gaps.fixedPatterns = pycogaps.Matrix(fixed_consensus)
        base_params.gaps.whichMatrixFixed = "A" if args.distributed_mode == "single-cell" else "P"

        logger.info("[RUN] Final stage | consensus_patterns=%s", fixed_consensus.shape[1])
        with multiprocessing.get_context("spawn").Pool(processes=pool_workers) as pool:
            final_stage_results = pool.map(run_subset_worker, final_tasks)

        final_stage_summary = pd.DataFrame([summarize_worker(entry) for entry in final_stage_results])
        final_stage_summary.to_csv(final_stage_summary_path, index=False)

        stitched = stitchTogether(
            [entry["result"] for entry in final_stage_results],
            [entry["result"] for entry in first_pass_results],
            base_params,
            explicit_sets,
            full_data.copy(),
        )
        result_adata = build_result_adata(full_data, stitched, final_stage_results)
        atomic_write_h5ad(result_adata, result_path)

        preprocessed_h5ad = Path(args.preprocessed_h5ad) if args.preprocessed_h5ad else None
        ifn_pattern, ifn_corr, top_genes = compute_ifn_pattern(
            result=result_adata,
            cogaps_input=full_data,
            preprocessed_cells_h5ad=preprocessed_h5ad,
            outdir=None,
            top_n=int(args.top_genes),
            stim_label=args.stim_label,
        )

        top_genes_table(result_adata, n_top=int(args.top_genes)).to_csv(top_genes_path, index=False)
        trace_dataframe(result_adata).to_csv(trace_path, index=False)

        diagnostics = {
            **summarize_result_uns(result_adata),
            "firstPassWorkers": int(len(first_pass_results)),
            "finalStageWorkers": int(len(final_stage_results)),
            "consensusClusters": int(len(matched["clustered"])),
        }
        atomic_write_json(diagnostics, diagnostics_path)

        runtime_sec = float(time.time() - t0)
        metrics = {
            "status": "ok",
            "language": "Python",
            "package": "PyCoGAPS",
            "distributedMode": args.distributed_mode,
            "nPatterns": int(args.k),
            "finalPatterns": int(result_adata.obs.shape[1]),
            "nIterations": int(args.n_iter),
            "seed": int(args.seed),
            "sparseOptimization": bool(sparse_opt),
            "nSets": int(args.n_sets),
            "cut": int(cut_value),
            "minNS": int(min_ns_value),
            "maxNS": int(max_ns_value),
            "poolWorkers": int(pool_workers),
            "explicitSetsJson": str(Path(args.explicit_sets_json)),
            "subsetSizeMin": int(min(len(x) for x in explicit_sets)),
            "subsetSizeMedian": float(np.median([len(x) for x in explicit_sets])),
            "subsetSizeMax": int(max(len(x) for x in explicit_sets)),
            "outputFrequency": int(args.output_frequency),
            "nSnapshots": int(args.n_snapshots),
            "snapshotPhase": args.snapshot_phase,
            "takePumpSamples": bool(args.take_pump_samples),
            "runtime_sec": runtime_sec,
            "ifn_pattern": ifn_pattern,
            "ifn_corr": float(ifn_corr),
            "top_genes": top_genes,
            "result_path": str(result_path),
            "consensus_path": str(consensus_path),
            "metrics_path": str(metrics_path),
            **diagnostics,
            "timestamp": datetime.now().isoformat(timespec="seconds"),
        }
        atomic_write_json(metrics, metrics_path)
        logger.info("[DONE] %s | OK | %.2f min | IFN %s (%.4f)", tag, runtime_sec / 60.0, ifn_pattern, ifn_corr)

    except Exception as exc:
        status = "failed"
        err_txt = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        runtime_sec = float(time.time() - t0)
        atomic_write_json(
            {
                "status": status,
                "language": "Python",
                "package": "PyCoGAPS",
                "distributedMode": args.distributed_mode,
                "nPatterns": int(args.k),
                "nIterations": int(args.n_iter),
                "seed": int(args.seed),
                "sparseOptimization": bool(sparse_opt),
                "nSets": int(args.n_sets),
                "runtime_sec": runtime_sec,
                "error": err_txt,
                "timestamp": datetime.now().isoformat(timespec="seconds"),
            },
            metrics_path,
        )
        logger.error("[ERROR] %s", err_txt)
        raise


if __name__ == "__main__":
    main()
