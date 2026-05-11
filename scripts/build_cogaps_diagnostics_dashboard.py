#!/usr/bin/env python3
"""Build diagnostic plots and summary tables for the CoGAPS comparison runs."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(".")
OUT = ROOT / "data" / "diagnostics_model_comparison"
FIG = OUT / "figures"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams["figure.dpi"] = 140
plt.rcParams["savefig.bbox"] = "tight"


RUNS: list[dict[str, Any]] = [
    {
        "id": "python_original",
        "label": "Original Python chosen",
        "language": "Python",
        "mode": "single-process",
        "diagnostics": "light",
        "metrics": "data/results_python/cogaps_K7_seed2_iter2000.metrics.json",
        "activity": "data/results_python/pattern_activity_by_celltype_condition.csv",
    },
    {
        "id": "r_sparse_chosen",
        "label": "R sparse chosen",
        "language": "R",
        "mode": "single-process",
        "diagnostics": "light",
        "metrics": "data/results_r/cogaps_K7_seed2_iter2000.metrics.json",
        "trace": "data/results_r/cogaps_K7_seed2_iter2000.trace.csv",
        "activity": "data/results_r/pattern_activity_by_celltype_condition.csv",
    },
    {
        "id": "r_sparse_instrumented",
        "label": "R sparse instrumented",
        "language": "R",
        "mode": "single-process",
        "diagnostics": "heavy",
        "metrics": "data/results_r_sparse_instrumented/cogaps_K7_seed2_iter2000.metrics.json",
        "trace": "data/results_r_sparse_instrumented/cogaps_K7_seed2_iter2000.trace.csv",
        "activity": "data/results_r_sparse_instrumented/pattern_activity_by_celltype_condition.csv",
        "snapshot": "data/results_r_sparse_instrumented/cogaps_K7_seed2_iter2000.snapshot_summary.csv",
        "uncertainty": "data/results_r_sparse_instrumented/cogaps_K7_seed2_iter2000.pattern_uncertainty_summary.csv",
        "correlations": "data/results_r_sparse_instrumented/pattern_correlations.csv",
    },
    {
        "id": "r_nosparse_instrumented",
        "label": "R no-sparse instrumented",
        "language": "R",
        "mode": "single-process",
        "diagnostics": "heavy",
        "metrics": "data/results_r_nosparse/cogaps_K7_seed2_iter2000.metrics.json",
        "trace": "data/results_r_nosparse/cogaps_K7_seed2_iter2000.trace.csv",
        "activity": "data/results_r_nosparse/pattern_activity_by_celltype_condition.csv",
        "snapshot": "data/results_r_nosparse/cogaps_K7_seed2_iter2000.snapshot_summary.csv",
        "uncertainty": "data/results_r_nosparse/cogaps_K7_seed2_iter2000.pattern_uncertainty_summary.csv",
        "correlations": "data/results_r_nosparse/pattern_correlations.csv",
    },
    {
        "id": "r_distributed_sparse_on",
        "label": "R distributed sparse-on",
        "language": "R",
        "mode": "distributed single-cell",
        "diagnostics": "distributed",
        "metrics": "data/results_r_distributed_sparse_on_full/cogaps_K7_seed2_iter2000_single_cell.metrics.json",
        "activity": "data/diagnostics_model_comparison/tables/r_distributed_sparse_on_activity_by_celltype_condition.csv",
        "consensus": "data/results_r_distributed_sparse_on_full/cogaps_K7_seed2_iter2000_single_cell.consensus_summary.csv",
        "workers": "data/results_r_distributed_sparse_on_full/cogaps_K7_seed2_iter2000_single_cell.first_pass_worker_summary.csv",
    },
    {
        "id": "r_distributed_sparse_off",
        "label": "R distributed sparse-off",
        "language": "R",
        "mode": "distributed single-cell",
        "diagnostics": "distributed",
        "metrics": "data/results_r_distributed_sparse_off_full/cogaps_K7_seed2_iter2000_single_cell.metrics.json",
        "activity": "data/diagnostics_model_comparison/tables/r_distributed_sparse_off_activity_by_celltype_condition.csv",
        "consensus": "data/results_r_distributed_sparse_off_full/cogaps_K7_seed2_iter2000_single_cell.consensus_summary.csv",
        "workers": "data/results_r_distributed_sparse_off_full/cogaps_K7_seed2_iter2000_single_cell.first_pass_worker_summary.csv",
    },
    {
        "id": "python_sparse_light_local",
        "label": "Python single sparse light-local",
        "language": "Python",
        "mode": "single-process",
        "diagnostics": "light-local",
        "metrics": "data/results_python_sparse_light_local/runs/cogaps_K7_seed2_iter2000.metrics.json",
        "trace": "data/results_python_sparse_light_local/runs/cogaps_K7_seed2_iter2000.trace.csv",
        "h5ad": "data/results_python_sparse_light_local/runs/cogaps_K7_seed2_iter2000.h5ad",
    },
    {
        "id": "python_distributed_sparse_on_light",
        "label": "Python distributed sparse-on light",
        "language": "Python",
        "mode": "distributed single-cell",
        "diagnostics": "distributed light",
        "metrics": "data/results_python_distributed_sparse_on_light_full/cogaps_K7_seed2_iter2000_single_cell.metrics.json",
        "trace": "data/results_python_distributed_sparse_on_light_full/cogaps_K7_seed2_iter2000_single_cell.trace.csv",
        "h5ad": "data/results_python_distributed_sparse_on_light_full/cogaps_K7_seed2_iter2000_single_cell.h5ad",
        "consensus": "data/results_python_distributed_sparse_on_light_full/cogaps_K7_seed2_iter2000_single_cell.consensus_summary.csv",
        "workers": "data/results_python_distributed_sparse_on_light_full/cogaps_K7_seed2_iter2000_single_cell.first_pass_worker_summary.csv",
        "final_workers": "data/results_python_distributed_sparse_on_light_full/cogaps_K7_seed2_iter2000_single_cell.final_stage_worker_summary.csv",
    },
]


def path_or_none(value: str | None) -> Path | None:
    if not value:
        return None
    p = Path(value)
    return p if p.exists() else None


def read_json(path: str | None) -> dict[str, Any]:
    p = path_or_none(path)
    if p is None:
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def pattern_sort_key(pattern: str) -> int:
    return int(str(pattern).replace("Pattern", "").replace("_", ""))


def normalize_trace(path: Path, run: dict[str, Any]) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.strip('"') for c in df.columns]
    if "relative_progress" not in df.columns:
        df["relative_progress"] = df["trace_index"] / df["trace_index"].max()
    df["phase"] = np.where(df["relative_progress"] <= 0.5, "equilibration", "sampling")
    df["run_id"] = run["id"]
    df["run"] = run["label"]
    df["mode"] = run["mode"]
    df["language"] = run["language"]
    df["delta_chisq"] = df["chisq"].diff()
    return df


def trace_stats(trace: pd.DataFrame) -> dict[str, Any]:
    chisq = trace["chisq"].astype(float)
    sampling = trace.loc[trace["phase"] == "sampling", "chisq"].astype(float)
    last5 = chisq.tail(5).to_numpy(dtype=float)
    if len(last5) > 1:
        slope = np.polyfit(np.arange(len(last5)), last5, 1)[0]
    else:
        slope = np.nan
    return {
        "trace_points": int(len(trace)),
        "chisq_start": float(chisq.iloc[0]),
        "chisq_end": float(chisq.iloc[-1]),
        "chisq_min": float(chisq.min()),
        "chisq_drop_pct": float((chisq.iloc[0] - chisq.iloc[-1]) / chisq.iloc[0] * 100),
        "chisq_increases": int((chisq.diff() > 0).sum()),
        "sampling_range_pct_of_end": float((sampling.max() - sampling.min()) / chisq.iloc[-1] * 100) if len(sampling) else np.nan,
        "last5_slope_pct_of_end": float(slope / chisq.iloc[-1] * 100) if len(last5) > 1 else np.nan,
        "atomsA_end": float(trace["atomsA"].iloc[-1]) if "atomsA" in trace else np.nan,
        "atomsP_end": float(trace["atomsP"].iloc[-1]) if "atomsP" in trace else np.nan,
    }


def pattern_columns(columns: pd.Index) -> list[str]:
    pats = [c for c in columns if str(c).startswith("Pattern")]
    return sorted(pats, key=pattern_sort_key)


def activity_from_h5ad(path: str, run: dict[str, Any]) -> pd.DataFrame:
    meta = ad.read_h5ad("data/processed/preprocessed_cells_hvg3000.h5ad", backed="r")
    obs = meta.obs[["cell_type", "condition", "replicate"]].copy()
    meta.file.close()

    result = ad.read_h5ad(path, backed="r")
    pats = pattern_columns(result.var.columns)
    P = result.var[pats].copy()
    result.file.close()

    df = obs.join(P, how="inner")
    rows = []
    for pat in pats:
        g = (
            df.groupby(["cell_type", "condition"], observed=True)[pat]
            .agg(mean_activity="mean", median_activity="median", n_cells="count")
            .reset_index()
        )
        g.insert(0, "pattern", pat)
        rows.append(g)
    out = pd.concat(rows, ignore_index=True)
    out["run_id"] = run["id"]
    out["run"] = run["label"]
    return out


metrics_rows: list[dict[str, Any]] = []
traces: list[pd.DataFrame] = []
activities: list[pd.DataFrame] = []
consensus_frames: list[pd.DataFrame] = []
worker_frames: list[pd.DataFrame] = []
snapshot_frames: list[pd.DataFrame] = []
uncertainty_frames: list[pd.DataFrame] = []
correlation_frames: list[pd.DataFrame] = []

for run in RUNS:
    metrics = read_json(run.get("metrics"))
    row: dict[str, Any] = {
        "run_id": run["id"],
        "run": run["label"],
        "language": run["language"],
        "mode": run["mode"],
        "diagnostics": run["diagnostics"],
        "status": metrics.get("status"),
        "sparseOptimization": metrics.get("sparseOptimization"),
        "runtime_min": metrics.get("runtime_sec", np.nan) / 60 if metrics.get("runtime_sec") else np.nan,
        "meanChiSq": metrics.get("meanChiSq"),
        "totalRunningTime_sec": metrics.get("totalRunningTime"),
        "totalUpdates": metrics.get("totalUpdates"),
        "ifn_pattern": metrics.get("ifn_pattern"),
        "ifn_corr": metrics.get("ifn_corr"),
        "nPatterns_requested": metrics.get("nPatterns") or metrics.get("K"),
        "seed_recorded": metrics.get("seed"),
        "outputFrequency": metrics.get("outputFrequency"),
        "checkpointInterval": metrics.get("checkpointInterval"),
        "nSnapshots": metrics.get("nSnapshots"),
        "takePumpSamples": metrics.get("takePumpSamples"),
    }

    trace_path = path_or_none(run.get("trace"))
    if trace_path is not None:
        trace = normalize_trace(trace_path, run)
        traces.append(trace)
        row.update(trace_stats(trace))

    act_path = path_or_none(run.get("activity"))
    if act_path is not None:
        act = pd.read_csv(act_path)
        act["run_id"] = run["id"]
        act["run"] = run["label"]
        activities.append(act)
    elif run.get("h5ad") and path_or_none(run.get("h5ad")) is not None:
        activities.append(activity_from_h5ad(run["h5ad"], run))

    consensus_path = path_or_none(run.get("consensus"))
    if consensus_path is not None:
        c = pd.read_csv(consensus_path)
        c["run_id"] = run["id"]
        c["run"] = run["label"]
        consensus_frames.append(c)

    for worker_key, stage in [("workers", "first_pass"), ("final_workers", "final_stage")]:
        worker_path = path_or_none(run.get(worker_key))
        if worker_path is not None:
            w = pd.read_csv(worker_path)
            w["run_id"] = run["id"]
            w["run"] = run["label"]
            w["stage"] = stage
            worker_frames.append(w)

    snapshot_path = path_or_none(run.get("snapshot"))
    if snapshot_path is not None:
        s = pd.read_csv(snapshot_path)
        s["run_id"] = run["id"]
        s["run"] = run["label"]
        snapshot_frames.append(s)

    uncertainty_path = path_or_none(run.get("uncertainty"))
    if uncertainty_path is not None:
        u = pd.read_csv(uncertainty_path)
        u["run_id"] = run["id"]
        u["run"] = run["label"]
        uncertainty_frames.append(u)

    corr_path = path_or_none(run.get("correlations"))
    if corr_path is not None:
        corr = pd.read_csv(corr_path)
        corr["run_id"] = run["id"]
        corr["run"] = run["label"]
        correlation_frames.append(corr)

    metrics_rows.append(row)

metrics_df = pd.DataFrame(metrics_rows)
metrics_df.to_csv(TAB / "model_diagnostic_summary.csv", index=False)

trace_df = pd.concat(traces, ignore_index=True) if traces else pd.DataFrame()
if not trace_df.empty:
    trace_df.to_csv(TAB / "combined_trace.csv", index=False)

activity_df = pd.concat(activities, ignore_index=True) if activities else pd.DataFrame()
if not activity_df.empty:
    activity_df.to_csv(TAB / "combined_activity_by_celltype_condition.csv", index=False)

consensus_df = pd.concat(consensus_frames, ignore_index=True) if consensus_frames else pd.DataFrame()
if not consensus_df.empty:
    consensus_df.to_csv(TAB / "combined_consensus_summary.csv", index=False)

worker_df = pd.concat(worker_frames, ignore_index=True) if worker_frames else pd.DataFrame()
if not worker_df.empty:
    worker_df.to_csv(TAB / "combined_worker_summary.csv", index=False)

snapshot_df = pd.concat(snapshot_frames, ignore_index=True) if snapshot_frames else pd.DataFrame()
if not snapshot_df.empty:
    snapshot_df.to_csv(TAB / "combined_snapshot_summary.csv", index=False)

uncertainty_df = pd.concat(uncertainty_frames, ignore_index=True) if uncertainty_frames else pd.DataFrame()
if not uncertainty_df.empty:
    uncertainty_df.to_csv(TAB / "combined_uncertainty_summary.csv", index=False)

correlation_df = pd.concat(correlation_frames, ignore_index=True) if correlation_frames else pd.DataFrame()
if not correlation_df.empty:
    correlation_df.to_csv(TAB / "combined_pattern_correlations.csv", index=False)


def savefig(name: str) -> str:
    path = FIG / name
    plt.savefig(path)
    plt.close()
    return f"figures/{name}"


fig_refs: list[tuple[str, str, str]] = []

plt.figure(figsize=(14, 7))
plot_df = metrics_df.dropna(subset=["runtime_min"]).copy()
plot_df["run_short"] = plot_df["run"].str.replace(" instrumented", " instr.", regex=False).str.replace(" light-local", " light", regex=False)
ax = sns.barplot(data=plot_df, x="runtime_min", y="run_short", hue="language", dodge=False)
ax.set_xlabel("Runtime (minutes)")
ax.set_ylabel("")
ax.set_title("Wall-clock runtime across saved CoGAPS runs")
for container in ax.containers:
    ax.bar_label(container, fmt="%.1f", padding=3, fontsize=9)
fig_refs.append(("Runtime", savefig("runtime_summary.png"), "Runtime separates biological reproducibility from performance. The current local Python sparse run remains slow even with light diagnostics."))

if not trace_df.empty:
    single_trace = trace_df[trace_df["mode"] == "single-process"].copy()
    plt.figure(figsize=(14, 7))
    ax = sns.lineplot(data=single_trace, x="relative_progress", y="chisq", hue="run", marker="o")
    ax.axvline(0.5, color="black", linestyle="--", linewidth=1)
    ax.set_title("Single-process chi-square traces")
    ax.set_xlabel("Relative progress (0.5 = end of equilibration)")
    ax.set_ylabel("Chi-square")
    fig_refs.append(("Single-Process Chi-Square", savefig("singleprocess_chisq_trace.png"), "Good traces drop steeply during equilibration and then fluctuate on a narrow plateau during sampling."))

    plt.figure(figsize=(14, 7))
    sampling_trace = single_trace[single_trace["phase"] == "sampling"].copy()
    ax = sns.lineplot(data=sampling_trace, x="relative_progress", y="chisq", hue="run", marker="o")
    ax.set_title("Sampling-phase chi-square zoom")
    ax.set_xlabel("Relative progress")
    ax.set_ylabel("Chi-square")
    fig_refs.append(("Sampling Zoom", savefig("sampling_phase_chisq_zoom.png"), "Sampling-phase wiggle is expected; persistent drift would be more concerning than small oscillation."))

    atom_long = single_trace.melt(
        id_vars=["run", "relative_progress", "phase"],
        value_vars=[c for c in ["atomsA", "atomsP"] if c in single_trace.columns],
        var_name="atom_matrix",
        value_name="atoms",
    )
    plt.figure(figsize=(14, 7))
    ax = sns.lineplot(data=atom_long, x="relative_progress", y="atoms", hue="run", style="atom_matrix")
    ax.axvline(0.5, color="black", linestyle="--", linewidth=1)
    ax.set_title("Single-process atom count traces")
    ax.set_xlabel("Relative progress")
    ax.set_ylabel("Active atoms")
    fig_refs.append(("Atom Counts", savefig("singleprocess_atom_traces.png"), "Atom counts show how complex the atomic representation becomes; late stabilization supports a settled posterior region."))

    dist_trace = trace_df[trace_df["mode"].str.contains("distributed", na=False)].copy()
    if not dist_trace.empty:
        plt.figure(figsize=(12, 6))
        ax = sns.lineplot(data=dist_trace, x="relative_progress", y="chisq", hue="run", marker="o")
        ax.axvline(0.5, color="black", linestyle="--", linewidth=1)
        ax.set_title("Distributed final-stage trace where available")
        ax.set_xlabel("Relative progress")
        ax.set_ylabel("Chi-square")
        fig_refs.append(("Distributed Trace", savefig("distributed_trace_available.png"), "Distributed traces are final-stage/worker-scale diagnostics, not directly comparable to full single-process chi-square values."))

plt.figure(figsize=(14, 7))
stab_cols = ["run", "chisq_drop_pct", "sampling_range_pct_of_end", "last5_slope_pct_of_end", "chisq_increases"]
stab = metrics_df.dropna(subset=["chisq_drop_pct"])[stab_cols].copy()
stab_long = stab.melt(id_vars=["run"], var_name="metric", value_name="value")
ax = sns.barplot(data=stab_long, x="value", y="run", hue="metric")
ax.set_title("Trace stability metrics")
ax.set_xlabel("Metric value")
ax.set_ylabel("")
fig_refs.append(("Trace Stability Metrics", savefig("trace_stability_metrics.png"), "Lower sampling range and near-zero late slope are reassuring; chi-square increases are normal in MCMC sampling."))

if not snapshot_df.empty:
    snap = snapshot_df.copy()
    snap["snapshot_iter"] = pd.to_numeric(snap["snapshot_label"], errors="coerce")
    plt.figure(figsize=(14, 8))
    ax = sns.lineplot(
        data=snap,
        x="snapshot_iter",
        y="mean_top_overlap",
        hue="run",
        style="matrix",
        markers=True,
        dashes=False,
    )
    ax.set_title("Snapshot top-gene overlap to final result")
    ax.set_xlabel("Snapshot iteration within phase")
    ax.set_ylabel("Mean top-gene overlap to final")
    fig_refs.append(("Snapshot Top-Gene Stability", savefig("snapshot_top_gene_overlap.png"), "Increasing overlap indicates that top genes stabilize as the sampler approaches and samples the posterior region."))

    plt.figure(figsize=(14, 8))
    ax = sns.lineplot(
        data=snap,
        x="snapshot_iter",
        y="correlation_to_final",
        hue="run",
        style="matrix",
        markers=True,
        dashes=False,
    )
    ax.set_title("Snapshot correlation to final result")
    ax.set_xlabel("Snapshot iteration within phase")
    ax.set_ylabel("Correlation to final matrix")
    fig_refs.append(("Snapshot Correlation", savefig("snapshot_correlation_to_final.png"), "Correlation-to-final is a broad stability check; low early values are expected, late instability is not."))

if not uncertainty_df.empty:
    unc = uncertainty_df.copy()
    unc["pattern_num"] = unc["pattern"].map(pattern_sort_key)
    plt.figure(figsize=(14, 7))
    ax = sns.barplot(data=unc, x="pattern", y="top_gene_median_loading_cv", hue="run")
    ax.set_title("Top-gene loading coefficient of variation")
    ax.set_xlabel("Pattern")
    ax.set_ylabel("Median loading CV among top genes")
    fig_refs.append(("Top-Gene Uncertainty", savefig("top_gene_loading_cv.png"), "Lower top-gene CV suggests more stable gene loading estimates for the program."))

    plt.figure(figsize=(14, 7))
    ax = sns.barplot(data=unc, x="pattern", y="mean_factor_sd", hue="run")
    ax.set_title("Mean sample-factor posterior SD")
    ax.set_xlabel("Pattern")
    ax.set_ylabel("Mean factor SD")
    fig_refs.append(("Activity Uncertainty", savefig("mean_factor_sd.png"), "Higher factor SD means more posterior uncertainty in cell/sample activity for that pattern."))

if not consensus_df.empty:
    cons = consensus_df.copy()
    cons["cluster"] = cons["cluster"].astype(str)
    cons_long = cons.melt(
        id_vars=["run", "cluster", "n_members"],
        value_vars=["min_corr_to_consensus", "mean_corr_to_consensus"],
        var_name="correlation_type",
        value_name="correlation",
    )
    plt.figure(figsize=(14, 7))
    ax = sns.barplot(data=cons_long, x="cluster", y="correlation", hue="run")
    ax.set_ylim(0.75, 1.01)
    ax.set_title("Consensus cluster correlation minima and means")
    ax.set_xlabel("Consensus cluster")
    ax.set_ylabel("Correlation to consensus")
    fig_refs.append(("Consensus Stability", savefig("consensus_cluster_correlations.png"), "Clusters with low minimum correlation or unusual member counts are the distributed patterns to inspect most carefully."))

if not worker_df.empty:
    worker = worker_df.copy()
    worker["worker_id"] = worker["worker_id"].astype(str)
    plt.figure(figsize=(14, 7))
    ax = sns.barplot(data=worker, x="worker_id", y="totalRunningTime", hue="run")
    ax.set_title("Distributed worker running times")
    ax.set_xlabel("Worker")
    ax.set_ylabel("CoGAPS-reported worker runtime (sec)")
    fig_refs.append(("Worker Runtime", savefig("worker_runtime.png"), "Distributed runtime is limited by the slowest worker plus matching/final-stage overhead."))

if not activity_df.empty:
    activity_df["pattern"] = activity_df["pattern"].astype(str).str.replace("Pattern_", "Pattern", regex=False)
    ifn_map = metrics_df.set_index("run_id")["ifn_pattern"].to_dict()
    ifn_activity = activity_df[activity_df.apply(lambda r: r["pattern"] == ifn_map.get(r["run_id"]), axis=1)].copy()
    ifn_activity.to_csv(TAB / "ifn_activity_by_celltype_condition.csv", index=False)
    cell_order = (
        ifn_activity.groupby("cell_type")["mean_activity"].max().sort_values(ascending=False).index.tolist()
    )
    n_runs = ifn_activity["run"].nunique()
    fig, axes = plt.subplots(n_runs, 1, figsize=(16, max(4, 3.4 * n_runs)), sharex=True)
    if n_runs == 1:
        axes = [axes]
    for ax, (run_name, sub) in zip(axes, ifn_activity.groupby("run", sort=False)):
        sns.barplot(data=sub, x="cell_type", y="mean_activity", hue="condition", order=cell_order, ax=ax)
        ax.set_title(run_name)
        ax.set_xlabel("")
        ax.set_ylabel("Mean IFN-pattern activity")
        ax.tick_params(axis="x", rotation=35)
    fig.suptitle("IFN-pattern activity by cell type and condition", y=1.01)
    fig_refs.append(("IFN Activity", savefig("ifn_activity_by_celltype_condition.png"), "Effective biological recovery should show strong stimulation-associated IFN activity across the expected immune cell groups, not just a numerical fit."))

    delta_rows = []
    for (run_id, run_name, pattern, cell_type), g in ifn_activity.groupby(["run_id", "run", "pattern", "cell_type"]):
        vals = dict(zip(g["condition"], g["mean_activity"]))
        if "stim" in vals and "ctrl" in vals:
            delta_rows.append(
                {
                    "run_id": run_id,
                    "run": run_name,
                    "pattern": pattern,
                    "cell_type": cell_type,
                    "stim_minus_ctrl": vals["stim"] - vals["ctrl"],
                }
            )
    delta = pd.DataFrame(delta_rows)
    delta.to_csv(TAB / "ifn_activity_stim_minus_ctrl.csv", index=False)
    if not delta.empty:
        pivot = delta.pivot(index="run", columns="cell_type", values="stim_minus_ctrl").fillna(0)
        pivot = pivot[cell_order]
        plt.figure(figsize=(16, max(5, 0.55 * len(pivot))))
        ax = sns.heatmap(pivot, cmap="vlag", center=0, linewidths=0.5)
        ax.set_title("IFN-pattern stim-minus-control mean activity")
        ax.set_xlabel("Cell type")
        ax.set_ylabel("")
        fig_refs.append(("IFN Stim-Control Delta", savefig("ifn_activity_stim_ctrl_delta_heatmap.png"), "This heatmap shows whether the IFN program is condition-responsive in the same biological direction across implementations."))

if not correlation_df.empty:
    plt.figure(figsize=(12, 6))
    ax = sns.barplot(data=correlation_df, x="pattern", y="corr_with_stim", hue="run")
    ax.axhline(0, color="black", linewidth=1)
    ax.set_title("Pattern correlations with stimulation")
    ax.set_xlabel("Pattern")
    ax.set_ylabel("Correlation with stim indicator")
    fig_refs.append(("Pattern-Stim Correlations", savefig("pattern_stim_correlations.png"), "This identifies whether the IFN-associated program is uniquely dominant or whether multiple programs track condition."))


def table_html(df: pd.DataFrame, max_rows: int = 20) -> str:
    return df.head(max_rows).to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}")


summary_cols = [
    "run",
    "language",
    "mode",
    "diagnostics",
    "runtime_min",
    "sparseOptimization",
    "ifn_pattern",
    "ifn_corr",
    "meanChiSq",
    "chisq_drop_pct",
    "sampling_range_pct_of_end",
    "last5_slope_pct_of_end",
]

html_parts = [
    "<!doctype html><html><head><meta charset='utf-8'><title>CoGAPS Model Diagnostics</title>",
    "<style>body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;max-width:1200px;margin:32px auto;line-height:1.45} img{max-width:100%;border:1px solid #ddd} table{border-collapse:collapse;font-size:13px} th,td{border:1px solid #ddd;padding:4px 6px} th{background:#f4f4f4} .note{background:#fff8df;padding:12px;border-left:4px solid #d6a100}</style>",
    "</head><body>",
    "<h1>CoGAPS Model Diagnostics</h1>",
    "<p class='note'>TensorBoard was not used here because these are completed tabular MCMC diagnostics rather than live training-event streams, and TensorBoard is not installed in the runtime. This report uses direct trace, posterior-summary, consensus, and activity diagnostics from the saved artifacts.</p>",
    "<h2>Diagnostic Summary Table</h2>",
    table_html(metrics_df[summary_cols].sort_values(["mode", "runtime_min"], na_position="last"), max_rows=30),
]

for title, rel, caption in fig_refs:
    html_parts.extend([
        f"<h2>{html.escape(title)}</h2>",
        f"<p>{html.escape(caption)}</p>",
        f"<img src='{html.escape(rel)}' alt='{html.escape(title)}'>",
    ])

html_parts.extend(
    [
        "<h2>Generated Tables</h2>",
        "<ul>",
        "<li><code>tables/model_diagnostic_summary.csv</code></li>",
        "<li><code>tables/combined_trace.csv</code></li>",
        "<li><code>tables/combined_activity_by_celltype_condition.csv</code></li>",
        "<li><code>tables/combined_consensus_summary.csv</code></li>",
        "<li><code>tables/combined_worker_summary.csv</code></li>",
        "<li><code>tables/ifn_activity_by_celltype_condition.csv</code></li>",
        "<li><code>tables/ifn_activity_stim_minus_ctrl.csv</code></li>",
        "</ul>",
        "</body></html>",
    ]
)

(OUT / "diagnostic_dashboard.html").write_text("\n".join(html_parts), encoding="utf-8")

md_parts = [
    "# CoGAPS Model Diagnostics",
    "",
    "This folder contains trace, stability, consensus, and activity diagnostics for the saved R and Python CoGAPS comparison runs.",
    "",
    "TensorBoard was considered but not used: the available artifacts are completed CSV/JSON/RDS/H5AD diagnostics, not live scalar event streams, and the runtime does not include TensorBoard. The HTML dashboard is a lower-dependency fit for this evidence.",
    "",
    "## Key Files",
    "",
    "- `diagnostic_dashboard.html`: browseable report with all figures.",
    "- `tables/model_diagnostic_summary.csv`: runtime, fit, trace-stability, and IFN summary metrics.",
    "- `tables/combined_trace.csv`: harmonized chi-square and atom traces.",
    "- `tables/combined_activity_by_celltype_condition.csv`: pattern activity summaries.",
    "- `tables/combined_consensus_summary.csv`: distributed consensus diagnostics.",
    "",
    "## Figures",
]
for title, rel, caption in fig_refs:
    md_parts.extend([f"### {title}", "", caption, "", f"![{title}]({rel})", ""])
(OUT / "README.md").write_text("\n".join(md_parts), encoding="utf-8")

print(f"Wrote {OUT / 'diagnostic_dashboard.html'}")
print(f"Wrote {OUT / 'README.md'}")
print(f"Wrote {TAB / 'model_diagnostic_summary.csv'}")
