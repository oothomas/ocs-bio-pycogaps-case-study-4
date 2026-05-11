#!/usr/bin/env python3
"""Write a results-driven CoGAPS diagnostic HTML/Markdown report."""

from __future__ import annotations

import csv
import html
from pathlib import Path


OUT = Path("data/diagnostics_model_comparison")
TAB = OUT / "tables"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def as_float(value: str | None) -> float | None:
    if value is None or value == "" or value == "NA" or value == "nan":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def fmt(value: float | None, digits: int = 3) -> str:
    if value is None:
        return "NA"
    return f"{value:,.{digits}f}"


def fmt_int(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{int(round(value)):,}"


summary = read_csv(TAB / "model_diagnostic_summary.csv")
chisq = {r["run"]: r for r in read_csv(TAB / "singleprocess_chisq_interpretive_summary.csv")}
consensus = read_csv(TAB / "combined_consensus_summary.csv")
ifn_deltas = read_csv(TAB / "ifn_activity_stim_minus_ctrl.csv")
uncertainty = read_csv(TAB / "combined_uncertainty_summary.csv")

by_run = {r["run"]: r for r in summary}


def row(run: str, key: str) -> str:
    return by_run.get(run, {}).get(key, "")


def val(run: str, key: str) -> float | None:
    return as_float(row(run, key))


def consensus_stats(run: str) -> dict[str, float | int | None]:
    rows = [r for r in consensus if r["run"] == run]
    if not rows:
        return {"clusters": None, "min_corr": None, "mean_min_corr": None, "weak_clusters": None}
    mins = [as_float(r["min_corr_to_consensus"]) for r in rows]
    means = [as_float(r["mean_corr_to_consensus"]) for r in rows]
    mins = [x for x in mins if x is not None]
    means = [x for x in means if x is not None]
    return {
        "clusters": len(rows),
        "min_corr": min(mins) if mins else None,
        "mean_min_corr": sum(mins) / len(mins) if mins else None,
        "weak_clusters": sum(1 for x in mins if x < 0.9),
        "mean_corr": sum(means) / len(means) if means else None,
    }


def mean_ifn_delta(run: str) -> float | None:
    vals = [as_float(r["stim_minus_ctrl"]) for r in ifn_deltas if r["run"] == run]
    vals = [x for x in vals if x is not None]
    return sum(vals) / len(vals) if vals else None


def ifn_delta_range(run: str) -> tuple[float | None, float | None]:
    vals = [as_float(r["stim_minus_ctrl"]) for r in ifn_deltas if r["run"] == run]
    vals = [x for x in vals if x is not None]
    return (min(vals), max(vals)) if vals else (None, None)


def ifn_uncertainty(run: str, pattern: str) -> tuple[float | None, float | None]:
    rows = [r for r in uncertainty if r["run"] == run and r["pattern"] == pattern]
    if not rows:
        return None, None
    r = rows[0]
    return as_float(r.get("top_gene_median_loading_cv")), as_float(r.get("mean_factor_sd"))


sparse_mean = val("R sparse instrumented", "meanChiSq")
nosparse_mean = val("R no-sparse instrumented", "meanChiSq")
mean_delta_pct = None
if sparse_mean is not None and nosparse_mean is not None:
    mean_delta_pct = (nosparse_mean - sparse_mean) / sparse_mean * 100

runtime_ratio_nosparse = None
if val("R sparse instrumented", "runtime_min") and val("R no-sparse instrumented", "runtime_min"):
    runtime_ratio_nosparse = val("R no-sparse instrumented", "runtime_min") / val("R sparse instrumented", "runtime_min")

runtime_ratio_python = None
if val("Original Python chosen", "runtime_min") and val("Python single sparse light-local", "runtime_min"):
    runtime_ratio_python = val("Python single sparse light-local", "runtime_min") / val("Original Python chosen", "runtime_min")


decision_rows = [
    {
        "question": "Best diagnostic anchor for final narrative?",
        "answer": "R sparse instrumented",
        "why": "Full 40-point trace, snapshots, pump/posterior uncertainty, stable IFN program, and fast runtime relative to no-sparse.",
    },
    {
        "question": "Does no-sparse improve fit enough to justify runtime?",
        "answer": "No",
        "why": f"Posterior meanChiSq is slightly worse than sparse by {fmt(mean_delta_pct, 4)}%, while runtime is {fmt(runtime_ratio_nosparse, 2)}x longer.",
    },
    {
        "question": "Does Python reproduce the sparse biology?",
        "answer": "Yes",
        "why": "Python single sparse light-local has the same trace, IFN pattern, IFN correlation, top genes, total updates, and meanChiSq as R sparse.",
    },
    {
        "question": "Does current local Python reproduce original 39 minute performance?",
        "answer": "No",
        "why": f"The light-local rerun took {fmt(runtime_ratio_python, 2)}x longer than the original saved Python run, even with pump/checkpoints/snapshots off.",
    },
    {
        "question": "Best acceleration candidate?",
        "answer": "R distributed sparse-on",
        "why": "Fastest completed full run, preserves the IFN program, but has weaker consensus for at least one secondary cluster.",
    },
]


model_rows = [
    "Original Python chosen",
    "R sparse instrumented",
    "R no-sparse instrumented",
    "R distributed sparse-on",
    "R distributed sparse-off",
    "Python single sparse light-local",
    "Python distributed sparse-on light",
]


def table(headers: list[str], rows: list[list[str]]) -> str:
    out = ["<table>", "<thead><tr>"]
    out.extend(f"<th>{html.escape(h)}</th>" for h in headers)
    out.append("</tr></thead><tbody>")
    for r in rows:
        out.append("<tr>")
        out.extend(f"<td>{html.escape(str(c))}</td>" for c in r)
        out.append("</tr>")
    out.append("</tbody></table>")
    return "\n".join(out)


model_table_rows = []
for run in model_rows:
    cstats = consensus_stats(run)
    model_table_rows.append(
        [
            run,
            row(run, "mode"),
            row(run, "sparseOptimization") or "NA",
            fmt(val(run, "runtime_min"), 2),
            row(run, "ifn_pattern") or "NA",
            fmt(val(run, "ifn_corr"), 4),
            fmt_int(val(run, "meanChiSq")),
            fmt(cstats["min_corr"], 3) if cstats["min_corr"] is not None else "NA",
            fmt(mean_ifn_delta(run), 3),
            row(run, "seed_recorded") or "NA",
        ]
    )

decision_table_rows = [[d["question"], d["answer"], d["why"]] for d in decision_rows]

chisq_table_rows = []
for run in ["R sparse instrumented", "R no-sparse instrumented", "Python single sparse light-local"]:
    r = chisq[run]
    chisq_table_rows.append(
        [
            run,
            fmt_int(as_float(r["final_logged_chisq"])),
            fmt_int(as_float(r["posterior_meanChiSq"])),
            fmt_int(as_float(r["sampling_range"])),
            fmt(as_float(r["runtime_min"]), 2),
        ]
    )

consensus_rows = []
for run in ["R distributed sparse-on", "R distributed sparse-off", "Python distributed sparse-on light"]:
    cstats = consensus_stats(run)
    consensus_rows.append(
        [
            run,
            cstats["clusters"],
            fmt(cstats["min_corr"], 3),
            fmt(cstats["mean_min_corr"], 3),
            cstats["weak_clusters"],
            fmt(cstats["mean_corr"], 3),
        ]
    )

ifn_rows = []
for run in model_rows:
    lo, hi = ifn_delta_range(run)
    ifn_rows.append([run, row(run, "ifn_pattern") or "NA", fmt(mean_ifn_delta(run), 3), fmt(lo, 3), fmt(hi, 3)])

uncertainty_rows = []
for run, pattern in [
    ("R sparse instrumented", "Pattern2"),
    ("R no-sparse instrumented", "Pattern2"),
]:
    cv, factor_sd = ifn_uncertainty(run, pattern)
    uncertainty_rows.append([run, pattern, fmt(cv, 4), fmt(factor_sd, 4)])


figures = [
    ("Interpretive Chi-Square Trace", "figures/singleprocess_chisq_trace_interpretive.png", "Use this instead of the earlier single-axis plot. It separates the full burn-in behavior from the sampling-phase zoom and adds posterior meanChiSq reference lines."),
    ("Sampling Trace Minus Posterior meanChiSq", "figures/singleprocess_chisq_minus_mean_sampling.png", "This shows why the green no-sparse trace looking lower is not the same as having a better posterior-summary fit."),
    ("Runtime Summary", "figures/runtime_summary.png", "Performance is part of model suitability here because no-sparse and local Python sparse impose large costs without improving the biological answer."),
    ("Trace Stability Metrics", "figures/trace_stability_metrics.png", "Use these as secondary checks: plateau behavior, sampling range, and late slope matter more than tiny differences in logged current-state chi-square."),
    ("Snapshot Top-Gene Stability", "figures/snapshot_top_gene_overlap.png", "For the instrumented R runs, this evaluates whether top genes stabilize during sampling."),
    ("Top-Gene Loading CV", "figures/top_gene_loading_cv.png", "Lower coefficient of variation means the leading genes for a pattern are more stable across posterior samples."),
    ("Consensus Cluster Correlations", "figures/consensus_cluster_correlations.png", "For distributed runs, this is the key stability diagnostic. Weak consensus clusters point to secondary-pattern instability."),
    ("IFN Activity Delta", "figures/ifn_activity_stim_ctrl_delta_heatmap.png", "This is the most direct biological suitability check: the IFN pattern should increase in stimulated cells across expected cell types."),
]


html_parts = [
    "<!doctype html><html><head><meta charset='utf-8'><title>CoGAPS Diagnostics: Interpreted Report</title>",
    "<style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;max-width:1180px;margin:32px auto;line-height:1.48;color:#222}",
    "h1,h2,h3{line-height:1.15} img{max-width:100%;border:1px solid #ddd;margin:8px 0 20px}",
    "table{border-collapse:collapse;font-size:13px;margin:12px 0 24px;width:100%} th,td{border:1px solid #ddd;padding:6px 8px;vertical-align:top} th{background:#f4f4f4;text-align:left}",
    ".callout{background:#eef7ff;border-left:5px solid #2474b5;padding:12px 16px;margin:16px 0}",
    ".warn{background:#fff8df;border-left:5px solid #d6a100;padding:12px 16px;margin:16px 0}",
    ".bad{background:#fff1f1;border-left:5px solid #c74343;padding:12px 16px;margin:16px 0}",
    "code{background:#f5f5f5;padding:1px 4px}",
    "</style></head><body>",
    "<h1>CoGAPS Diagnostics: Interpreted Report</h1>",
    "<div class='callout'><strong>Bottom line.</strong> The sparse single-process model is the best diagnostic anchor. R distributed sparse-on is the best acceleration candidate. R no-sparse is not meaningfully a better fit, and current local Python is biologically reproducible but not performance-competitive.</div>",
    "<h2>Decision Summary</h2>",
    table(["Question", "Answer", "Why"], decision_table_rows),
    "<h2>Model-Level Evidence</h2>",
    table(["Run", "Mode", "Sparse", "Runtime min", "IFN pattern", "IFN corr", "meanChiSq", "Min consensus corr", "Mean IFN stim-control delta", "Recorded seed"], model_table_rows),
    "<div class='warn'><strong>Seed caveat.</strong> PyCoGAPS outputs record seed <code>0</code> for the reruns even though the launcher passed seed 2. Because the trace, top genes, total updates, and meanChiSq match the R sparse run exactly, this looks like a PyCoGAPS metadata/reporting issue, but it should be noted.</div>",
    "<h2>What The Chi-Square Results Actually Mean</h2>",
    "<p>The logged trace is a current-state sampler diagnostic. The final CoGAPS result is a posterior summary accumulated during sampling, so the last trace point is not the whole model-quality story.</p>",
    table(["Run", "Final logged chi-square", "Posterior meanChiSq", "Sampling range", "Runtime min"], chisq_table_rows),
    f"<p>The no-sparse trace sits lower in the logged current-state plot, but its posterior <code>meanChiSq</code> is slightly higher than sparse by {fmt(mean_delta_pct, 4)}%. That difference is effectively negligible, while no-sparse took {fmt(runtime_ratio_nosparse, 2)}x longer.</p>",
    "<h2>Distributed Stability</h2>",
    "<p>Distributed runs require a different stability check: consensus clusters. A good distributed result should preserve the main biology and have high within-cluster correlations. Weak clusters identify secondary programs that are less reproducible across subsets.</p>",
    table(["Run", "Consensus clusters", "Minimum cluster corr", "Mean of cluster minima", "Clusters below 0.90", "Mean consensus corr"], consensus_rows),
    "<p>R distributed sparse-on and Python distributed sparse-on both preserve the IFN program, but each has one weak consensus cluster around 0.78 minimum correlation. That means the dominant IFN biology is robust, while at least one secondary program is less stable under distributed subset matching.</p>",
    "<h2>Biological Suitability</h2>",
    "<p>The strongest biological criterion is whether the IFN-associated pattern has canonical IFN genes and stimulation-increased activity across expected immune cell types. This is preserved across the completed sparse and distributed runs despite label permutation.</p>",
    table(["Run", "IFN pattern", "Mean IFN stim-control delta", "Minimum cell-type delta", "Maximum cell-type delta"], ifn_rows),
    "<h2>Posterior Uncertainty For The IFN Program</h2>",
    "<p>Uncertainty summaries are available for the instrumented R runs. The IFN program has low top-gene loading CV in both sparse and no-sparse runs, supporting a stable leading IFN gene set.</p>",
    table(["Run", "IFN pattern", "Top-gene median loading CV", "Mean factor SD"], uncertainty_rows),
    "<h2>What Is Important And Why</h2>",
    "<ol>",
    "<li><strong>Trace plateau:</strong> shows the sampler reached a stable region. Sparse R and Python sparse traces are identical and reassuring.</li>",
    "<li><strong>Posterior meanChiSq:</strong> better summarizes the final posterior result than the final logged current-state trace. It does not support no-sparse as meaningfully better.</li>",
    "<li><strong>Top-gene stability and IFN activity:</strong> more important for this case study than tiny chi-square differences, because the narrative depends on biological program recovery.</li>",
    "<li><strong>Consensus stability:</strong> essential for distributed runs. Distributed sparse-on is suitable for runtime/main-biology claims, but not as the cleanest exact factorization anchor.</li>",
    "<li><strong>Runtime:</strong> no-sparse and local Python slowdowns are not justified by a corresponding fit or biology gain.</li>",
    "</ol>",
]

for title, rel, caption in figures:
    html_parts.extend([f"<h2>{html.escape(title)}</h2>", f"<p>{html.escape(caption)}</p>", f"<img src='{html.escape(rel)}' alt='{html.escape(title)}'>"])

html_parts.extend([
    "<h2>Recommended Interpretation For The Case Study</h2>",
    "<p>Use the sparse single-process model as the primary diagnostic and biological reference. Mention that no-sparse was tested and did not improve posterior-summary fit or biological interpretation enough to justify its runtime. Use R distributed sparse-on as an acceleration experiment that preserves the dominant IFN program but introduces consensus-level caveats for secondary patterns. Treat Python as a reproducibility check for biology, not as the performance implementation in the current local environment.</p>",
    "<h2>Generated Tables</h2>",
    "<ul>",
    "<li><code>tables/model_diagnostic_summary.csv</code></li>",
    "<li><code>tables/singleprocess_chisq_interpretive_summary.csv</code></li>",
    "<li><code>tables/combined_trace.csv</code></li>",
    "<li><code>tables/combined_consensus_summary.csv</code></li>",
    "<li><code>tables/ifn_activity_stim_minus_ctrl.csv</code></li>",
    "<li><code>tables/combined_uncertainty_summary.csv</code></li>",
    "</ul>",
    "</body></html>",
])

(OUT / "diagnostic_dashboard.html").write_text("\n".join(html_parts), encoding="utf-8")

md = [
    "# CoGAPS Diagnostics: Interpreted Report",
    "",
    "## Bottom Line",
    "",
    "The sparse single-process model is the best diagnostic anchor. R distributed sparse-on is the best acceleration candidate. R no-sparse is not meaningfully a better fit, and current local Python is biologically reproducible but not performance-competitive.",
    "",
    "## Decision Summary",
]
for item in decision_rows:
    md.append(f"- **{item['question']}** {item['answer']}. {item['why']}")
md.extend([
    "",
    "## Important Result-Specific Interpretations",
    "",
    f"- No-sparse has a lower-looking logged current-state trace, but posterior `meanChiSq` is slightly worse than sparse by {fmt(mean_delta_pct, 4)}%.",
    f"- No-sparse took {fmt(runtime_ratio_nosparse, 2)}x longer than R sparse instrumented.",
    "- Python single sparse light-local reproduces R sparse trace, IFN pattern, IFN correlation, top genes, total updates, and meanChiSq, but is slow in the current local Docker environment.",
    "- Distributed sparse-on preserves the dominant IFN program, but consensus diagnostics flag one weaker secondary cluster.",
    "- The IFN program is biologically stable across implementations: canonical IFN genes and positive stim-control activity deltas are preserved.",
    "",
    "## Updated Figures",
])
for title, rel, caption in figures:
    md.extend([f"### {title}", "", caption, "", f"![{title}]({rel})", ""])

(OUT / "diagnostic_interpretation.md").write_text("\n".join(md), encoding="utf-8")

print(OUT / "diagnostic_dashboard.html")
print(OUT / "diagnostic_interpretation.md")
