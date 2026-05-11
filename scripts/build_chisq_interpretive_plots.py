#!/usr/bin/env python3
"""Generate clearer single-process chi-square diagnostics."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


OUT = Path("data/diagnostics_model_comparison")
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.bbox"] = "tight"

trace = pd.read_csv(OUT / "tables" / "combined_trace.csv")
summary = pd.read_csv(OUT / "tables" / "model_diagnostic_summary.csv")

focus_runs = [
    "R sparse instrumented",
    "R no-sparse instrumented",
    "Python single sparse light-local",
]

colors = {
    "R sparse instrumented": "#4C72B0",
    "R no-sparse instrumented": "#55A868",
    "Python single sparse light-local": "#C44E52",
}
linestyles = {
    "R sparse instrumented": "-",
    "R no-sparse instrumented": "-",
    "Python single sparse light-local": "--",
}

df = trace[trace["run"].isin(focus_runs)].copy()
df["chisq_millions"] = df["chisq"] / 1_000_000
df["phase"] = np.where(df["relative_progress"] <= 0.5, "equilibration", "sampling")

means = (
    summary[summary["run"].isin(focus_runs)][["run", "meanChiSq"]]
    .dropna()
    .set_index("run")["meanChiSq"]
    .to_dict()
)
df["meanChiSq"] = df["run"].map(means)
df["chisq_minus_mean_millions"] = (df["chisq"] - df["meanChiSq"]) / 1_000_000
df["pct_above_run_min"] = (df["chisq"] - df.groupby("run")["chisq"].transform("min")) / df.groupby("run")["chisq"].transform("min") * 100


def plot_line(ax, data: pd.DataFrame, y: str, title: str, ylabel: str) -> None:
    for run in focus_runs:
        sub = data[data["run"] == run]
        ax.plot(
            sub["relative_progress"],
            sub[y],
            marker="o",
            markersize=4,
            linewidth=2.2,
            color=colors[run],
            linestyle=linestyles[run],
            label=run,
        )
    ax.axvline(0.5, color="black", linestyle="--", linewidth=1)
    ax.set_title(title)
    ax.set_xlabel("Relative progress (0.5 = end of equilibration)")
    ax.set_ylabel(ylabel)


fig, axes = plt.subplots(1, 2, figsize=(18, 7), gridspec_kw={"width_ratios": [1, 1.05]})
plot_line(
    axes[0],
    df,
    "chisq_millions",
    "Full trace",
    "Logged current-state chi-square (millions)",
)

sampling = df[df["phase"] == "sampling"].copy()
plot_line(
    axes[1],
    sampling,
    "chisq_millions",
    "Sampling-phase zoom",
    "Logged current-state chi-square (millions)",
)
for run, mean in means.items():
    axes[1].axhline(mean / 1_000_000, color=colors[run], linestyle=":", linewidth=1.8)
axes[1].text(
    0.51,
    axes[1].get_ylim()[0],
    "dotted lines = posterior meanChiSq",
    fontsize=10,
    va="bottom",
)
handles, labels = axes[1].get_legend_handles_labels()
axes[0].legend_.remove() if axes[0].legend_ else None
axes[1].legend(handles, labels, loc="best", fontsize=10)
fig.suptitle("Single-process chi-square: full trace and sampling zoom", y=1.02)
fig.savefig(FIG / "singleprocess_chisq_trace_interpretive.png")
plt.close(fig)


fig, ax = plt.subplots(figsize=(14, 7))
plot_line(
    ax,
    sampling,
    "chisq_minus_mean_millions",
    "Sampling trace relative to each run's posterior meanChiSq",
    "Logged chi-square minus meanChiSq (millions)",
)
ax.axhline(0, color="black", linewidth=1)
ax.legend(loc="best", fontsize=10)
fig.savefig(FIG / "singleprocess_chisq_minus_mean_sampling.png")
plt.close(fig)


fig, ax = plt.subplots(figsize=(14, 7))
plot_line(
    ax,
    df,
    "pct_above_run_min",
    "Trace shape after normalizing each run to its own minimum",
    "Percent above run-specific minimum chi-square",
)
ax.legend(loc="best", fontsize=10)
fig.savefig(FIG / "singleprocess_chisq_normalized_to_min.png")
plt.close(fig)


comparison = []
for run in focus_runs:
    sub = df[df["run"] == run]
    comparison.append(
        {
            "run": run,
            "final_logged_chisq": float(sub["chisq"].iloc[-1]),
            "minimum_logged_chisq": float(sub["chisq"].min()),
            "posterior_meanChiSq": float(means.get(run, np.nan)),
            "sampling_range": float(sampling[sampling["run"] == run]["chisq"].max() - sampling[sampling["run"] == run]["chisq"].min()),
            "runtime_min": float(summary.loc[summary["run"] == run, "runtime_min"].iloc[0]),
        }
    )
pd.DataFrame(comparison).to_csv(OUT / "tables" / "singleprocess_chisq_interpretive_summary.csv", index=False)

print(FIG / "singleprocess_chisq_trace_interpretive.png")
print(FIG / "singleprocess_chisq_minus_mean_sampling.png")
print(FIG / "singleprocess_chisq_normalized_to_min.png")
print(OUT / "tables" / "singleprocess_chisq_interpretive_summary.csv")
