#!/usr/bin/env python3
"""Build learner-facing figures for the selected R-primary K6 CoGAPS model."""

from __future__ import annotations

import json
from pathlib import Path
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = Path(
    os.environ.get(
        "COGAPS_SELECTED_RESULTS_DIR",
        ROOT / "data/results_r_k6_sparse_mt_t4_heavy",
    )
)
IMG_DIR = ROOT / "img"
METRICS_JSON = RESULTS_DIR / "cogaps_K6_seed2_iter2000.metrics.json"
CELL_ACTIVITIES = RESULTS_DIR / "pattern_cell_activities_with_metadata.csv.gz"
DIRECTIONALITY_HEATMAP = RESULTS_DIR / "pattern_gene_directionality_heatmap.png"


def pattern_columns(df: pd.DataFrame) -> list[str]:
    return sorted(
        [c for c in df.columns if str(c).startswith("Pattern")],
        key=lambda value: int(str(value).replace("Pattern", "")),
    )


def selected_ifn_pattern() -> str:
    with METRICS_JSON.open() as f:
        metrics = json.load(f)
    return str(metrics["ifn_pattern"])


def save_ifn_by_celltype_split(df: pd.DataFrame, ifn_pattern: str) -> None:
    cell_types = (
        df.groupby("cell_type", observed=True)[ifn_pattern]
        .median()
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    conditions = ["ctrl", "stim"]
    colors = {"ctrl": "#6B7280", "stim": "#2B8CBE"}
    positions = []
    labels = []
    data = []
    color_list = []
    for i, cell_type in enumerate(cell_types):
        base = i * 3
        for j, condition in enumerate(conditions):
            vals = df.loc[
                (df["cell_type"] == cell_type) & (df["condition"] == condition),
                ifn_pattern,
            ].dropna()
            positions.append(base + j)
            labels.append(cell_type if j == 0 else "")
            data.append(vals.to_numpy())
            color_list.append(colors[condition])

    fig, ax = plt.subplots(figsize=(12, 5.8))
    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.7,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.2},
    )
    for patch, color in zip(bp["boxes"], color_list):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
        patch.set_edgecolor("#111827")

    ax.set_xticks([i * 3 + 0.5 for i in range(len(cell_types))])
    ax.set_xticklabels(cell_types, rotation=35, ha="right")
    ax.set_ylabel(f"{ifn_pattern} activity")
    ax.set_title("Selected R K6 IFN-associated pattern activity by cell type and condition")
    handles = [
        plt.Line2D([0], [0], color=colors["ctrl"], lw=8, label="ctrl"),
        plt.Line2D([0], [0], color=colors["stim"], lw=8, label="stim"),
    ]
    ax.legend(handles=handles, title="condition", frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(IMG_DIR / "ifn_pattern_by_celltype_split_condition.png", dpi=180)
    plt.close(fig)


def save_ifn_by_celltype_all(df: pd.DataFrame, ifn_pattern: str) -> None:
    cell_types = (
        df.groupby("cell_type", observed=True)[ifn_pattern]
        .median()
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    data = [
        df.loc[df["cell_type"] == cell_type, ifn_pattern].dropna().to_numpy()
        for cell_type in cell_types
    ]
    fig, ax = plt.subplots(figsize=(11, 5.4))
    bp = ax.boxplot(
        data,
        widths=0.65,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.2},
    )
    for patch in bp["boxes"]:
        patch.set_facecolor("#2B8CBE")
        patch.set_alpha(0.72)
        patch.set_edgecolor("#111827")
    ax.set_xticklabels(cell_types, rotation=35, ha="right")
    ax.set_ylabel(f"{ifn_pattern} activity")
    ax.set_title("Selected R K6 IFN-associated pattern activity by cell type")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(IMG_DIR / "ifn_pattern_by_celltype_all.png", dpi=180)
    plt.close(fig)


def save_mean_pattern_by_condition(df: pd.DataFrame) -> None:
    patterns = pattern_columns(df)
    condition_order = ["ctrl", "stim"]
    means = (
        df.groupby("condition", observed=True)[patterns]
        .mean()
        .reindex(condition_order)
    )
    values = means.to_numpy()
    fig, ax = plt.subplots(figsize=(9, 3.8))
    im = ax.imshow(values, aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(len(patterns)))
    ax.set_xticklabels(patterns, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(condition_order)))
    ax.set_yticklabels(condition_order)
    ax.set_title("Selected R K6 mean pattern activity by condition")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            ax.text(j, i, f"{values[i, j]:.2f}", ha="center", va="center", color="white", fontsize=8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Mean activity")
    fig.tight_layout()
    fig.savefig(IMG_DIR / "mean_pattern_by_condition.png", dpi=180)
    plt.close(fig)


def main() -> None:
    IMG_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(CELL_ACTIVITIES)
    ifn_pattern = selected_ifn_pattern()
    save_ifn_by_celltype_split(df, ifn_pattern)
    save_ifn_by_celltype_all(df, ifn_pattern)
    save_mean_pattern_by_condition(df)
    if DIRECTIONALITY_HEATMAP.exists():
        target = IMG_DIR / "pattern_gene_directionality_heatmap.png"
        target.write_bytes(DIRECTIONALITY_HEATMAP.read_bytes())
    print(f"Selected K6 IFN-associated pattern: {ifn_pattern}")
    print(f"Wrote figures to {IMG_DIR}")


if __name__ == "__main__":
    main()
