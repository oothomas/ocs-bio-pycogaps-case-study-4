#!/usr/bin/env python3
"""Assess stim-vs-ctrl directionality for CoGAPS pattern genes.

This script answers a reviewer-facing question that CoGAPS alone does not:
the CoGAPS A matrix identifies genes that define each pattern, but expression
direction needs to be estimated from the underlying count data.

The analysis uses replicate-aware pseudobulk summaries:

- extract top genes from the chosen CoGAPS result
- aggregate raw counts by replicate + condition, globally and within cell type
- library-size normalize pseudobulk counts to CPM
- compute paired log2(stim CPM + pseudocount) - log2(ctrl CPM + pseudocount)
- summarize direction and consistency across donors
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
PREPROCESSED_H5AD = ROOT / "data/processed/preprocessed_cells_hvg3000.h5ad"
CHOSEN_RESULT_H5AD = ROOT / "data/results/cogaps_K7_seed2_iter2000.h5ad"
OUT_DIR = ROOT / "data/results"

TOP_N_GENES = 15
MIN_CELLS_PER_PSEUDOBULK = 20
CPM_PSEUDOCOUNT = 0.5
STRONG_ABS_LOG2FC = 0.5
MODERATE_ABS_LOG2FC = 0.25
HIGH_CONSISTENCY = 0.75


@dataclass(frozen=True)
class DirectionSummary:
    direction: str
    evidence_strength: str


def pattern_columns(df: pd.DataFrame) -> list[str]:
    return sorted(
        [c for c in df.columns if str(c).startswith("Pattern")],
        key=lambda s: int(str(s).replace("Pattern", "")),
    )


def classify_direction(mean_log2fc: float, sign_consistency: float, n_pairs: int) -> DirectionSummary:
    if pd.isna(mean_log2fc) or n_pairs == 0:
        return DirectionSummary("not_tested", "insufficient_pairs")

    if mean_log2fc > 0:
        direction = "up_in_stim"
    elif mean_log2fc < 0:
        direction = "down_in_stim"
    else:
        direction = "mixed_or_unclear"

    abs_fc = abs(mean_log2fc)
    if abs_fc >= STRONG_ABS_LOG2FC and sign_consistency >= HIGH_CONSISTENCY:
        strength = "strong"
    elif abs_fc >= MODERATE_ABS_LOG2FC and sign_consistency >= 0.625:
        strength = "moderate"
    elif abs_fc < 0.1 or sign_consistency < 0.625:
        strength = "weak_or_mixed"
    else:
        strength = "suggestive"

    if direction == "mixed_or_unclear":
        strength = "weak_or_mixed"
    return DirectionSummary(direction, strength)


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjusted p-values, preserving missing values."""
    out = pd.Series(np.nan, index=pvalues.index, dtype=float)
    valid = pvalues.dropna().astype(float)
    if valid.empty:
        return out

    order = np.argsort(valid.to_numpy())
    ordered_p = valid.to_numpy()[order]
    n = len(ordered_p)
    adjusted = ordered_p * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    valid_index = valid.index.to_numpy()[order]
    out.loc[valid_index] = adjusted
    return out


def paired_wilcoxon_pvalue(values: np.ndarray) -> float:
    if len(values) < 2 or np.allclose(values, 0):
        return float("nan")
    try:
        return float(stats.wilcoxon(values, alternative="two-sided", zero_method="wilcox").pvalue)
    except ValueError:
        return float("nan")


def get_counts_matrix(adata: ad.AnnData):
    if "counts" in adata.layers:
        return adata.layers["counts"]
    return adata.X


def subset_counts(X, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    block = X[rows][:, cols]
    if sparse.issparse(block):
        return np.asarray(block.sum(axis=0)).ravel()
    return np.asarray(block).sum(axis=0)


def group_library_size(X, rows: np.ndarray) -> float:
    block = X[rows]
    return float(block.sum())


def extract_top_pattern_genes(result: ad.AnnData, top_n: int) -> pd.DataFrame:
    rows = []
    patterns = pattern_columns(result.obs)
    for pattern in patterns:
        weights = result.obs[pattern].sort_values(ascending=False).head(top_n)
        for rank, (gene, weight) in enumerate(weights.items(), start=1):
            rows.append(
                {
                    "pattern": pattern,
                    "gene": str(gene),
                    "rank": rank,
                    "cogaps_weight": float(weight),
                }
            )
    return pd.DataFrame(rows)


def make_pseudobulk(
    adata: ad.AnnData,
    genes: list[str],
    group_cols: list[str],
    min_cells: int,
) -> pd.DataFrame:
    missing = sorted(set(genes) - set(adata.var_names.astype(str)))
    if missing:
        raise ValueError(f"Pattern genes missing from processed AnnData: {missing}")

    X = get_counts_matrix(adata)
    gene_indices = np.asarray([adata.var_names.get_loc(gene) for gene in genes])
    obs = adata.obs.copy()
    obs["_row_index"] = np.arange(adata.n_obs)

    rows = []
    for group_values, group_df in obs.groupby(group_cols, observed=True):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        n_cells = len(group_df)
        if n_cells < min_cells:
            continue

        cell_indices = group_df["_row_index"].to_numpy()
        gene_counts = subset_counts(X, cell_indices, gene_indices)
        library_size = group_library_size(X, cell_indices)

        row = dict(zip(group_cols, group_values))
        row["n_cells"] = int(n_cells)
        row["library_size_total_counts"] = library_size
        row["library_size_pattern_genes"] = float(gene_counts.sum())
        for gene, count in zip(genes, gene_counts):
            row[gene] = float(count)
        rows.append(row)

    return pd.DataFrame(rows)


def add_cpm(pseudobulk: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    out = pseudobulk.copy()
    count_matrix = out[genes].to_numpy(dtype=float)
    library_sizes = out["library_size_total_counts"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        cpm = np.divide(
            count_matrix,
            library_sizes[:, None],
            out=np.zeros_like(count_matrix, dtype=float),
            where=library_sizes[:, None] > 0,
        ) * 1_000_000
    for i, gene in enumerate(genes):
        out[f"{gene}__cpm"] = cpm[:, i]
    return out


def paired_directionality(
    pseudobulk_cpm: pd.DataFrame,
    genes: list[str],
    pair_cols: list[str],
    context_cols: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return per-pair log2FC rows and one summary row per gene."""
    pair_rows = []
    summary_rows = []

    group_cols = [*pair_cols]
    grouped = pseudobulk_cpm.groupby(group_cols, observed=True, dropna=False)
    for pair_values, group in grouped:
        if not isinstance(pair_values, tuple):
            pair_values = (pair_values,)
        by_condition = group.set_index("condition", drop=False)
        if "ctrl" not in by_condition.index or "stim" not in by_condition.index:
            continue
        pair_info = dict(zip(pair_cols, pair_values))
        ctrl = by_condition.loc["ctrl"]
        stim = by_condition.loc["stim"]
        for gene in genes:
            ctrl_cpm = float(ctrl[f"{gene}__cpm"])
            stim_cpm = float(stim[f"{gene}__cpm"])
            log2fc = np.log2(stim_cpm + CPM_PSEUDOCOUNT) - np.log2(ctrl_cpm + CPM_PSEUDOCOUNT)
            pair_rows.append(
                {
                    **context_cols,
                    **pair_info,
                    "gene": gene,
                    "ctrl_cpm": ctrl_cpm,
                    "stim_cpm": stim_cpm,
                    "log2fc_stim_vs_ctrl": float(log2fc),
                    "ctrl_n_cells": int(ctrl["n_cells"]),
                    "stim_n_cells": int(stim["n_cells"]),
                }
            )

    pair_df = pd.DataFrame(pair_rows)
    if pair_df.empty:
        return pair_df, pd.DataFrame()

    summary_group_cols = ["gene", *[c for c in pair_cols if c != "replicate"]]
    for group_values, gene_df in pair_df.groupby(summary_group_cols, observed=True, dropna=False):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        summary_group = dict(zip(summary_group_cols, group_values))
        gene = summary_group["gene"]
        log2fc_values = gene_df["log2fc_stim_vs_ctrl"].to_numpy(dtype=float)
        mean_log2fc = float(np.mean(log2fc_values))
        median_log2fc = float(np.median(log2fc_values))
        n_pairs = int(len(log2fc_values))
        n_up = int(np.sum(log2fc_values > 0))
        n_down = int(np.sum(log2fc_values < 0))
        n_zero = int(np.sum(log2fc_values == 0))
        sign_consistency = max(n_up, n_down, n_zero) / n_pairs
        direction = classify_direction(mean_log2fc, sign_consistency, n_pairs)
        pvalue = paired_wilcoxon_pvalue(log2fc_values)
        summary_rows.append(
            {
                **context_cols,
                **summary_group,
                "mean_log2fc_stim_vs_ctrl": mean_log2fc,
                "median_log2fc_stim_vs_ctrl": median_log2fc,
                "n_pairs": n_pairs,
                "n_up_pairs": n_up,
                "n_down_pairs": n_down,
                "n_zero_pairs": n_zero,
                "sign_consistency": float(sign_consistency),
                "wilcoxon_pvalue": pvalue,
                "direction": direction.direction,
                "evidence_strength": direction.evidence_strength,
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    if not summary_df.empty:
        summary_df["wilcoxon_fdr"] = benjamini_hochberg(summary_df["wilcoxon_pvalue"])
    return pair_df, summary_df


def merge_pattern_annotations(summary: pd.DataFrame, top_genes: pd.DataFrame) -> pd.DataFrame:
    merged = top_genes.merge(summary, on="gene", how="left")
    sort_cols = [c for c in ["pattern", "rank", "cell_type"] if c in merged.columns]
    return merged.sort_values(sort_cols).reset_index(drop=True)


def make_pattern_heatmap(global_summary: pd.DataFrame, output_path: Path) -> None:
    ordered = global_summary.sort_values(["pattern", "rank"])
    labels = [f"{row.pattern}:{row.gene}" for row in ordered.itertuples()]
    values = ordered["mean_log2fc_stim_vs_ctrl"].to_numpy(dtype=float)[:, None]

    height = max(8, len(labels) * 0.18)
    fig, ax = plt.subplots(figsize=(5, height))
    vmax = max(1.0, float(np.nanmax(np.abs(values))))
    im = ax.imshow(values, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xticks([0])
    ax.set_xticklabels(["stim vs ctrl\nmean paired log2FC"])
    ax.set_title("Directionality of top CoGAPS pattern genes")
    fig.colorbar(im, ax=ax, fraction=0.08, pad=0.04, label="mean paired log2FC")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def make_pattern_direction_summary(global_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pattern, df in global_summary.groupby("pattern", observed=True):
        direction_counts = df["direction"].value_counts().to_dict()
        strength_counts = df["evidence_strength"].value_counts().to_dict()
        up_genes = (
            df.loc[df["direction"] == "up_in_stim"]
            .sort_values("mean_log2fc_stim_vs_ctrl", ascending=False)
            .head(5)["gene"]
        )
        down_genes = (
            df.loc[df["direction"] == "down_in_stim"]
            .sort_values("mean_log2fc_stim_vs_ctrl", ascending=True)
            .head(5)["gene"]
        )
        rows.append(
            {
                "pattern": pattern,
                "n_top_genes": int(len(df)),
                "n_up_in_stim": int(direction_counts.get("up_in_stim", 0)),
                "n_down_in_stim": int(direction_counts.get("down_in_stim", 0)),
                "n_mixed_or_unclear": int(direction_counts.get("mixed_or_unclear", 0)),
                "n_strong": int(strength_counts.get("strong", 0)),
                "n_moderate": int(strength_counts.get("moderate", 0)),
                "n_suggestive": int(strength_counts.get("suggestive", 0)),
                "n_weak_or_mixed": int(strength_counts.get("weak_or_mixed", 0)),
                "median_mean_log2fc": float(df["mean_log2fc_stim_vs_ctrl"].median()),
                "top_up_genes": ", ".join(up_genes),
                "top_down_genes": ", ".join(down_genes),
            }
        )
    return pd.DataFrame(rows).sort_values("pattern")


def write_readme(
    output_path: Path,
    n_genes: int,
    global_summary: pd.DataFrame,
    pattern_summary: pd.DataFrame,
) -> None:
    p2 = pattern_summary.loc[pattern_summary["pattern"] == "Pattern2"].iloc[0]
    p7 = pattern_summary.loc[pattern_summary["pattern"] == "Pattern7"].iloc[0]

    lines = [
        "# Pattern Gene Directionality Analysis",
        "",
        "This folder contains a replicate-aware directionality check for the genes that define the chosen CoGAPS model patterns.",
        "",
        "CoGAPS weights identify genes associated with each pattern, but they do not by themselves indicate whether those genes are up- or down-regulated. Direction was estimated separately from the underlying count data by comparing IFN-beta stimulated cells to control cells.",
        "",
        "## Inputs",
        "",
        f"- processed AnnData: `{PREPROCESSED_H5AD.relative_to(ROOT)}`",
        f"- chosen CoGAPS result: `{CHOSEN_RESULT_H5AD.relative_to(ROOT)}`",
        f"- top genes per pattern: `{TOP_N_GENES}`",
        f"- total unique genes tested: `{n_genes}`",
        "",
        "## Method",
        "",
        "- extracted top CoGAPS genes from the chosen model A matrix",
        "- aggregated raw counts from `.layers[\"counts\"]` into pseudobulk profiles",
        "- computed paired donor-level `log2(stim CPM + 0.5) - log2(ctrl CPM + 0.5)`",
        "- added paired Wilcoxon signed-rank p-values and Benjamini-Hochberg FDR values as targeted descriptive support",
        "- summarized direction globally across donors and within cell types where enough cells were available",
        "",
        "## Main Result",
        "",
        f"- `Pattern2`: {int(p2.n_up_in_stim)} of {int(p2.n_top_genes)} top genes are up in stimulation by the global paired pseudobulk summary.",
        f"- `Pattern7`: {int(p7.n_up_in_stim)} of {int(p7.n_top_genes)} top genes are up in stimulation by the global paired pseudobulk summary.",
        "- This supports treating Pattern2 as an IFN-stimulated program, with Pattern7 as a secondary stimulation-associated candidate.",
        "",
        "## Caveats",
        "",
        "- This is a targeted directionality check for pattern genes, not a genome-wide differential-expression analysis.",
        "- The by-cell-type summaries are descriptive and depend on cell counts per donor/condition/cell-type stratum.",
        "- The current dataset is an IFN-beta stimulation PBMC dataset; any dengue interpretation should be framed as biological relevance to IFN/inflammatory programs, not direct dengue infection evidence.",
        "",
        "## Output Files",
        "",
        "- `pattern_top_genes.csv`: top genes extracted from each CoGAPS pattern",
        "- `pattern_gene_directionality_global.csv`: global paired pseudobulk direction by gene",
        "- `pattern_gene_directionality_by_celltype.csv`: cell-type paired pseudobulk direction by gene",
        "- `pattern_direction_summary.csv`: one-row summary per pattern",
        "- `pattern_gene_directionality_pairs_global.csv`: donor-level global log2FC values",
        "- `pattern_gene_directionality_pairs_by_celltype.csv`: donor-level cell-type log2FC values",
        "- `pattern_gene_directionality_heatmap.png`: heatmap of global top-gene directionality",
        "",
    ]
    output_path.write_text("\n".join(lines))


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)

    adata = ad.read_h5ad(PREPROCESSED_H5AD)
    result = ad.read_h5ad(CHOSEN_RESULT_H5AD)

    top_genes = extract_top_pattern_genes(result, TOP_N_GENES)
    genes = sorted(top_genes["gene"].unique())

    global_pb = add_cpm(
        make_pseudobulk(
            adata=adata,
            genes=genes,
            group_cols=["replicate", "condition"],
            min_cells=MIN_CELLS_PER_PSEUDOBULK,
        ),
        genes,
    )
    global_pairs, global_gene_summary = paired_directionality(
        global_pb,
        genes,
        pair_cols=["replicate"],
        context_cols={"contrast_scope": "global"},
    )
    global_summary = merge_pattern_annotations(global_gene_summary, top_genes)

    by_cell_pb = add_cpm(
        make_pseudobulk(
            adata=adata,
            genes=genes,
            group_cols=["replicate", "condition", "cell_type"],
            min_cells=MIN_CELLS_PER_PSEUDOBULK,
        ),
        genes,
    )
    cell_pairs, cell_gene_summary = paired_directionality(
        by_cell_pb,
        genes,
        pair_cols=["replicate", "cell_type"],
        context_cols={"contrast_scope": "by_cell_type"},
    )
    cell_summary = (
        top_genes.merge(cell_gene_summary, on="gene", how="left")
        .sort_values(["pattern", "rank", "cell_type"])
        .reset_index(drop=True)
    )

    pattern_summary = make_pattern_direction_summary(global_summary)

    top_genes.to_csv(OUT_DIR / "pattern_top_genes.csv", index=False)
    global_summary.to_csv(OUT_DIR / "pattern_gene_directionality_global.csv", index=False)
    cell_summary.to_csv(OUT_DIR / "pattern_gene_directionality_by_celltype.csv", index=False)
    pattern_summary.to_csv(OUT_DIR / "pattern_direction_summary.csv", index=False)
    global_pairs.to_csv(OUT_DIR / "pattern_gene_directionality_pairs_global.csv", index=False)
    cell_pairs.to_csv(OUT_DIR / "pattern_gene_directionality_pairs_by_celltype.csv", index=False)
    global_pb.to_csv(OUT_DIR / "pseudobulk_counts_global.csv", index=False)
    by_cell_pb.to_csv(OUT_DIR / "pseudobulk_counts_by_celltype.csv", index=False)

    make_pattern_heatmap(global_summary, OUT_DIR / "pattern_gene_directionality_heatmap.png")
    write_readme(OUT_DIR / "README.md", len(genes), global_summary, pattern_summary)

    print(f"Wrote outputs to: {OUT_DIR}")
    print(pattern_summary.to_string(index=False))
    print("\nPattern2 global gene directions:")
    print(
        global_summary.loc[global_summary["pattern"] == "Pattern2", [
            "pattern",
            "rank",
            "gene",
            "mean_log2fc_stim_vs_ctrl",
            "n_up_pairs",
            "n_down_pairs",
            "wilcoxon_fdr",
            "direction",
            "evidence_strength",
        ]].to_string(index=False)
    )
    print("\nPattern7 global gene directions:")
    print(
        global_summary.loc[global_summary["pattern"] == "Pattern7", [
            "pattern",
            "rank",
            "gene",
            "mean_log2fc_stim_vs_ctrl",
            "n_up_pairs",
            "n_down_pairs",
            "wilcoxon_fdr",
            "direction",
            "evidence_strength",
        ]].to_string(index=False)
    )


if __name__ == "__main__":
    main()
