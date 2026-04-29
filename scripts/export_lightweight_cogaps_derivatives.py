from __future__ import annotations

from pathlib import Path
import os

import anndata as ad
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PREPROCESSED_H5AD = ROOT / "data/processed/preprocessed_cells_hvg3000.h5ad"
CHOSEN_RESULT_H5AD = ROOT / "data/results/cogaps_K7_seed2_iter2000.h5ad"
OUT_DIR = ROOT / "data/results"

PATTERN_GENE_WEIGHTS_TSV_GZ = OUT_DIR / "pattern_gene_weights.tsv.gz"
PATTERN_CELL_ACTIVITIES_TSV_GZ = OUT_DIR / "pattern_cell_activities.tsv.gz"
PATTERN_CELL_ACTIVITIES_WITH_METADATA_CSV_GZ = OUT_DIR / "pattern_cell_activities_with_metadata.csv.gz"
PATTERN_TOP_GENES_CSV = OUT_DIR / "pattern_top_genes.csv"
PATTERN_CORRELATIONS_CSV = OUT_DIR / "pattern_correlations.csv"
PATTERN_ACTIVITY_BY_CELLTYPE_CONDITION_CSV = OUT_DIR / "pattern_activity_by_celltype_condition.csv"
PATTERN_ACTIVITY_BY_REPLICATE_CONDITION_CSV = OUT_DIR / "pattern_activity_by_replicate_condition.csv"
PATTERN_SUMMARY_CSV = OUT_DIR / "pattern_summary.csv"

STIM_LABEL = "stim"


def pattern_columns(df: pd.DataFrame) -> list[str]:
    return sorted(
        [c for c in df.columns if str(c).startswith("Pattern")],
        key=lambda value: int(str(value).replace("Pattern", "")),
    )


def top_genes_from_matrix(matrix: pd.DataFrame, pattern: str, n: int = 50) -> pd.DataFrame:
    top = matrix[pattern].sort_values(ascending=False).head(n)
    top_df = top.rename("weight").reset_index()
    top_df.columns = ["gene", "weight"]
    top_df.insert(0, "rank", range(1, len(top_df) + 1))
    top_df.insert(0, "pattern", pattern)
    return top_df


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    adata_cells = ad.read_h5ad(PREPROCESSED_H5AD)
    result = ad.read_h5ad(CHOSEN_RESULT_H5AD)

    pattern_names = pattern_columns(result.obs)

    # In the saved CoGAPS result, genes live in obs and cells live in var.
    A = result.obs[pattern_names].copy()
    A.index.name = "gene"

    P = result.var[pattern_names].copy()
    P.index = result.var_names.astype(str)
    P.index.name = "cell_barcode"

    cell_meta_cols = [c for c in ["condition", "cell_type", "replicate", "label"] if c in adata_cells.obs.columns]
    cell_meta = adata_cells.obs.loc[P.index, cell_meta_cols].copy()
    cell_meta.index.name = "cell_barcode"

    if "X_umap" in adata_cells.obsm:
        umap = pd.DataFrame(
            adata_cells.obsm["X_umap"],
            index=adata_cells.obs_names.astype(str),
            columns=["umap_1", "umap_2"],
        )
        cell_meta = cell_meta.join(umap.loc[P.index], how="left")

    if "X_pca" in adata_cells.obsm:
        pca = pd.DataFrame(
            adata_cells.obsm["X_pca"][:, :2],
            index=adata_cells.obs_names.astype(str),
            columns=["pca_1", "pca_2"],
        )
        cell_meta = cell_meta.join(pca.loc[P.index], how="left")

    cell_activities_with_metadata = cell_meta.join(P, how="left").reset_index()

    stim_indicator = (cell_meta["condition"].astype(str) == STIM_LABEL).astype(int)
    pattern_corrs = (
        pd.Series(
            {pattern: P[pattern].corr(stim_indicator) for pattern in pattern_names},
            name="corr_with_stim",
        )
        .sort_values(ascending=False)
        .reset_index()
        .rename(columns={"index": "pattern"})
    )

    top_genes = pd.concat(
        [top_genes_from_matrix(A, pattern, n=50) for pattern in pattern_names],
        ignore_index=True,
    )

    activity_by_celltype_condition = pd.concat(
        [
            cell_activities_with_metadata
            .groupby(["cell_type", "condition"], observed=True)[pattern]
            .agg(["mean", "median", "count"])
            .reset_index()
            .rename(columns={"mean": "mean_activity", "median": "median_activity", "count": "n_cells"})
            .assign(pattern=pattern)
            for pattern in pattern_names
        ],
        ignore_index=True,
    )[
        ["pattern", "cell_type", "condition", "mean_activity", "median_activity", "n_cells"]
    ]

    activity_by_replicate_condition = pd.concat(
        [
            cell_activities_with_metadata
            .groupby(["replicate", "condition"], observed=True)[pattern]
            .agg(["mean", "median", "count"])
            .reset_index()
            .rename(columns={"mean": "mean_activity", "median": "median_activity", "count": "n_cells"})
            .assign(pattern=pattern)
            for pattern in pattern_names
        ],
        ignore_index=True,
    )[
        ["pattern", "replicate", "condition", "mean_activity", "median_activity", "n_cells"]
    ]

    pattern_summary_rows: list[dict[str, object]] = []
    corr_lookup = pattern_corrs.set_index("pattern")["corr_with_stim"]
    for pattern in pattern_names:
        top8 = ", ".join(top_genes.query("pattern == @pattern").head(8)["gene"].tolist())
        top15 = ", ".join(top_genes.query("pattern == @pattern").head(15)["gene"].tolist())
        dominant = (
            cell_activities_with_metadata
            .groupby("cell_type", observed=True)[pattern]
            .mean()
            .sort_values(ascending=False)
        )
        condition_means = (
            cell_activities_with_metadata
            .groupby("condition", observed=True)[pattern]
            .mean()
            .to_dict()
        )
        pattern_summary_rows.append(
            {
                "pattern": pattern,
                "corr_with_stim": float(corr_lookup[pattern]),
                "dominant_cell_type": dominant.index[0],
                "dominant_cell_type_mean_activity": float(dominant.iloc[0]),
                "mean_activity_ctrl": float(condition_means.get("ctrl", float("nan"))),
                "mean_activity_stim": float(condition_means.get("stim", float("nan"))),
                "top_genes_top8": top8,
                "top_genes_top15": top15,
            }
        )

    pattern_summary = pd.DataFrame(pattern_summary_rows).sort_values("corr_with_stim", ascending=False)

    A.reset_index().to_csv(PATTERN_GENE_WEIGHTS_TSV_GZ, sep="\t", index=False, compression="gzip")
    P.reset_index().to_csv(PATTERN_CELL_ACTIVITIES_TSV_GZ, sep="\t", index=False, compression="gzip")
    cell_activities_with_metadata.to_csv(PATTERN_CELL_ACTIVITIES_WITH_METADATA_CSV_GZ, index=False, compression="gzip")
    top_genes.to_csv(PATTERN_TOP_GENES_CSV, index=False)
    pattern_corrs.to_csv(PATTERN_CORRELATIONS_CSV, index=False)
    activity_by_celltype_condition.to_csv(PATTERN_ACTIVITY_BY_CELLTYPE_CONDITION_CSV, index=False)
    activity_by_replicate_condition.to_csv(PATTERN_ACTIVITY_BY_REPLICATE_CONDITION_CSV, index=False)
    pattern_summary.to_csv(PATTERN_SUMMARY_CSV, index=False)

    for path in [
        PATTERN_GENE_WEIGHTS_TSV_GZ,
        PATTERN_CELL_ACTIVITIES_TSV_GZ,
        PATTERN_CELL_ACTIVITIES_WITH_METADATA_CSV_GZ,
        PATTERN_TOP_GENES_CSV,
        PATTERN_CORRELATIONS_CSV,
        PATTERN_ACTIVITY_BY_CELLTYPE_CONDITION_CSV,
        PATTERN_ACTIVITY_BY_REPLICATE_CONDITION_CSV,
        PATTERN_SUMMARY_CSV,
    ]:
        size_mb = os.path.getsize(path) / (1024 * 1024)
        print(f"{path.name}: {size_mb:.2f} MB")


if __name__ == "__main__":
    main()
