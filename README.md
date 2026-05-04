# Bio-OCS CoGAPS PBMC Case Study Draft

This is a first-draft Bio-OCS/OTTR Quarto project for a case study on latent immune-response programs in PBMC single-cell RNA-seq using CoGAPS.

The project was created from the Bio-OCS latent-variable template and filled with draft content based on the local CoGAPS sweep, chosen model, and reviewer-triggered pattern-gene directionality analysis.

## Draft Status

This is not yet a publication-ready case study. It is a full first draft intended for review, revision, and rendering tests.

Known items to finalize:

- repository landing page and public documentation polish for `https://github.com/oothomas/ocs-bio-pycogaps-case-study-4`
- author list and citation text
- package/environment instructions beyond the current runtime image `othomas2/pycogaps-runtime-guide:0.2.0`
- whether saved data files should remain committed or be distributed separately
- final decision about mentioning dengue in the motivation
- optional enrichment analysis for Pattern2 and Pattern7
- final GitHub Actions policy for lightweight-versus-full rendering

## Main Learner Path

The default learner workflow is:

1. load `data/processed/preprocessed_cells_hvg3000.h5ad`
1. inspect AnnData metadata, condition balance, cell types, and sparsity
1. load or construct the CoGAPS input
1. load the saved chosen `K = 7` CoGAPS model in Python or R
1. extract gene weights and cell activities
1. interpret Pattern2 as the main IFN-associated program
1. inspect pattern-gene directionality
1. summarize Pattern2 across cell types

The full SLURM sweep is documented only as optional instructor/advanced material.

## Dual-language result layout

The repository now distinguishes Python and R saved results:

- `data/results_python/` for Python/PyCoGAPS saved artifacts
- `data/results_r/` for R/Bioconductor CoGAPS saved artifacts

The legacy `data/results/` directory is still present as a compatibility source for the existing Python full-result `.h5ad`, but the case study is moving toward explicit language-specific result folders.

## Lightweight Render Path

The repository now supports a lightweight render mode for GitHub or classroom contexts where the two large CoGAPS `.h5ad` files are not available:

- `data/processed/cogaps_input_genesxcells_hvg3000_float64.h5ad`
- `data/results/cogaps_K7_seed2_iter2000.h5ad`

When those large files are absent, the case study falls back to smaller exported pattern artifacts so the main interpretation sections can still render:

- `data/results_python/pattern_gene_weights.tsv.gz`
- `data/results_python/pattern_cell_activities.tsv.gz`
- `data/results_python/pattern_cell_activities_with_metadata.csv.gz`
- `data/results_python/pattern_top_genes.csv`
- `data/results_python/pattern_correlations.csv`
- `data/results_python/pattern_activity_by_celltype_condition.csv`
- `data/results_python/pattern_activity_by_replicate_condition.csv`
- `data/results_python/pattern_summary.csv`
- the matching `data/results_r/` files for the R implementation
- additional R directionality support files such as `pattern_top_genes_directionality.csv`,
  donor-pair directionality tables, pseudobulk count tables, and a small heatmap PNG

## Important Data Files

| Path | Purpose |
| --- | --- |
| `data/processed/preprocessed_cells_hvg3000.h5ad` | Processed cells x genes AnnData |
| `data/processed/cogaps_input_genesxcells_hvg3000_float64.h5ad` | Optional CoGAPS input |
| `data/results_python/cogaps_K7_seed2_iter2000.h5ad` | Saved chosen Python CoGAPS model |
| `data/results_python/cogaps_K7_seed2_iter2000.metrics.json` | Python chosen-model metadata |
| `data/results_r/cogaps_K7_seed2_iter2000.metrics.json` | R chosen-model metadata |
| `data/results_python/summary_by_K.csv` | Sweep summary used for rank justification |
| `data/results_python/pattern_gene_weights.tsv.gz` | Python lightweight gene-weight (`A`) matrix |
| `data/results_r/pattern_gene_weights.tsv.gz` | R lightweight gene-weight (`A`) matrix |
| `data/results_python/pattern_cell_activities.tsv.gz` | Python lightweight cell-activity (`P`) matrix |
| `data/results_r/pattern_cell_activities.tsv.gz` | R lightweight cell-activity (`P`) matrix |
| `data/results_python/pattern_top_genes.csv` | Python top genes per pattern |
| `data/results_r/pattern_top_genes.csv` | R top genes per pattern |
| `data/results_python/pattern_correlations.csv` | Python pattern correlations with stimulation |
| `data/results_r/pattern_correlations.csv` | R pattern correlations with stimulation |
| `data/results_python/pattern_activity_by_celltype_condition.csv` | Python pattern activity summaries by cell type and condition |
| `data/results_r/pattern_activity_by_celltype_condition.csv` | R pattern activity summaries by cell type and condition |
| `data/results_python/pattern_summary.csv` | Python compact per-pattern summary for lightweight rendering |
| `data/results_r/pattern_summary.csv` | R compact per-pattern summary for lightweight rendering |
| `data/results_python/pattern_gene_directionality_global.csv` | Python pattern-gene directionality across donors |
| `data/results_r/pattern_gene_directionality_global.csv` | R pattern-gene directionality across donors |
| `data/results_r/pattern_top_genes_directionality.csv` | R top-gene subset used for directionality summaries |
| `data/results_r/pattern_gene_directionality_pairs_global.csv` | R donor-paired global log2FC table |
| `data/results_r/pattern_gene_directionality_pairs_by_celltype.csv` | R donor-paired cell-type log2FC table |

## Rendering

From this folder, try:

```bash
quarto render
```

This draft now uses both Python and R chunks. The intended local environment is the Docker image `othomas2/pycogaps-runtime-guide:0.2.0`, which preinstalls:

- Python `PyCoGAPS`, `anndata`, `pandas`, `numpy`, `scipy`, `matplotlib`, and Jupyter
- R/Bioconductor `CoGAPS`, `zellkonverter`, `SingleCellExperiment`, and the supporting R packages used by the R alternative chunks

## Scientific Caveats

- The dataset tests IFN-beta stimulation in PBMCs, not direct dengue infection.
- Pattern2 is strongly supported as IFN-associated.
- Python and R CoGAPS outputs should be treated as parallel implementations, not as guaranteed numerically identical runs.
- Secondary patterns should remain cautiously interpreted.
- Directionality results are targeted to top pattern genes, not genome-wide differential expression.
- Cell-type differences are descriptive unless additional donor-aware modeling is added.
