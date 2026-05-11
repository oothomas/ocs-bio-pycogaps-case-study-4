# Bio-OCS CoGAPS PBMC Case Study Draft

This is a first-draft Bio-OCS/OTTR Quarto project for a case study on latent immune-response programs in PBMC single-cell RNA-seq using CoGAPS.

The project was created from the Bio-OCS latent-variable template and filled with draft content based on the local CoGAPS sweep, chosen model, and reviewer-triggered pattern-gene directionality analysis.

## Draft Status

This is not yet a publication-ready case study. It is a full first draft intended for review, revision, and rendering tests.

Known items to finalize:

- repository landing page and public documentation polish for `https://github.com/oothomas/ocs-bio-pycogaps-case-study-4`
- author list and citation text
- package/environment instructions beyond the current OpenMP runtime image `othomas2/pycogaps-runtime-guide:0.3.0`
- final policy for distributing optional full model objects separately from lightweight GitHub artifacts
- final decision about mentioning dengue in the motivation
- optional enrichment analysis for the selected K6 IFN-associated pattern and secondary stimulation-associated pattern
- final GitHub Actions policy for lightweight-versus-full rendering

## Main Learner Path

The default learner workflow is:

1. load `data/processed/preprocessed_cells_hvg3000.h5ad`
1. inspect metadata, condition balance, cell types, and sparsity in R, with Python shown as a secondary implementation
1. load or construct the CoGAPS input
1. load selected R `K = 6` sparse-on CoGAPS artifacts
1. extract gene weights and cell activities
1. interpret the selected K6 IFN-associated pattern as the main IFN-associated program
1. inspect pattern-gene directionality
1. summarize the IFN-associated pattern across cell types

The full SLURM sweep is documented only as optional instructor/advanced material.

## Dual-language result layout

The repository now distinguishes R-primary and Python-secondary saved results:

- `data/results_r_k6_sparse_mt_t4_heavy/` for the selected R K6 sparse-on artifacts used by the main learner path
- `data/results_selected/` for selected K6 Python/PyCoGAPS-compatible artifacts used in the secondary tab
- older `data/results_python/`, `data/results_r/`, and `data/results/` folders for K7 reference or compatibility artifacts

The R and Python outputs should be treated as parallel implementations rather than guaranteed numerically identical objects. Pattern labels are run-specific; the selected R K6 IFN-associated pattern is `Pattern4`, while the selected Python K6 output labels the corresponding signal `Pattern5`.

## Lightweight Render Path

The repository supports a lightweight render mode for GitHub or classroom contexts where full CoGAPS model files are not available:

- optional R full model: `data/results_r_k6_sparse_mt_t4_heavy/cogaps_K6_seed2_iter2000.rds`
- optional Python full model: `data/results_selected/cogaps_K6_seed2_iter2000.h5ad`
- optional Python dense input: `data/processed/cogaps_input_genesxcells_hvg3000_float64.h5ad`

When those larger files are absent, the case study falls back to smaller exported pattern artifacts so the main interpretation sections can still render. The R-primary cooking-show artifacts are:

- `data/results_r_k6_sparse_mt_t4_heavy/cogaps_K6_seed2_iter2000.metrics.json`
- `data/results_r_k6_sparse_mt_t4_heavy/pattern_gene_weights.tsv.gz`
- `data/results_r_k6_sparse_mt_t4_heavy/pattern_cell_activities.tsv.gz`
- `data/results_r_k6_sparse_mt_t4_heavy/pattern_cell_activities_with_metadata.csv.gz`
- `data/results_r_k6_sparse_mt_t4_heavy/pattern_top_genes.csv`
- `data/results_r_k6_sparse_mt_t4_heavy/pattern_correlations.csv`
- `data/results_r_k6_sparse_mt_t4_heavy/pattern_activity_by_celltype_condition.csv`
- `data/results_r_k6_sparse_mt_t4_heavy/pattern_activity_by_replicate_condition.csv`
- `data/results_r_k6_sparse_mt_t4_heavy/pattern_summary.csv`
- `data/results_r_k6_sparse_mt_t4_heavy/pattern_gene_directionality_global.csv`
- `data/results_r_k6_sparse_mt_t4_heavy/pattern_gene_directionality_by_celltype.csv`
- `data/results_r_k6_sparse_mt_t4_heavy/pattern_direction_summary.csv`
- additional R directionality support files such as `pattern_top_genes_directionality.csv`,
  donor-pair directionality tables, pseudobulk count tables, and a small heatmap PNG

The secondary Python implementation keeps the same full-or-lightweight pattern using `data/results_selected/`.
Its full `.h5ad` result and dense genes x cells input are optional local artifacts, not GitHub-required files. When they are absent, the Python tabs load the smaller exported pattern matrices and summaries from `data/results_selected/`, just as the R tabs load the smaller R exports from `data/results_r_k6_sparse_mt_t4_heavy/`.

## Important Data Files

| Path | Purpose |
| --- | --- |
| `data/processed/preprocessed_cells_hvg3000.h5ad` | Processed cells x genes AnnData |
| `data/processed/cogaps_input_genesxcells_hvg3000_float64.h5ad` | Optional Python genes x cells CoGAPS input |
| `data/results_r_k6_sparse_mt_t4_heavy/cogaps_K6_seed2_iter2000.metrics.json` | Selected R K6 sparse-on model metadata |
| `data/model_selection/phase1b_lowerk_long_summary.csv` | Sweep summary used for rank justification |
| `data/results_r_k6_sparse_mt_t4_heavy/pattern_gene_weights.tsv.gz` | Selected R K6 lightweight gene-weight (`A`) matrix |
| `data/results_r_k6_sparse_mt_t4_heavy/pattern_cell_activities.tsv.gz` | Selected R K6 lightweight cell-activity (`P`) matrix |
| `data/results_r_k6_sparse_mt_t4_heavy/pattern_top_genes.csv` | Selected R K6 top genes per pattern |
| `data/results_r_k6_sparse_mt_t4_heavy/pattern_correlations.csv` | Selected R K6 pattern correlations with stimulation |
| `data/results_r_k6_sparse_mt_t4_heavy/pattern_activity_by_celltype_condition.csv` | Selected R K6 pattern activity summaries by cell type and condition |
| `data/results_r_k6_sparse_mt_t4_heavy/pattern_summary.csv` | Selected R K6 compact per-pattern summary for lightweight rendering |
| `data/results_r_k6_sparse_mt_t4_heavy/pattern_gene_directionality_global.csv` | Selected R K6 pattern-gene directionality across donors |
| `data/results_selected/pattern_gene_weights.tsv.gz` | Secondary Python K6 lightweight gene-weight matrix |
| `data/results_selected/pattern_cell_activities.tsv.gz` | Secondary Python K6 lightweight cell-activity matrix |
| `data/results_selected/pattern_top_genes.csv` | Secondary Python K6 top genes per pattern |
| `data/results_selected/pattern_correlations.csv` | Secondary Python K6 pattern correlations with stimulation |
| `data/results_selected/pattern_activity_by_celltype_condition.csv` | Secondary Python K6 pattern activity summaries by cell type and condition |
| `data/results_selected/pattern_direction_summary.csv` | Secondary Python K6 pattern-gene directionality summary |

## Rendering

From this folder, try:

```bash
quarto render
```

This draft uses R-first / Python-second tabbed chunks. The intended local validation environment is the OpenMP-enabled Docker image `othomas2/pycogaps-runtime-guide:0.3.0`, which preinstalls:

- R/Bioconductor `CoGAPS`, `zellkonverter`, `SingleCellExperiment`, and the supporting R packages used by the primary R chunks
- Python `PyCoGAPS`, `anndata`, `pandas`, `numpy`, `scipy`, `matplotlib`, and Jupyter for the secondary Python chunks

## Scientific Caveats

- The dataset tests IFN-beta stimulation in PBMCs, not direct dengue infection.
- The selected R K6 IFN-associated pattern is strongly supported by canonical interferon-stimulated genes and directionality.
- Python and R CoGAPS outputs should be treated as parallel implementations, not as guaranteed numerically identical runs.
- Secondary patterns should remain cautiously interpreted.
- Directionality results are targeted to top pattern genes, not genome-wide differential expression.
- Cell-type differences are descriptive unless additional donor-aware modeling is added.
