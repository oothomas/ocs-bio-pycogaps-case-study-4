# Bio-OCS CoGAPS PBMC Case Study Draft

This is a first-draft Bio-OCS/OTTR Quarto project for a case study on latent immune-response programs in PBMC single-cell RNA-seq using CoGAPS.

The project was created from the Bio-OCS latent-variable template and filled with draft content based on the local CoGAPS sweep, chosen model, and reviewer-triggered pattern-gene directionality analysis.

## Draft Status

This is not yet a publication-ready case study. It is a full first draft intended for review, revision, and rendering tests.

Known items to finalize:

- final repository URL
- author list and citation text
- package/environment instructions
- whether saved data files should remain committed or be distributed separately
- final decision about mentioning dengue in the motivation
- optional enrichment analysis for Pattern2 and Pattern7
- final GitHub Actions policy for lightweight-versus-full rendering

## Main Learner Path

The default learner workflow is:

1. load `data/processed/preprocessed_cells_hvg3000.h5ad`
1. inspect AnnData metadata, condition balance, cell types, and sparsity
1. load or construct the CoGAPS genes x cells input
1. load the saved chosen `K = 7` CoGAPS model
1. extract gene weights and cell activities
1. interpret Pattern2 as the main IFN-associated program
1. inspect pattern-gene directionality
1. summarize Pattern2 across cell types

The full SLURM sweep is documented only as optional instructor/advanced material.

## Lightweight Render Path

The repository now supports a lightweight render mode for GitHub or classroom contexts where the two large CoGAPS `.h5ad` files are not available:

- `data/processed/cogaps_input_genesxcells_hvg3000_float64.h5ad`
- `data/results/cogaps_K7_seed2_iter2000.h5ad`

When those large files are absent, the case study falls back to smaller exported pattern artifacts so the main interpretation sections can still render:

- `data/results/pattern_gene_weights.tsv.gz`
- `data/results/pattern_cell_activities.tsv.gz`
- `data/results/pattern_cell_activities_with_metadata.csv.gz`
- `data/results/pattern_top_genes.csv`
- `data/results/pattern_correlations.csv`
- `data/results/pattern_activity_by_celltype_condition.csv`
- `data/results/pattern_activity_by_replicate_condition.csv`
- `data/results/pattern_summary.csv`

## Important Data Files

| Path | Purpose |
| --- | --- |
| `data/processed/preprocessed_cells_hvg3000.h5ad` | Processed cells x genes AnnData |
| `data/processed/cogaps_input_genesxcells_hvg3000_float64.h5ad` | Optional CoGAPS input |
| `data/results/cogaps_K7_seed2_iter2000.h5ad` | Saved chosen CoGAPS model |
| `data/results/cogaps_K7_seed2_iter2000.metrics.json` | Chosen model metadata |
| `data/results/summary_by_K.csv` | Sweep summary used for rank justification |
| `data/results/pattern_gene_weights.tsv.gz` | Lightweight gene-weight (`A`) matrix |
| `data/results/pattern_cell_activities.tsv.gz` | Lightweight cell-activity (`P`) matrix |
| `data/results/pattern_top_genes.csv` | Top genes per pattern |
| `data/results/pattern_correlations.csv` | Pattern correlations with stimulation |
| `data/results/pattern_activity_by_celltype_condition.csv` | Pattern activity summaries by cell type and condition |
| `data/results/pattern_summary.csv` | Compact per-pattern summary for lightweight rendering |
| `data/results/pattern_gene_directionality_global.csv` | Pattern-gene directionality across donors |
| `data/results/pattern_gene_directionality_by_celltype.csv` | Pattern-gene directionality by cell type |

## Rendering

From this folder, try:

```bash
quarto render
```

This draft uses Python chunks, so the Python environment also needs Jupyter available to Quarto. The core case study environment needs `anndata`, `pandas`, `numpy`, `scipy`, `matplotlib`, and `jupyter`. `PyCoGAPS` and `scanpy` are only needed for optional advanced scripts that rerun or re-aggregate models.

## Scientific Caveats

- The dataset tests IFN-beta stimulation in PBMCs, not direct dengue infection.
- Pattern2 is strongly supported as IFN-associated.
- Secondary patterns should remain cautiously interpreted.
- Directionality results are targeted to top pattern genes, not genome-wide differential expression.
- Cell-type differences are descriptive unless additional donor-aware modeling is added.
