# a1_p1_qc figures module

This module recreates the four ordered figures previously present as timestamped
`out/figures/annotate-*.png` Marimo exports.  It writes stable semantic PNG names
plus `out/figures/manifest.json` linking each output to the legacy reference.

Run from `a1_p1_qc/` with Snakemake target `figures_all`, or run the script:

```bash
python workflow/modules/figures/scripts/make_figures.py \
  --h5ad out/integration_annotation/all_studies/myeloid.h5ad \
  --metadata-csv config/liver_scRNAseq_metadata_everything.csv
```

Ordered figures:

1. `01_umap_study_accession.png` — cell-level UMAP from `.obsm['X_umap']`, colored by `.obs['study_accession']`.
2. `02_unique_samples_by_fibrosis_stage.png` — unique sample counts by fibrosis stage from metadata summaries or metadata CSV.
3. `03_unique_samples_by_study_accession.png` — unique sample counts by study accession.
4. `04_cell_type_study_composition.png` — cell-type-by-study composition from the annotated H5AD obs table.

The H5AD object is downstream of integration/annotation; metadata summaries are
cohort diagnostics derived from the configured metadata table.
