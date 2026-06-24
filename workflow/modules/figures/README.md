# a1_p1_qc figures module

This module generates four ordered figures from canonical pipeline objects.  It
writes stable semantic filenames plus `out/figures/manifest.json` recording the
source object, structure, and interpretation for each image.

Each figure is written as a **publication-ready vector `.pdf`** (embedded
editable TrueType fonts) plus a **`.png` preview** (the Snakemake-tracked
output). Styling uses [SciencePlots](https://github.com/garrettj403/SciencePlots)
(`science` + `nature`, `no-latex`) when installed, falling back to a built-in
300-dpi rcParams style so figures always render.

Run from `a1_p1_qc/` with Snakemake target `figures_all`, or run the script:

```bash
python workflow/modules/figures/scripts/make_figures.py \
  --h5ad out/integration_annotation/all_studies/all_annotated.h5ad \
  --metadata-csv config/liver_scRNAseq_metadata_everything.csv
```

Ordered figures:

1. `01_umap_study_accession.png` — cell-level UMAP from `.obsm['X_umap']`, colored by `.obs['study_accession']`.
2. `02_unique_samples_by_fibrosis_stage.png` — included (`include` truthy) unique sample counts by fibrosis stage from the declared metadata CSV.
3. `03_unique_samples_by_study_accession.png` — included (`include` truthy) unique sample counts by study accession.
4. `04_cell_type_study_composition.png` — canonical `.obs['tier_1_cell_type']` by `.obs['study_accession']` composition from the annotated H5AD obs table.

The figure contract requires the canonical `study_accession` and
`tier_1_cell_type` obs columns. The default H5AD object is the canonical
all-study integration/annotation handoff
(`out/integration_annotation/all_studies/all_annotated.h5ad`). Override
`figures.h5ad` in Snakemake config or pass `--h5ad` directly only when a
different annotated object is intentionally being reviewed. Metadata cohort
diagnostics are derived from the declared metadata CSV input after applying the
same include-truthy cohort filter used by upstream standardization.
