# integration_annotation module runbook

This project uses Snakemake. This runbook documents how to run the integration + annotation module.

Everything below assumes you run from the `a1_p1_qc/` directory.

## Quick start

Dry-run the full module:

```bash
snakemake --profile profiles/local -s workflow/modules/integration_annotation/Snakefile -n -p integration_annotation_all
```

## Primary config

```text
config/integration_annotation.yaml
```

This worktree is preconfigured to:

- build a Liver Cell Atlas **myeloid** reference atlas from
  - `rawData_human.zip`
  - `annot_humanMyeloid.csv`
- run **SysVI** integration
- derive `technology` from `platform + chemistry`
- use **GPU** execution for integration + SCANVI label transfer
- run **SCANVI** annotation against the configured reference atlas

## Targets

### Full all-studies integration + annotation + evaluation

```bash
snakemake --profile profiles/local -s workflow/modules/integration_annotation/Snakefile integration_annotation_all
```

Primary outputs:

- `out/integration_annotation/all_studies/all_integrated.h5ad`
- `out/integration_annotation/all_studies/all_annotated.h5ad`
- `out/integration_annotation/all_studies/evaluation_metrics.json`

### Per-study integration + annotation + evaluation

```bash
snakemake --profile profiles/local -s workflow/modules/integration_annotation/Snakefile \
  out/integration_annotation/<study>/.<study>_integration_annotation_done
```

Per-study outputs include:

- `out/integration_annotation/<study>/<study>_integrated.h5ad`
- `out/integration_annotation/<study>/<study>_annotated.h5ad`
- `out/integration_annotation/<study>/<study>_evaluation_metrics.json`

### Reference atlas only

```bash
snakemake --profile profiles/local -s workflow/modules/integration_annotation/Snakefile \
  data/ref/.reference_myeloid_ready
```

Reference outputs:

- `data/ref/lca_human_myeloid_reference.h5ad`
- `data/ref/reference_build_report.myeloid.json`

## Notes

- The module consumes QC outputs such as `out/all_studies/all_pre_integration.h5ad`.
- SCANVI is configured without a batch key because the selected myeloid reference metadata exposes `sample`, not `sample_id`.
- If you change the reference atlas again, update `reference_atlas.output_h5ad`, `reference_atlas.report`, and `reference_atlas.ready_marker` together.
