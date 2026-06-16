# integration_annotation module runbook

This project uses Snakemake. This runbook documents how to run the integration + annotation module.

Everything below assumes you run from the `a1_p1_qc/` directory.

## Quick start

Dry-run the full module with default/template config:

```bash
snakemake --profile profiles/local -s workflow/modules/integration_annotation/Snakefile -n -p integration_annotation_all
```

Dry-run the broad tier-1 profile:

```bash
snakemake --profile profiles/lca_tier_1_cpu_sysvi_scanvi -s workflow/modules/integration_annotation/Snakefile -n -p integration_annotation_all
```

Dry-run the myeloid tier-2 profile:

```bash
snakemake --profile profiles/lca_tier_2_gpu_sysvi_scanvi -s workflow/modules/integration_annotation/Snakefile -n -p integration_annotation_all
```

Dry-run the a1_p2 macrophage/DC reintegration profile:

```bash
snakemake --profile profiles/newint -s workflow/modules/integration_annotation/Snakefile -n -p integration_annotation_all
```

## Primary config

```text
config/integration_annotation.yaml
```

The default/template config is preconfigured to:

- build a Liver Cell Atlas **broad tier-1** reference atlas from
  - `rawData_human.zip`
  - `annot_humanAll.csv`
- run **SysVI** integration
- derive `technology` from `platform + chemistry`
- use **CPU** execution by default for local/dry-run-friendly configuration
- run **SCANVI** annotation against the configured reference atlas
- write `obs["tier_1_cell_type"]` and compatibility alias `obs["scanvi_pred"]`

Run-specific atlas/compute presets live under `profiles/`. If a profile sets
`profile_config_dir` and contains `integration_annotation.yaml`, that profile-local
YAML is used instead of `config/integration_annotation.yaml`.

Available label-transfer profiles include:

- `profiles/lca_tier_1_cpu_sysvi_scanvi` — broad LCA tier-1 labels, CPU SysVI/SCANVI.
- `profiles/lca_tier_2_gpu_sysvi_scanvi` — myeloid LCA tier-2 labels, GPU SysVI/SCANVI training.
- `profiles/newint` — reads `../a1_p2_def/out/trimint/mac.h5ad`, reruns SysVI/SCANVI against the myeloid LCA tier-2 reference, and writes `../a1_p2_def/out/trimint/all_annotated.h5ad` for a1_p2 downstream scripts.

## Targets

### Full all-studies integration + annotation + evaluation

```bash
snakemake --profile profiles/local -s workflow/modules/integration_annotation/Snakefile integration_annotation_all
```

Primary outputs:

- `out/integration_annotation/all_studies/all_integrated.h5ad`
- `out/integration_annotation/all_studies/all_annotated.h5ad`
- `out/integration_annotation/all_studies/evaluation_metrics.json`

With `profiles/newint`, the same target is redirected to a1_p2 paths:

- `../a1_p2_def/out/trimint/all_integrated.h5ad`
- `../a1_p2_def/out/trimint/all_annotated.h5ad`
- `../a1_p2_def/out/trimint/evaluation_metrics.json`

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
- SCANVI is configured without a batch key because the selected LCA reference metadata may expose `sample`, not `sample_id`.
- If you change a reference atlas profile, update `reference_atlas.output_h5ad`, `reference_atlas.report`, and `reference_atlas.ready_marker` together so tiers do not overwrite each other.
