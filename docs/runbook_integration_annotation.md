# integration_annotation module

Runs reference build, integration, label transfer, and evaluation. Run from the repository root.

## Commands

```bash
# default CPU tier-1 dry-run
snakemake --profile profiles/local -s workflow/modules/integration_annotation/Snakefile -n -p integration_annotation_all

# tier-1 CPU profile
snakemake --profile profiles/lca_tier_1_cpu_sysvi_scanvi -s workflow/modules/integration_annotation/Snakefile integration_annotation_all

# tier-2 GPU profile
snakemake --profile profiles/lca_tier_2_gpu_sysvi_scanvi -s workflow/modules/integration_annotation/Snakefile integration_annotation_all
```

## Config

Default file: `config/integration_annotation.yaml`.
Profile override: `profiles/<profile>/integration_annotation.yaml`.

Important knobs:

- `reference_atlas.*`: LCA download/build target.
- `integration.method`: `sysvi`, `scvi`, or `harmony`.
- `integration.backend`: `cpu`, `gpu`, or `auto` for graph/preprocessing backend.
- `integration.use_gpu`, `integration.require_gpu`, `integration.accelerator`, `integration.devices`: training placement.
- `annotation.method`: normally `scanvi`.
- `evaluation.scib.n_cores`: scIB parallelism.
- `resources.*`: Snakemake threads/memory/GPU scheduling.

## Main outputs

Default all-study outputs:

```text
out/integration_annotation/all_studies/all_integrated.h5ad
out/integration_annotation/all_studies/all_annotated.h5ad
out/integration_annotation/all_studies/evaluation_metrics.json
```

