# Snakemake profiles

Profiles are run presets: cores, conda, rerun policy, and optional module-local config overrides.

## Config lookup

Module config resolution order:

1. explicit `--config <module>_config=<path>`;
2. `profiles/<profile>/<module>.yaml` when `profile_config_dir` is set;
3. `config/<module>.yaml`.

Standardization specs can also be overridden with `profiles/<profile>/standardize/*.yaml` or `--config standardize_config_dir=<dir>`.

A profile does not have to live in this repository. Keeping one in the project that uses
this repository as a submodule leaves this code untouched:

```bash
snakemake --workflow-profile /path/to/project/profile -s workflow/modules/standardize/Snakefile -n
```

with `profile_config_dir=/path/to/project/profile` in its `config:` list.

## Required profile shape

```yaml
cores: <max threads this profile should allow>
printshellcmds: true
rerun-incomplete: true
keep-going: true
use-conda: true
config:
  - profile_config_dir=profiles/<profile-name>
```

## Current profiles

- `local/`: safe local defaults; falls back to `config/*.yaml`.
- `lca_tier_1_cpu_sysvi_scanvi/`: 16-core CPU tier-1 labels; writes `tier_1_cell_type`.
- `lca_tier_2_gpu_sysvi_scanvi/`: 16-core + 1-GPU tier-2 labels; writes `tier_2_cell_type`.
- `validation/`, `validation_tier2/`: held-out liver cohorts on a machine with no GPU (Harmony and CPU scANVI).
- `jruiz_first_qc_run/`: standardize-only pass for the cervical cancer cohorts.
- `newint/`: 40-core + 1-GPU a1_p2 macrophage/DC reintegration; reads `../a1_p2_def/out/trimint/mac.h5ad` and writes back to `../a1_p2_def/out/trimint/`.

## Example

```bash
snakemake --profile profiles/lca_tier_2_gpu_sysvi_scanvi \
  -s workflow/modules/integration_annotation/Snakefile \
  -n -p integration_annotation_all
```
