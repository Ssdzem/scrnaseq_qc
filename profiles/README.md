# Snakemake profiles

Profiles are run presets: cores, conda, rerun policy, and optional module-local config overrides.

## Config lookup

Module config resolution order:

1. explicit `--config <module>_config=<path>`;
2. `profiles/<profile>/<module>.yaml` when `profile_config_dir` is set;
3. `config/<module>.yaml`.

Standardization specs can also be overridden with `profiles/<profile>/standardize/*.yaml` or `--config standardize_config_dir=<dir>`.

Only templates live here, and git ignores everything else in this folder. A real run's
profile holds its studies, paths and thresholds, so keep it in the project that uses this
repository as a submodule, starting from a copy of a template:

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

## Templates

- `local/`: safe local defaults; falls back to `config/*.yaml`.
- `lca_tier_1_cpu_sysvi_scanvi/`: 16-core CPU tier-1 labels; writes `tier_1_cell_type`.
- `lca_tier_2_gpu_sysvi_scanvi/`: 16-core + 1-GPU tier-2 labels; writes `tier_2_cell_type`.

## Example

```bash
snakemake --profile profiles/lca_tier_2_gpu_sysvi_scanvi \
  -s workflow/modules/integration_annotation/Snakefile \
  -n -p integration_annotation_all
```
