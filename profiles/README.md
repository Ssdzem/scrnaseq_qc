# Snakemake profiles and profile-local module configs

Profiles are reproducible run presets. The default/template YAMLs stay under
`config/`. A profile can override any module config by placing a YAML with the
same filename under that profile directory.

Resolution order for module configs is:

1. explicit Snakemake config key, e.g. `--config integration_annotation_config=...`
2. `profiles/<profile>/MODULE.yaml` when the profile sets `profile_config_dir`
3. `config/MODULE.yaml`

For standardization specs, a profile can override the whole study-spec directory
with `profiles/<profile>/standardize/*.yaml`; otherwise `config/standardize/*.yaml`
is used.

Each profile `config.yaml` should point back to its directory:

```yaml
cores: 1
printshellcmds: true
rerun-incomplete: true
keep-going: true
use-conda: true
config:
  - profile_config_dir=profiles/<profile-name>
```

Example dry-runs from `a1_p1_qc/`:

```bash
snakemake --profile profiles/lca_tier_1_cpu_sysvi_scanvi \
  -s workflow/modules/integration_annotation/Snakefile \
  -n -p integration_annotation_all

snakemake --profile profiles/lca_tier_2_gpu_sysvi_scanvi \
  -s workflow/modules/integration_annotation/Snakefile \
  -n -p integration_annotation_all
```

## Current label-transfer profiles

- `lca_tier_1_cpu_sysvi_scanvi/`
  - broad Liver Cell Atlas labels from `annot_humanAll.csv`
  - CPU SysVI/SCANVI settings
  - writes `obs["tier_1_cell_type"]`
  - also writes compatibility alias `obs["scanvi_pred"]`
- `lca_tier_2_gpu_sysvi_scanvi/`
  - myeloid Liver Cell Atlas labels from `annot_humanMyeloid.csv`
  - GPU SysVI/SCANVI training settings
  - writes `obs["tier_2_cell_type"]`
  - also writes compatibility alias `obs["scanvi_predicted_cell_type"]`
