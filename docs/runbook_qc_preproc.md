# qc_preproc module

Attaches metadata, writes pre-filter plots, filters/normalizes each study, and creates all-study pre-integration objects. Run from `a1_p1_qc/`.

## Commands

```bash
# dry-run full QC
snakemake --profile profiles/local -s workflow/modules/qc_preproc/Snakefile -n -p qc_all

# execute full QC
snakemake --profile profiles/local -s workflow/modules/qc_preproc/Snakefile qc_all

# one study through pre-integration
snakemake --profile profiles/local -s workflow/modules/qc_preproc/Snakefile \
  out/<study>/qc/.study_pre_integration_done
```

## Config

Primary file: `config/qc.yaml`.

Important knobs:

- `inputs.studies`: optional study allowlist; empty means infer from metadata.
- `qc.min_genes`, `qc.min_counts`, `qc.max_mito_frac`: filtering thresholds.
- `qc.backend`: `cpu`, `gpu`, or `auto` for RAPIDS-backed preprocessing.
- `qc.require_gpu`: fail instead of falling back when GPU is requested.
- `resources.*`: Snakemake threads/memory/GPU scheduling for heavy QC rules.

## Outputs

- `out/<study>/raw_metadata.h5ad`
- `out/<study>/qc/<study>_pre_integration.h5ad`
- `out/all_studies/all_raw.h5ad`
- `out/all_studies/all_pre_integration.h5ad`
