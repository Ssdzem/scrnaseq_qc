# qc_preproc module

Attaches metadata, writes pre-filter plots, normalizes each study, and creates all-study pre-integration objects. Run from the repository root.

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

Primary file: `config/qc.yaml`. Cells are filtered earlier, in standardize, through
each study's `input.qc` block; this module does not drop cells.

Important knobs:

- `inputs.studies`: optional study allowlist; empty means infer from metadata.
- `qc.hvg_n_top_genes`, `qc.hvg_batch_key`: highly variable gene selection.
- `qc.backend`: `cpu`, `gpu`, or `auto` for RAPIDS-backed preprocessing.
- `qc.require_gpu`: fail instead of falling back when GPU is requested.
- `resources.*`: Snakemake threads/memory/GPU scheduling for heavy QC rules.

## Outputs

- `out/<study>/raw_metadata.h5ad`
- `out/<study>/qc/<study>_pre_integration.h5ad`
- `out/all_studies/all_raw.h5ad`
- `out/all_studies/all_pre_integration.h5ad`
