# standardize module

Converts study-specific raw inputs to canonical `raw.h5ad` files. Run from `a1_p1_qc/`.

## Commands

```bash
# dry-run all configured studies
snakemake --profile profiles/local -s workflow/modules/standardize/Snakefile -n -p standardize_all

# execute all configured studies
snakemake --profile profiles/local -s workflow/modules/standardize/Snakefile standardize_all

# execute one study
snakemake --profile profiles/local -s workflow/modules/standardize/Snakefile out/<study>/raw.h5ad
```

## Config

Study specs live in `config/standardize/*.yaml`.

Profile override: put specs in `profiles/<profile>/standardize/*.yaml` or pass `--config standardize_config_dir=<dir>`.

Important per-study knobs:

- `pipeline`: cleaner sequence.
- `input.*`: raw files to unpack/standardize.
- `input.workers`: Snakemake threads and script workers.
- `input.mem_mb`: optional Snakemake memory resource.

## Outputs

- `out/<study>/raw.h5ad`
- `out/<study>/trimming_tables.csv`
