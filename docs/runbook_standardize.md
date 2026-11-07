# standardize module

Converts study-specific raw inputs to canonical `raw.h5ad` files. Run from the repository root.

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

Each study has one spec. Copy `config/templates/standardize_study.yaml` to
`<your profile>/standardize/<study>.yaml`, or pass `--config standardize_config_dir=<dir>`.

Important per-study knobs:

- `pipeline`: cleaner sequence.
- `input.*`: raw files to unpack/standardize.
- `input.qc`: per-sample MAD filtering (`k`, `mito_upper_k`), miQC and doublet scoring.
  Optional `min_genes` and `max_pct_counts_mt` add absolute limits on top of MAD.
  Leave the block out to keep every cell.
- `input.workers`: Snakemake threads and script workers.
- `input.mem_mb`: optional Snakemake memory resource.

## Outputs

- `out/<study>/raw.h5ad`
- `out/<study>/trimming_tables.csv`
- `out/<study>/qc_raw/raw_qc_metrics.csv` — per-cell total_counts, n_genes_by_counts,
  pct_counts_mt (computed from `raw.h5ad` as-is, independent of whether the
  spec's `input.qc` block filtered anything).
- `out/<study>/qc_raw/pictures/raw_{log1p_total_counts,log1p_n_genes_by_counts,pct_counts_mt}.png`
  — per-sample distribution histograms.
- `out/<study>/qc_raw/pictures/raw_counts_vs_genes_by_mito.png` — counts vs.
  genes scatter, colored by %mito.

Use these to pick cutoffs before deciding whether/what to set in a spec's
`input.qc` block (or a study with none, like a raw-data-review pass, just
gets these plots with no cells dropped).
