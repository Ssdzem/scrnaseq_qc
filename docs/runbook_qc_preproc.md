# qc_preproc module runbook

This project uses Snakemake. This runbook mirrors the style of `docs/runbook_fetch.md` and documents how to run the QC preprocessing module.

Everything below assumes you run from the `a1_p1_qc/` directory.

## Quick start

Dry-run the full QC module:

```bash
snakemake --profile profiles/local -s workflow/modules/qc_preproc/Snakefile -n -p qc_all
```

## Snakemake profile

Use the local profile directly:

```bash
snakemake --profile profiles/local -s workflow/modules/qc_preproc/Snakefile <target>
```

## Configuration inputs

Primary config:

```text
config/qc.yaml
```

Study selection behavior:

- If `inputs.studies` is non-empty, QC runs only those studies.
- Otherwise, studies are inferred from unique `study_accession` values in `inputs.metadata_csv`.

## Targets

### Run the complete QC module

Target:

```text
qc_all
```

Command:

```bash
snakemake --profile profiles/local -s workflow/modules/qc_preproc/Snakefile qc_all
```

Done marker:

```text
out/all_studies/.qc_all_done
```

### Run per-study QC stages for one study

#### 1) Pre-trim metrics and figures

Target:

```text
out/<study>/qc/.pretrim_pictures_done
```

Command:

```bash
snakemake --profile profiles/local -s workflow/modules/qc_preproc/Snakefile \
  out/<study>/qc/.pretrim_pictures_done
```

#### 2) Pre-integration filtered object

Target:

```text
out/<study>/qc/.study_pre_integration_done
```

Command:

```bash
snakemake --profile profiles/local -s workflow/modules/qc_preproc/Snakefile \
  out/<study>/qc/.study_pre_integration_done
```

#### 3) Study explore object

Target:

```text
out/<study>/qc/.study_explore_done
```

Command:

```bash
snakemake --profile profiles/local -s workflow/modules/qc_preproc/Snakefile \
  out/<study>/qc/.study_explore_done
```

### Run all-studies concatenation/integration stages explicitly

Targets:

```text
out/all_studies/.all_raw_done
out/all_studies/.all_pre_integration_done
out/all_studies/.all_explore_done
```

Commands:

```bash
snakemake --profile profiles/local -s workflow/modules/qc_preproc/Snakefile out/all_studies/.all_raw_done
snakemake --profile profiles/local -s workflow/modules/qc_preproc/Snakefile out/all_studies/.all_pre_integration_done
snakemake --profile profiles/local -s workflow/modules/qc_preproc/Snakefile out/all_studies/.all_explore_done
```

## Output highlights

Per study:

- `out/<study>/qc/pre_filter_metrics.csv`
- `out/<study>/qc/pictures/pre_filter_genes_per_cell.png`
- `out/<study>/qc/pictures/pre_filter_counts_per_cell.png`
- `out/<study>/qc/pictures/pre_filter_mito_pct.png`
- `out/<study>/qc/<study>_pre_integration.h5ad`
- `out/<study>/qc/<study>_explore.h5ad`

All studies:

- `out/all_studies/all_raw.h5ad`
- `out/all_studies/all_pre_integration.h5ad`
- `out/all_studies/all_explore.h5ad`

## Notes

- Per-study raw inputs are expected at `out/<study>/raw.h5ad` (typically from the `standardize` module).
- Per-study explore mode is controlled by `qc.study_explore_with_integration` in `config/qc.yaml`:
  - `true` -> integration mode
  - `false` -> no-integration mode
