# standardize module runbook

This project uses Snakemake. This runbook mirrors the style of `docs/runbook_fetch.md` and documents how to run the study standardization module.

Everything below assumes you run from the `a1_p1_qc/` directory.

## Quick start

Dry-run all standardization outputs discovered from `config/standardize/*.yaml`:

```bash
snakemake --profile profiles/local -s workflow/modules/standardize/Snakefile -n -p standardize_all
```

## Snakemake profile

Use the local profile directly:

```bash
snakemake --profile profiles/local -s workflow/modules/standardize/Snakefile <target>
```

## Study specifications

The module discovers studies from:

```text
config/standardize/*.yaml
```

For each study spec `<study>.yaml`, canonical output is:

```text
out/<study>/raw.h5ad
```

## Targets

### Standardize all discovered studies

Target:

```text
standardize_all
```

Command:

```bash
snakemake --profile profiles/local -s workflow/modules/standardize/Snakefile standardize_all
```

### Standardize a single study

Target:

```text
out/<study>/raw.h5ad
```

Command:

```bash
snakemake --profile profiles/local -s workflow/modules/standardize/Snakefile out/<study>/raw.h5ad
```

### Run unpack/validation stages explicitly (debugging)

Unpack target:

```text
out/<study>/.unpack_done
```

Validation target:

```text
out/<study>/.standardize_inputs_validated_done
```

Commands:

```bash
snakemake --profile profiles/local -s workflow/modules/standardize/Snakefile out/<study>/.unpack_done
snakemake --profile profiles/local -s workflow/modules/standardize/Snakefile out/<study>/.standardize_inputs_validated_done
```

## Execution flow for each study

1. `standardize_unpack_study`
   - Requires `out/<study>/.fetch_done` and `config/standardize/<study>.yaml`.
   - Runs `workflow/modules/standardize/scripts/unpack_inputs.py`.
2. `standardize_validate_inputs`
   - Verifies all required input file paths from the study spec exist.
3. `standardize_study`
   - Runs `workflow/modules/standardize/scripts/standardize.py`.
   - Produces `out/<study>/raw.h5ad`.

## Conda environment selection

`standardize_study` chooses env by cleaners listed in `pipeline`:

- if any cleaner is `env="r"` in the cleaner registry -> `workflow/modules/standardize/env/standardize_r.yaml`
- otherwise -> `workflow/modules/standardize/env/standardize_py.yaml`

## Notes

- Ensure `fetch` has completed for each study (`out/<study>/.fetch_done`) before running standardization.
- Keep one YAML per study under `config/standardize/` so `standardize_all` can discover it automatically.
