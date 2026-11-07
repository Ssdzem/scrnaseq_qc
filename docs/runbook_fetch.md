# fetch module

Downloads files listed in `config/samples.csv`.
Run from the repository root.

## Commands

```bash
# dry-run all studies
snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile -n -p fetch_all

# fetch all studies
snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile fetch_all

# fetch one study
# !! DO NOT USE the form below — it is a SILENT NO-OP for studies with no FASTQs.
#    Use `fetch_all` above. See the warning after this block.
#    snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile \
#      out/<study_accession>/.fetch_done

# fetch one run
snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile \
  out/<study_accession>/fastqs/<donor_id>/<sample_id>/<run_id>/.fetch_done
```

> **`out/<study>/.fetch_done` downloads nothing for studies with no FASTQs.** That relative
> path reaches a rule with no inputs that only touches the sentinel. The rule that fetches
> emits the same file at an absolute path, so only `fetch_all` reaches it. Because
> standardize uses this sentinel as its fetch dependency, running a single
> `out/<study>/raw.h5ad` can also pass without fetching, unless the study is listed in a
> profile's `extract:` map. Run `fetch_all` first. This is a known open defect.

## Config

Primary file: `config/fetch.yaml`.

Important knobs:

- `paths.metadata_csv`: source table.
- `paths.out_root`: output root, normally `out`.
- `aspera.enabled`: use Aspera for GSA/CNCB-style downloads.
- `ASCP_PATH` / `ASCP_KEY`: optional env overrides for Aspera.

## Outputs

- Study done marker: `out/<study>/.fetch_done`
- Run done marker: `out/<study>/fastqs/<donor>/<sample>/<run>/.fetch_done`
