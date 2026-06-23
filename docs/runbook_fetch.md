# fetch module

Downloads files listed in `config/liver_scRNAseq_metadata_everything.csv`.
Run from `a1_p1_qc/`.

## Commands

```bash
# dry-run all studies
snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile -n -p fetch_all

# fetch all studies
snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile fetch_all

# fetch one study
snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile \
  out/<study_accession>/.fetch_done

# fetch one run
snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile \
  out/<study_accession>/fastqs/<donor_id>/<sample_id>/<run_id>/.fetch_done
```

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
