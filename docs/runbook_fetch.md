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
# !! DO NOT USE the form below — it is a SILENT NO-OP for GEO studies (P1-FETCH-01).
#    Use `fetch_all` above. See the warning after this block.
#    snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile \
#      out/<study_accession>/.fetch_done

# fetch one run
snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile \
  out/<study_accession>/fastqs/<donor_id>/<sample_id>/<run_id>/.fetch_done
```

> **P1-FETCH-01 — `out/<study>/.fetch_done` is a silent no-op for GEO studies.**
>
> That relative path routes to `fetch_done_data` (`fetch.smk:305`), whose input list is
> `TRIPLES_BY_STUDY_FASTQS.get(study, set())` — **empty for every non-FASTQ study** — and
> whose shell body is `"true"`. It touches the sentinel and downloads nothing. The rule
> that actually fetches, `fetch_done` (`fetch.smk:218`), emits the same file at an
> *absolute* path, so only `fetch_all` reaches `download-study`.
>
> This matters beyond the runbook: `standardize.smk:138-139` returns
> `out/{study}/.fetch_done` as the fetch dependency for any study whose inputs sit under
> `out/<study>/raw` or `/qconly`, so `snakemake out/GSE202379/raw.h5ad` satisfies its
> fetch dependency **without fetching**. The validation profile escapes it, because
> `fetch_done_data`'s wildcard constraint (`fetch.smk:314`) excludes `extract:` studies.
>
> **Still open.** The fix is a code change (make `OUT_ROOT` relative, or drop
> `fetch_done_data`'s no-op branch and have it depend on `fetch_done` for studies with no
> FASTQ triples). Deliberately not made here: it alters a download path that cannot be
> exercised on this machine, and an untested change to fetching is worse than a documented
> defect. This runbook previously documented the broken form as the way to fetch a study.

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
