# a1_p1_qc workflow wiring

Run from `a1_p1_qc/`. The root `Snakefile` includes:

```text
fetch -> fastq_prepros -> standardize -> qc_preproc -> integration_annotation
```

The default `all` target currently drives fetch, FASTQ preprocessing, and QC. Standardization and integration/annotation are explicit module targets because they are often run/debugged separately.

## Module entrypoints

| Module | Snakefile | Main target | Primary output |
| --- | --- | --- | --- |
| fetch | `workflow/modules/fetch/Snakefile` | `fetch_all` | `out/<study>/.fetch_done` |
| fastq_prepros | `workflow/modules/fastq_prepros/Snakefile` | `fastq_prepros_all` | `out/<study>/counts/.../.count_done` |
| standardize | `workflow/modules/standardize/Snakefile` | `standardize_all` | `out/<study>/raw.h5ad` |
| qc_preproc | `workflow/modules/qc_preproc/Snakefile` | `qc_all` | `out/all_studies/all_pre_integration.h5ad` |
| integration_annotation | `workflow/modules/integration_annotation/Snakefile` | `integration_annotation_all` | `out/integration_annotation/all_studies/all_annotated.h5ad` |

## Profile/config resolution

Profiles set `profile_config_dir=profiles/<profile>`. Module configs resolve in this order:

1. explicit Snakemake config key;
2. profile-local YAML, if present;
3. default YAML under `config/`.

See `profiles/README.md` for CPU/GPU profiles and `docs/runbook_*.md` for exact commands.
