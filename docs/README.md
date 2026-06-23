# a1_p1_qc runbooks

Run commands from `a1_p1_qc/`. Use Snakemake profiles for cores, conda, reruns, and module-local config overrides.

| Module | Main target | Runbook |
| --- | --- | --- |
| fetch | `fetch_all` | `runbook_fetch.md` |
| fastq_prepros | `fastq_prepros_all` | `runbook_fastq_prepros.md` |
| standardize | `standardize_all` | `runbook_standardize.md` |
| qc_preproc | `qc_all` | `runbook_qc_preproc.md` |
| integration_annotation | `integration_annotation_all` | `runbook_integration_annotation.md` |

Default profile:

```bash
snakemake --profile profiles/local -s <module>/Snakefile -n -p <target>
```

Profile behavior and CPU/GPU presets are summarized in `../profiles/README.md`.
