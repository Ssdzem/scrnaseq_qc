# scrnaseq_qc

Turns single-cell count data into QC-checked, integrated and annotated AnnData objects.
Originally extracted from the lSMIF project.

```text
fetch -> fastq_prepros -> standardize -> qc_preproc -> integration_annotation
```

| Module | Does | Main output |
|---|---|---|
| `fetch` | Downloads study files | `out/<study>/.fetch_done` |
| `fastq_prepros` | Counts FASTQs with CeleScope | Count tree under `out/<study>/counts/` |
| `standardize` | Converts each study to one object and filters cells if `input.qc` is set | `out/<study>/raw.h5ad` |
| `qc_preproc` | Attaches metadata, normalizes, selects variable genes | `out/<study>/qc/<study>_pre_integration.h5ad` |
| `integration_annotation` | Integrates (sysVI, scVI or Harmony) and transfers labels (scANVI or Symphony) | `*_annotated.h5ad` |

Run from the repository root, one module at a time:

```bash
snakemake --profile profiles/local -s workflow/modules/standardize/Snakefile -n -p standardize_all
```

Nothing here is tied to a project. Everything specific to a run belongs in a profile kept
outside this repository: studies, thresholds, reference atlas, memory, CPU or GPU. Start
from `config/templates/samples.csv`, `config/templates/standardize_study.yaml` and the
templates in `profiles/`. See `profiles/README.md` for how profiles work and `docs/README.md`
for each module's runbook. Conda environments are not fully pinned, so record the exact
solve of each run.
