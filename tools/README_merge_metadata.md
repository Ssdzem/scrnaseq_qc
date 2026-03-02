# Merge study metadata into main metadata

Goal: enrich `config/liver_scRNAseq_metadata.csv` with selected donor-level fields
(currently: `age_years`, `sex`, `bmi`) from per-study metadata CSVs in
`config/`.

## How it works

- Each study in `liver_scRNAseq_metadata.csv:study_accession` is expected to
  have a corresponding file `config/{study_accession}_metadata.csv`.
- Some studies may also have `config/{study_accession}_sup_metadata.csv` (used
  as an additional source of rows).
- The join key is derived from `liver_scRNAseq_metadata.csv:donor_id` and is
  study-specific.

## Usage

Write a new merged CSV (recommended):

```bash
python projects/a1p1_qc/tools/merge_study_metadata_into_main.py \
  --out projects/a1p1_qc/config/liver_scRNAseq_metadata.merged.csv
```

Overwrite the main CSV in place:

```bash
python projects/a1p1_qc/tools/merge_study_metadata_into_main.py --in-place
```

Only merge specific studies:

```bash
python projects/a1p1_qc/tools/merge_study_metadata_into_main.py \
  --studies GSE202379 HRA007511 \
  --out /tmp/liver_scRNAseq_metadata.partial.csv
```

By default, existing non-empty values in the main CSV are preserved. To force
overwriting:

```bash
python projects/a1p1_qc/tools/merge_study_metadata_into_main.py \
  --overwrite-existing \
  --out /tmp/liver_scRNAseq_metadata.overwrite.csv
```

