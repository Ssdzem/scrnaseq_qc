# fastq_prepros module runbook

This project uses Snakemake. This runbook mirrors the style of `docs/runbook_fetch.md` and documents how to run the FASTQ preprocessing module (CeleScope counting + optional reference build).

Everything below assumes you run from the `a1_p1_qc/` directory.

## Quick start

Dry-run all configured preprocessing targets:

```bash
snakemake --profile profiles/local -s workflow/modules/fastq_prepros/Snakefile -n -p fastq_prepros_all
```

## Snakemake profile

Use the local profile directly:

```bash
snakemake --profile profiles/local -s workflow/modules/fastq_prepros/Snakefile <target>
```

## Configuration inputs

Main config file:

```text
config/fastq_prepros.yaml
```

Important knobs:

- `fastq_prepros.studies`: studies to process.
- `fastq_prepros.target_level`: `run` or `sample`.
- `fastq_prepros.only_existing_inputs`: if `true`, only process runs/samples with existing FASTQ inputs in `out/<study>/fastqs/...`.
- `fastq_prepros.celescope.star_index`: required STAR index directory (produced by CeleScope `mkref`).

## Targets

### Run preprocessing for all configured studies

Target:

```text
fastq_prepros_all
```

Command:

```bash
snakemake --profile profiles/local -s workflow/modules/fastq_prepros/Snakefile fastq_prepros_all
```

Output pattern (one `.count_done` per run or sample, based on `target_level`):

```text
out/<study>/counts/singleton/<target_level>/<donor_id>/<sample_id>/<target_id>/.count_done
```

### Run preprocessing for a single run-level target

Use this when `fastq_prepros.target_level: run`.

Target:

```text
out/<study>/counts/singleton/run/<donor_id>/<sample_id>/<run_id>/.count_done
```

Command:

```bash
snakemake --profile profiles/local -s workflow/modules/fastq_prepros/Snakefile \
  out/<study>/counts/singleton/run/<donor_id>/<sample_id>/<run_id>/.count_done
```

### Run preprocessing for a single sample-level target

Use this when `fastq_prepros.target_level: sample`.

Target:

```text
out/<study>/counts/singleton/sample/<donor_id>/<sample_id>/<sample_id>/.count_done
```

Command:

```bash
snakemake --profile profiles/local -s workflow/modules/fastq_prepros/Snakefile \
  out/<study>/counts/singleton/sample/<donor_id>/<sample_id>/<sample_id>/.count_done
```

## CeleScope reference build (optional helper in same module)

The module also exposes a reference sub-workflow that downloads Ensembl FASTA/GTF and runs CeleScope `mkgtf` + `mkref`.

Top target:

```text
celescope_ref_all
```

Command:

```bash
snakemake --profile profiles/local -s workflow/modules/fastq_prepros/Snakefile celescope_ref_all
```

Done marker:

```text
data/aux/ref/hs_ensembl_110/.celescope_ref_done
```

You can customize reference URLs and output directory under:

```yaml
fastq_prepros:
  celescope:
    reference:
      ref_dir: ...
      fasta_url: ...
      gtf_url: ...
```

## Notes

- This module currently enforces Singleton platform inputs.
- FASTQs are expected under `out/<study>/fastqs/<donor_id>/<sample_id>/<run_id>/`.
- Supported paired FASTQ filename patterns include:
  - `*_f1.fastq.gz` and `*_r2.fastq.gz`
  - `*_R1_*.fastq.gz` and `*_R2_*.fastq.gz`
