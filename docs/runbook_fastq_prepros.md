# fastq_prepros module

Runs CeleScope counting for configured singleton FASTQ studies. Run from `a1_p1_qc/`.

## Commands

```bash
# dry-run configured targets
snakemake --profile profiles/local -s workflow/modules/fastq_prepros/Snakefile -n -p fastq_prepros_all

# execute configured targets
snakemake --profile profiles/local -s workflow/modules/fastq_prepros/Snakefile fastq_prepros_all

# optional reference build
snakemake --profile profiles/local -s workflow/modules/fastq_prepros/Snakefile celescope_ref_all
```

## Config

Primary file: `config/fastq_prepros.yaml`.

Important knobs:

- `fastq_prepros.studies`: studies to count.
- `fastq_prepros.target_level`: `run` or `sample`.
- `fastq_prepros.only_existing_inputs`: skip missing FASTQ folders when `true`.
- `fastq_prepros.celescope.threads`: Snakemake threads and CeleScope threads.
- `fastq_prepros.celescope.mem_gb`: Snakemake `mem_mb` resource.
- `fastq_prepros.celescope.star_index`: required STAR index.

## Outputs

```text
out/<study>/counts/singleton/<run|sample>/<donor>/<sample>/<target>/.count_done
```
