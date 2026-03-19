# Fetch module runbook

This project uses Snakemake. To keep an audit trail of commands (especially on remote/HPC machines), prefer running via:

- a Snakemake profile (`profiles/...`) to record *execution settings* (cores, conda, rerun policy, etc.)
- a small Makefile (`Makefile`) to record *common entrypoints/targets*

Everything below assumes you run from the `a1_p1_qc/` directory.

## Quick start

Show available Makefile targets:

```bash
make
```

## Snakemake profiles

The default local profile is `profiles/local/config.yaml`.

Use it directly:

```bash
snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile -n -p <target>
```

## Fetch targets (study / run / all)

### Fetch a whole study

This downloads all assets referenced in the metadata CSV for that study.

Target file:

```text
out/<study_accession>/.fetch_done
```

Command:

```bash
snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile \
  out/<study_accession>/.fetch_done
```

Makefile wrapper:

```bash
make fetch-study STUDY=<study_accession> CORES=6
```

### Fetch a single run

This is intended for testing and for parallelization (e.g. `--cores 6` to fetch 6 runs simultaneously).

Target file:

```text
out/<study_accession>/fastqs/<donor_id>/<sample_id>/<run_id>/.fetch_done
```

Command:

```bash
snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile \
  out/<study_accession>/fastqs/<donor_id>/<sample_id>/<run_id>/.fetch_done
```

Makefile wrapper:

```bash
make fetch-run STUDY=<study_accession> DONOR=<donor_id> SAMPLE=<sample_id> RUN=<run_id> CORES=1
```

### Fetch all studies in the metadata

Command:

```bash
snakemake --profile profiles/local -s workflow/modules/fetch/Snakefile fetch_all
```

Makefile wrapper:

```bash
make fetch-all CORES=6
```

## Curl vs Aspera (ascp)

The fetch module supports two download modes:

- Default: `curl` (URLs in metadata CSV)
- Optional: Aspera Connect `ascp` for CNCB/GSA studies like `HRA...`.

To enable Aspera:

1) Ensure `ascp` is available, via one of:

- Put `ascp` on your `PATH` (recommended), or
- Set `ASCP_PATH` to the full path of the `ascp` binary, or
- Set `aspera.ascp_path` in `config/fetch.yaml`.

The resolution precedence is:

```text
ASCP_PATH env var > aspera.ascp_path config > `which ascp`
```

2) Provide the Aspera SSH key, via one of:

- Set `ASCP_KEY` to the full path of the key file (recommended on shared systems), or
- Set `aspera.key_path` in `config/fetch.yaml`.

The resolution precedence is:

```text
ASCP_KEY env var > aspera.key_path config
```

3) Edit `config/fetch.yaml` and set:

```yaml
aspera:
  enabled: true
  # Optional if `ascp` is on PATH or ASCP_PATH is set.
  ascp_path: ""
  # Optional if ASCP_KEY is set.
  key_path: ""
```

4) Then run the same Snakemake targets as above.

Notes:
- Aspera mode downloads per-run directories from:
  `aspera.host:aspera.remote_root/<study>/<run_id>/`
- Destination directories are created before `ascp` runs.
