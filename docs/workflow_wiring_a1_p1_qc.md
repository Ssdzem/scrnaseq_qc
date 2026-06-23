# a1_p1_qc modular Snakemake wiring: rules, inputs/outputs, and script workings

## 1) Top-level workflow wiring

- Root `Snakefile` includes modules in this order: `fetch`, `fastq_prepros`, `standardize`, `qc_preproc`.
- Default target `rule all` requires:
  - `fetch_all.input`
  - `fastq_prepros_all.input`
  - `qc_all.input`
- Single-study target `qc_preintegration_study` points to `out/{study}/qc/.study_explore_done`.

Practical branch behavior:
- The **fetch** module always runs.
- The **fastq_prepros** module only produces targets for studies listed in `config/fastq_prepros.yaml` and for rows with available fastq directories / metadata after filtering.
- The **standardize** module is not in `rule all`, but `qc_preproc` consumes `out/{study}/raw.h5ad`; in practice standardization must complete first for each study being QC’d.
- If a study has no FASTQs, you can omit it from `fastq_prepros.studies`; then it effectively follows fetch -> standardize -> qc_preproc.

---

## 2) Module: fetch

Config + path assumptions:
- Loads `config/fetch.yaml` (project-relative) and requires execution from `a1_p1_qc/`.
- Reads metadata CSV via `paths.metadata_csv`.
- Writes outputs under `paths.out_root`.

### Rules

1. `fetch_all`
- **Input:** expanded `out/{study}/.fetch_done` for all studies parsed from metadata.
- **Output:** none directly (aggregator target).
- **Work:** drives per-study completion.

2. `fetch_done`
- **Input:** metadata CSV, `workflow/modules/fetch/scripts/fetch_data.py`.
- **Output:** `out/{study}/.fetch_done` (touch).
- **Script/work:** executes `fetch_data.py ... download-study {study}` (curl mode or Aspera mode depending on config).

3. `fetch_study`
- **Input:** `out/{study}/.fetch_done`.
- **Output:** none directly.
- **Work:** convenience target alias.

4. `fetch_run`
- **Input:** metadata CSV, `fetch_data.py`.
- **Output:** `out/{study}/fastqs/{run_id}/.fetch_done` (touch).
- **Script/work:** executes `fetch_data.py ... download-run {study} {run_id}`.

5. `fetch_done_data`
- **Input:** metadata CSV, `fetch_data.py`.
- **Output:** `out/{study}/.fetch_done` (touch, explicit `out/` pattern variant).
- **Script/work:** same as `fetch_done`, forces `--out-root out`.

6. `fetch_run_data`
- **Input:** metadata CSV, `fetch_data.py`.
- **Output:** `out/{study}/fastqs/{donor_id}/{sample_id}/{run_id}/.fetch_done`.
- **Script/work:** same run download command, explicit nested `out/` target shape.

### `fetch_data.py` behavior summary

- Parses metadata rows and validates required columns: `study_accession`, `donor_id`, `sample_id`, `run_id`.
- Supports commands: `list-studies`, `list-studies-json`, `download-study`, `download-run`.
- Download sources by metadata columns:
  - FASTQs: `fastq_r1`, `fastq_r2` -> `fastqs/`
  - Raw counts: `raw_counts` -> `raw/`
  - Processed counts: `processed_counts` -> `qconly/`
- Curl mode:
  - Resolves destination hierarchy level (run/sample/donor/study) from URL token matching + fallback uniqueness heuristic.
  - Downloads with robust curl resume/retry flags.
- Aspera mode:
  - Validates ascp/key/host/remote-root.
  - Downloads per-run directories with retry schedule and stall detection/termination logic.

---

## 3) Module: fastq_prepros

Config + assumptions:
- Loads `config/fastq_prepros.yaml`, requires execution from `a1_p1_qc/`.
- Filters metadata to `include` truthy rows and configured studies.
- Currently supports platform `singleton` only.
- Requires `fastq_prepros.celescope.star_index`.

### Rules (`rules/fastq_prepros.smk`)

1. `fastq_prepros_all`
- **Input:** all expected `.count_done` files for computed targets (run-level or sample-level).
- **Output:** none directly.
- **Work:** aggregator.

2. `celescope_count`
- **Input:** paired FASTQs from `_input_fastqs` resolver (from `out/{study}/fastqs/{donor}/{sample}/{run}`)
- **Output:** `out/{study}/counts/singleton/{target_level}/{donor}/{sample}/{target_id}/.count_done`
- **Script/work:** runs `run_celescope_multi_rna.py` with outdir/sample/fastqs/genome-dir/threads.

### Reference build rules (`rules/ref_singleton_hs_ensembl_110.smk`)

1. `celescope_ref_all`
- **Input:** reference done marker.

2. `download_ensembl_fasta_gtf`
- **Output:** reference FASTA+GTF gz files.
- **Script:** `ref_download_ensembl.py` (URL download with temp file).

3. `gunzip_ensembl`
- **Input:** gz FASTA+GTF.
- **Output:** uncompressed FASTA+GTF.
- **Script:** `ref_gunzip.py`.

4. `ensure_mt_gene_list`
- **Output:** `mt_gene_list.txt`.
- **Script:** `ref_write_mt_gene_list.py` (writes `MT-` and `mt-`).

5. `celescope_mkgtf`
- **Input:** GTF.
- **Output:** filtered GTF.
- **Script:** `ref_mkgtf.py` (`celescope utils mkgtf`).

6. `celescope_mkref`
- **Input:** FASTA, filtered GTF, mt list.
- **Output:** reference done marker.
- **Script:** `ref_mkref.py` (`celescope rna mkref`).

### `run_celescope_multi_rna.py` behavior summary

- Normalizes FASTQ filenames (`*_f1.fastq.gz`/`*_r2.fastq.gz`) to CeleScope-friendly `*_1.fq.gz`/`*_2.fq.gz` in place.
- Writes one-line `rna.mapfile`.
- Executes `multi_rna` with supplied genome dir and threads.
- In shell mode, executes generated `shell/{sample}.sh`.

---

## 4) Module: standardize

Config/selection:
- Auto-discovers studies from `config/standardize/*.yaml`.
- Canonical output per study: `out/{study}/raw.h5ad`.
- Cleaner pipeline is read from study spec `pipeline`/`Pipeline` and validated against registry.

### Rules

1. `standardize_all`
- **Input:** `out/{study}/raw.h5ad` for all discovered studies.

2. `standardize_unpack_study`
- **Input:** study spec yaml, `out/{study}/.fetch_done`.
- **Output:** `out/{study}/.unpack_done`.
- **Script/work:** `unpack_inputs.py` unpacks tar/tgz/gz under `out/{study}/raw` and `out/{study}/qconly` (skips fastqs).

3. `standardize_validate_inputs`
- **Input:** spec yaml + unpack done marker.
- **Output:** `out/{study}/.standardize_inputs_validated_done`.
- **Work:** Python `run:` block computes required input paths for active cleaners and fails if any are missing.

4. `standardize_study`
- **Input:** spec yaml + validated marker.
- **Output:** `out/{study}/raw.h5ad`.
- **Script/work:** `standardize.py` orchestrates cleaner pipeline, step-wise intermediate `.tmp` h5ads for non-final steps, chooses R or Python conda env depending on cleaner registry.

### Standardize orchestration scripts

- `cleaner_registry.py`: declares cleaner names, required path keys, and env type (`py`/`r`).
- `unpack_inputs.py`: safe tar extraction guard (prevents path traversal), optional symlink convenience for nested tar layouts, gzip detection + decompression.
- `standardize.py`:
  - Parses study spec.
  - Builds pipeline order.
  - Dispatches to cleaner functions.
  - Handles special local spec shaping for `metadata_join` and `tenx_h5_multi` chaining.
- `rtlog.py`: timestamped streaming logger + subprocess streamer for long-running steps.

### Cleaner scripts and what they do

1. `count_matrix_10x.py`
- Validates triplet paths (`matrix`, `genes`, `barcodes`).
- Normalizes to gzipped 10x files (`matrix.mtx.gz`, `genes.tsv.gz`, `barcodes.tsv.gz`), including CSV->TSV conversions and header-count disambiguation.
- Ensures both `genes.tsv.gz` and `features.tsv.gz` compatibility.
- Loads via `scanpy.read_10x_mtx`, optionally joins metadata by barcode/cell column, writes h5ad.

2. `tenx_h5_multi.py`
- Accepts list of samples each with `h5` or `h5ad`.
- Optional per-sample QC (MAD-based thresholds, mitochondrial cap, miQC-like GMM filter, scrublet doublet annotation).
- Concatenates samples, optional metadata join, writes merged h5ad + qc metadata in `.uns`.

3. `metadata_join.py`
- Loads input h5ad + metadata CSV.
- Joins metadata on configured `join_column` with optional `obs_key` override.
- Writes updated h5ad.

4. `seurat_rds.py`
- Runs temporary R script (`Seurat` + `zellkonverter`) converting RDS/SCE to h5ad.
- Reads converted file in Python, annotates `.uns`, writes output h5ad.

5. `csv_handler.py`
- Detects `*-raw_counts.csv(.gz)` files.
- Converts each to AnnData then to intermediate 10x-style H5.
- Stores generated sample list into spec input for downstream processing and writes placeholder h5ad.

6. `count_matrix_10x_multi_runs.py`
- Loads a run table from CSV and filters included rows for a study.
- Standardizes each run via `count_matrix_10x` and prepares `input.samples` for downstream `tenx_h5_multi`.
- Supports configurable run-table column names and triplet layout template.

---

## 5) Module: qc_preproc

Config + study list:
- Uses `config/qc.yaml`.
- Study list from `inputs.studies` or inferred from metadata CSV.

### Rules

1. `qc_all`
- **Input:** per-study done markers (`.raw_metadata_done`, `.pretrim_pictures_done`, `.study_explore_done`) and all-studies markers (`.all_raw_done`, `.all_explore_done`) plus qc config.
- **Output:** `out/all_studies/.qc_all_done`.

2. `qc_attach_metadata`
- **Input:** `out/{study}/raw.h5ad`, qc config.
- **Output:** `out/{study}/raw_metadata.h5ad`, `out/{study}/qc/.raw_metadata_done`.
- **Script:** `attach_metadata.py`.

3. `qc_pretrim_pictures`
- **Input:** `raw_metadata.h5ad`, qc config, `.raw_metadata_done`.
- **Output:** pre-filter metrics CSV + three histogram PNGs + `.pretrim_pictures_done`.
- **Script:** `qc_figures.py`.

4. `qc_study_pre_integration`
- **Input:** `raw_metadata.h5ad`, qc config, `.pretrim_pictures_done`.
- **Output:** `{study}_pre_integration.h5ad`, `.study_pre_integration_done`.
- **Script:** `explore_scanpy.py --mode pre`.

5. `qc_study_explore`
- **Input:** `{study}_pre_integration.h5ad`, qc config, `.study_pre_integration_done`.
- **Output:** `{study}_explore.h5ad`, `.study_explore_done`.
- **Script:** `explore_scanpy.py --mode integration|no_integration` (from config flag).

6. `qc_concat_all_raw`
- **Input:** qc config + all studies `out/{study}/raw_metadata.h5ad`.
- **Output:** `out/all_studies/all_raw.h5ad`, `.all_raw_done`.
- **Script:** `concat_raw.py`.

7. `qc_all_pre_integration`
- **Input:** `all_raw.h5ad`, qc config, `.all_raw_done`.
- **Output:** `all_pre_integration.h5ad`, `.all_pre_integration_done`.
- **Script:** `explore_scanpy.py --mode pre`.

8. `qc_all_explore`
- **Input:** `all_pre_integration.h5ad`, qc config, `.all_pre_integration_done`.
- **Output:** `all_explore.h5ad`, `.all_explore_done`.
- **Script:** `explore_scanpy.py --mode integration`.

### QC script behavior summary

- `attach_metadata.py`: merges metadata CSV columns into `adata.obs` by `sample_id`.
- `qc_figures.py`: computes QC metrics (`n_genes_by_counts`, `total_counts`, `pct_counts_mt`), writes one-row metrics summary, and optional pre-filter histograms.
- `explore_scanpy.py`:
  - `pre` mode: normalization/log1p/HVG preparation and bookkeeping.
  - `integration` mode: scVI or Harmony integration + neighbors/UMAP/Leiden.
  - `no_integration` mode: uses PCA representation for downstream embedding/clustering.
- `concat_raw.py`: concatenates study-level `raw_metadata.h5ad` into all-studies matrix.
- `qc_scanpy.py` exists (basic filter workflow) but is not wired by current `qc_preproc.smk` rules.

---

## 6) End-to-end dependency chain (what to run)

For your requested sequence (fetch -> conditional fastq_prepros -> standardize -> qc_preproc):

1. Fetch study assets (`fetch_all` / `fetch_study`).
2. If FASTQs exist and study is configured for fastq pre-processing, run `fastq_prepros_all` (and optionally `celescope_ref_all` once per reference).
3. Run `standardize_all` (or per-study `standardize_study`) to produce `out/{study}/raw.h5ad`.
4. Run `qc_all` (or per-study chain ending in `out/{study}/qc/.study_explore_done`).

Suggested command pattern:
- `snakemake -s Snakefile fetch_all`
- `snakemake -s Snakefile fastq_prepros_all` (only for applicable studies)
- `snakemake -s Snakefile standardize_all`
- `snakemake -s Snakefile qc_all`

Or run selectively per study with explicit targets if you want strict branch control.

---

## 7) Concrete run recipe for studies GSE202379, GSE244832, GSE256398

### Parameters to tweak first

1. `config/fetch.yaml`
- Set `aspera.enabled: false` (your noted workaround).

2. `config/qc.yaml`
- Set:
  ```yaml
  inputs:
    studies: ["GSE202379", "GSE244832", "GSE256398"]
  ```
  This is required so QC aggregation (`qc_concat_all_raw`, `qc_all_pre_integration`, `qc_all_explore`) only uses these three studies.

3. `config/qc.yaml` integration knobs (optional)
- Keep/adjust:
  - `qc.integration_method` (`scvi` or `harmony`)
  - `qc.study_explore_with_integration` (per-study explore behavior)
  - `qc.scvi_max_epochs`, `qc.hvg_n_top_genes`, `qc.n_pcs`, `qc.n_neighbors`

4. `config/fastq_prepros.yaml`
- No change needed for these three GEO studies in the current repo wiring.
- Keep `fastq_prepros.studies` as-is (currently `HRA007511`) unless you explicitly want to run FASTQ->counting for another study.

### Commands to run

Run from repo root `a1_p1_qc/`.

1. Fetch only the three studies:
```bash
snakemake -s Snakefile --use-conda --cores 8 \
  out/GSE202379/.fetch_done \
  out/GSE244832/.fetch_done \
  out/GSE256398/.fetch_done
```

2. Standardize to per-study `raw.h5ad`:
```bash
snakemake -s Snakefile --use-conda --cores 8 \
  out/GSE202379/raw.h5ad \
  out/GSE244832/raw.h5ad \
  out/GSE256398/raw.h5ad
```

3. QC + concatenate + pre-integration + integration on concatenated object:
```bash
snakemake -s Snakefile --use-conda --cores 8 out/all_studies/all_explore.h5ad
```

If you want to stop exactly at concatenated pre-integration (before integration embedding), use:
```bash
snakemake -s Snakefile --use-conda --cores 8 out/all_studies/all_pre_integration.h5ad
```
