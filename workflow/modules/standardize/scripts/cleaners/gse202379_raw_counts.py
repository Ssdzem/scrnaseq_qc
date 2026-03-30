from __future__ import annotations

import re
from pathlib import Path
from typing import TextIO

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from rtlog import log_kv, log_line

from cleaners.tenx_h5_multi import run_tenx_h5_multi

RAW_COUNTS_RE = re.compile(r"^(?P<gsm>GSM\d+)_B-(?P<run_id>[^-]+)-raw_counts\.csv(?:\.gz)?$")


def _require_path(value: str | None, name: str) -> Path:
    if not value:
        raise SystemExit(f"Missing required input '{name}'")
    path = Path(value)
    if not path.exists():
        raise SystemExit(f"Input '{name}' not found: {path}")
    return path


def _pick_raw_count_files(raw_dir: Path) -> list[tuple[str, Path]]:
    grouped: dict[str, list[Path]] = {}
    for p in sorted(raw_dir.iterdir()):
        if not p.is_file():
            continue
        match = RAW_COUNTS_RE.match(p.name)
        if not match:
            continue
        run_id = match.group("run_id")
        grouped.setdefault(run_id, []).append(p)

    selected: list[tuple[str, Path]] = []
    for run_id in sorted(grouped):
        options = grouped[run_id]
        plain_csv = [p for p in options if p.suffix == ".csv"]
        chosen = plain_csv[0] if plain_csv else options[0]
        selected.append((run_id, chosen))
    return selected


def _looks_like_gene_symbol(value: str) -> bool:
    v = value.strip()
    if not v:
        return False
    if v.upper().startswith("ENSG"):
        return True
    if v.upper().startswith("MT-"):
        return True
    return any(c.isalpha() for c in v)


def _infer_needs_transpose(df: pd.DataFrame) -> bool:
    if df.shape[0] > df.shape[1]:
        return True
    if df.shape[0] < df.shape[1]:
        return False

    row_examples = [str(v) for v in df.index[:20]]
    col_examples = [str(v) for v in df.columns[:20]]
    row_gene_like = sum(_looks_like_gene_symbol(v) for v in row_examples)
    col_gene_like = sum(_looks_like_gene_symbol(v) for v in col_examples)
    return row_gene_like >= col_gene_like


def _csv_to_adata(*, csv_path: Path, sample_id: str, log: TextIO) -> ad.AnnData:
    log_line(log, f"[{sample_id}] Reading raw counts CSV: {csv_path}")
    df = pd.read_csv(csv_path, index_col=0)
    if df.empty:
        raise SystemExit(f"Raw counts CSV is empty: {csv_path}")

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    if _infer_needs_transpose(df):
        cells = df.columns.astype(str)
        genes = df.index.astype(str)
        matrix = sparse.csr_matrix(df.to_numpy(dtype=np.float32).T)
    else:
        cells = df.index.astype(str)
        genes = df.columns.astype(str)
        matrix = sparse.csr_matrix(df.to_numpy(dtype=np.float32))

    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=pd.Index(cells, name="barcode")),
        var=pd.DataFrame(index=pd.Index(genes, name="gene_symbol")),
    )
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    log_line(log, f"[{sample_id}] Built AnnData in-memory: n_obs={adata.n_obs} n_vars={adata.n_vars}")
    return adata


def _write_10x_h5_from_adata(*, adata: ad.AnnData, out_h5: Path, sample_id: str, log: TextIO) -> None:
    out_h5.parent.mkdir(parents=True, exist_ok=True)

    # 10x expects genes x cells, CSC compressed.
    matrix = sparse.csc_matrix(adata.X.T)
    if matrix.nnz and not np.issubdtype(matrix.data.dtype, np.integer):
        matrix.data = np.rint(matrix.data).astype(np.int64, copy=False)
    elif matrix.nnz:
        matrix.data = matrix.data.astype(np.int64, copy=False)

    barcodes = np.asarray(adata.obs_names.astype(str), dtype="S")
    genes = np.asarray(adata.var_names.astype(str), dtype="S")
    feature_type = np.asarray(["Gene Expression"] * adata.n_vars, dtype="S")
    genome = np.asarray(["unknown"] * adata.n_vars, dtype="S")

    with h5py.File(out_h5, "w") as h5:
        matrix_grp = h5.create_group("matrix")
        matrix_grp.create_dataset("barcodes", data=barcodes)
        matrix_grp.create_dataset("data", data=matrix.data)
        matrix_grp.create_dataset("indices", data=matrix.indices.astype(np.int32, copy=False))
        matrix_grp.create_dataset("indptr", data=matrix.indptr.astype(np.int64, copy=False))
        matrix_grp.create_dataset(
            "shape",
            data=np.asarray([adata.n_vars, adata.n_obs], dtype=np.int64),
        )

        features_grp = matrix_grp.create_group("features")
        features_grp.create_dataset("id", data=genes)
        features_grp.create_dataset("name", data=genes)
        features_grp.create_dataset("feature_type", data=feature_type)
        features_grp.create_dataset("genome", data=genome)

    log_line(log, f"[{sample_id}] Wrote intermediate 10x h5: {out_h5}")


def run_gse202379_raw_counts(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    inp = study_spec.get("input", {})
    raw_dir = _require_path(inp.get("raw_dir"), "raw_dir")
    tmp_h5_dir = output_h5ad.parent / ".tmp" / f"{study}_raw_count_h5"
    tmp_h5_dir.mkdir(parents=True, exist_ok=True)

    run_files = _pick_raw_count_files(raw_dir)
    if not run_files:
        raise SystemExit(
            f"No files matched '*-raw_counts.csv(.gz)' pattern under {raw_dir}"
        )

    log_kv(log, "gse202379_raw_counts.start", study=study, raw_dir=raw_dir, samples=len(run_files))

    sample_entries: list[dict[str, str]] = []
    for run_id, csv_path in run_files:
        sample_id = run_id
        adata = _csv_to_adata(csv_path=csv_path, sample_id=sample_id, log=log)
        sample_h5 = tmp_h5_dir / f"{sample_id}.h5"
        _write_10x_h5_from_adata(adata=adata, out_h5=sample_h5, sample_id=sample_id, log=log)
        sample_entries.append({"sample_id": sample_id, "h5": str(sample_h5)})

    delegated_spec = {"input": {"samples": sample_entries, "qc": inp.get("qc"), "metadata": inp.get("metadata")}}
    run_tenx_h5_multi(study=study, study_spec=delegated_spec, output_h5ad=output_h5ad, log=log)
    log_kv(log, "gse202379_raw_counts.complete", study=study, output=output_h5ad)
