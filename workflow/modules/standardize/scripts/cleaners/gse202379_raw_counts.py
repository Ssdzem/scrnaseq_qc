from __future__ import annotations

import re
from pathlib import Path
from typing import TextIO

import anndata as ad
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


def _csv_to_h5ad(*, csv_path: Path, sample_id: str, out_h5ad: Path, log: TextIO) -> None:
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

    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad)
    log_line(log, f"[{sample_id}] Wrote intermediate h5ad: {out_h5ad}")


def run_gse202379_raw_counts(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    inp = study_spec.get("input", {})
    raw_dir = _require_path(inp.get("raw_dir"), "raw_dir")
    tmp_root = output_h5ad.parent / ".tmp" / f"{study}_raw_count_h5ad"
    tmp_root.mkdir(parents=True, exist_ok=True)

    run_files = _pick_raw_count_files(raw_dir)
    if not run_files:
        raise SystemExit(
            f"No files matched '*-raw_counts.csv(.gz)' pattern under {raw_dir}"
        )

    log_kv(log, "gse202379_raw_counts.start", study=study, raw_dir=raw_dir, samples=len(run_files))

    sample_entries: list[dict[str, str]] = []
    for run_id, csv_path in run_files:
        sample_id = run_id
        sample_h5ad = tmp_root / f"{sample_id}.h5ad"
        _csv_to_h5ad(csv_path=csv_path, sample_id=sample_id, out_h5ad=sample_h5ad, log=log)
        sample_entries.append({"sample_id": sample_id, "h5ad": str(sample_h5ad)})

    delegated_spec = {
        "input": {
            "samples": sample_entries,
            "qc": inp.get("qc"),
            "metadata": inp.get("metadata"),
        }
    }
    run_tenx_h5_multi(study=study, study_spec=delegated_spec, output_h5ad=output_h5ad, log=log)
    log_kv(log, "gse202379_raw_counts.complete", study=study, output=output_h5ad)
