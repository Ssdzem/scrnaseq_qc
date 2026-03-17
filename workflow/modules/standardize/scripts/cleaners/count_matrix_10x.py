from __future__ import annotations

from pathlib import Path
from typing import TextIO

import anndata as ad
import pandas as pd
import scanpy as sc


def _require_path(value: str | None, name: str) -> Path:
    if not value:
        raise SystemExit(f"Missing required input '{name}'")
    path = Path(value)
    if not path.exists():
        raise SystemExit(f"Input '{name}' not found: {path}")
    return path


def run_count_matrix_10x(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    inp = study_spec.get("input", {})
    matrix = _require_path(inp.get("matrix"), "matrix")
    genes = _require_path(inp.get("genes"), "genes")
    barcodes = _require_path(inp.get("barcodes"), "barcodes")
    metadata_path = inp.get("metadata")

    log.write(f"[{study}] Reading MTX triplet\n")

    # scanpy expects 10x-style files; we pass explicit filenames.
    adata = sc.read_mtx(str(matrix)).T

    genes_df = pd.read_csv(genes)
    barcodes_df = pd.read_csv(barcodes)

    # Try to be robust to common column naming.
    gene_col = "gene" if "gene" in genes_df.columns else genes_df.columns[0]
    cell_col = "cell" if "cell" in barcodes_df.columns else barcodes_df.columns[0]
    adata.var_names = genes_df[gene_col].astype(str).values
    adata.obs_names = barcodes_df[cell_col].astype(str).values

    if metadata_path:
        md = _require_path(metadata_path, "metadata")
        log.write(f"[{study}] Merging metadata: {md}\n")
        md_df = pd.read_csv(md)
        # Best-effort merge by barcode/cell.
        key = None
        for candidate in ["cell", "barcode", "Cell", "Barcode"]:
            if candidate in md_df.columns:
                key = candidate
                break
        if key is not None:
            md_df = md_df.set_index(key)
            md_df.index = md_df.index.astype(str)
            adata.obs = adata.obs.join(md_df, how="left")
        else:
            log.write(f"[{study}] WARNING: metadata has no obvious barcode column; skipping join\n")

    adata.uns["standardize"] = {"study": study, "process": "count_matrix_10x"}

    log.write(f"[{study}] Writing: {output_h5ad}\n")
    adata.write_h5ad(output_h5ad)

