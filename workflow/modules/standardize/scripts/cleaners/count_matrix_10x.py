from __future__ import annotations

import gzip
import shutil
from pathlib import Path
from typing import TextIO

import pandas as pd
import scanpy as sc


def _require_path(value: str | None, name: str) -> Path:
    if not value:
        raise SystemExit(f"Missing required input '{name}'")
    path = Path(value)
    if not path.exists():
        raise SystemExit(f"Input '{name}' not found: {path}")
    return path


def _ensure_gz_tsv(src: Path, dst: Path) -> None:
    """Write `src` into `dst` as gzipped text.

    If `src` is already gzipped, copy as-is.
    """

    if src.suffix == ".gz":
        shutil.copyfile(src, dst)
        return
    with src.open("rb") as fin, gzip.open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as f:
        return f.read(2) == b"\x1f\x8b"


def _ensure_gz_bytes(src: Path, dst: Path) -> None:
    """Ensure `dst` is gzipped bytes representing `src` content.

    Some files may end with `.gz` but not actually be gzipped.
    """

    if src.suffix == ".gz" and _is_gzip(src):
        shutil.copyfile(src, dst)
        return
    with src.open("rb") as fin, gzip.open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def _normalize_barcodes_to_tsv_gz(src: Path, dst: Path) -> None:
    """Write a 10x-compatible `barcodes.tsv.gz`.

    Supports:
    - existing `barcodes.tsv(.gz)` (copied/gzipped)
    - CSV with a single column of barcodes
    - CSV with a barcode column named like `barcode`/`cell`
    """

    if src.name.startswith("barcodes.tsv"):
        _ensure_gz_tsv(src, dst)
        return

    if src.suffix != ".csv":
        raise SystemExit(
            "Unsupported barcodes format. Expected barcodes.tsv(.gz) or a CSV. "
            f"Got: {src}"
        )

    df = pd.read_csv(src)
    if df.shape[1] == 1:
        series = df.iloc[:, 0]
    else:
        key = None
        for candidate in ["barcode", "barcodes", "cell", "cells", "Cell", "Barcode"]:
            if candidate in df.columns:
                key = candidate
                break
        if key is None:
            raise SystemExit(
                f"Barcode CSV {src} must have a single column or a column named like barcode/cell. "
                f"Columns: {list(df.columns)}"
            )
        series = df[key]

    barcodes = series.astype(str)
    with gzip.open(dst, "wt", encoding="utf-8") as fout:
        for val in barcodes:
            fout.write(f"{val}\n")


def _normalize_genes_to_tsv_gz(src: Path, dst: Path) -> None:
    """Write a 10x v2-style `genes.tsv.gz`.

    Supports:
    - existing `genes.tsv(.gz)`/`features.tsv(.gz)` (copied/gzipped to genes.tsv.gz)
    - CSV with one column of gene symbols
    - CSV with gene symbol/id columns (best effort)
    """

    if src.name.startswith("genes.tsv") or src.name.startswith("features.tsv"):
        _ensure_gz_tsv(src, dst)
        return

    if src.suffix != ".csv":
        raise SystemExit(
            "Unsupported genes format. Expected genes.tsv(.gz)/features.tsv(.gz) or a CSV. "
            f"Got: {src}"
        )

    df = pd.read_csv(src)
    if df.shape[1] == 1:
        gene_ids = df.iloc[:, 0].astype(str)
        gene_symbols = df.iloc[:, 0].astype(str)
    else:
        id_col = None
        symbol_col = None
        for candidate in ["gene_id", "gene_ids", "ensembl_id", "GeneID", "id", "ID"]:
            if candidate in df.columns:
                id_col = candidate
                break
        for candidate in ["gene", "gene_name", "symbol", "gene_symbol", "Gene", "GeneName", "GeneSymbol"]:
            if candidate in df.columns:
                symbol_col = candidate
                break
        if symbol_col is None:
            symbol_col = df.columns[0]
        if id_col is None:
            id_col = symbol_col
        gene_ids = df[id_col].astype(str)
        gene_symbols = df[symbol_col].astype(str)

    with gzip.open(dst, "wt", encoding="utf-8") as fout:
        for gid, gsym in zip(gene_ids, gene_symbols, strict=False):
            fout.write(f"{gid}\t{gsym}\n")


def _normalize_10x_triplet_inplace(*, matrix: Path, genes: Path, barcodes: Path, log: TextIO) -> None:
    """Normalize input files to a Scanpy-friendly 10x directory.

    Writes the following files next to the provided inputs:
    - matrix.mtx.gz
    - genes.tsv.gz
    - barcodes.tsv.gz

    Also converts common non-10x inputs (CSV) into 10x TSVs.
    """

    if matrix.parent != genes.parent or matrix.parent != barcodes.parent:
        raise SystemExit(
            "10x inputs must be in the same directory. "
            f"Got matrix={matrix.parent}, genes={genes.parent}, barcodes={barcodes.parent}"
        )

    out_dir = matrix.parent
    out_matrix = out_dir / "matrix.mtx.gz"
    out_genes = out_dir / "genes.tsv.gz"
    out_barcodes = out_dir / "barcodes.tsv.gz"

    # Matrix must be a MatrixMarket file.
    if ".mtx" not in matrix.name:
        raise SystemExit(f"Expected MatrixMarket input matrix.mtx(.gz). Got: {matrix}")

    # Make compression consistent: always gz, and normalize naming.
    _ensure_gz_bytes(matrix, out_matrix)
    _normalize_genes_to_tsv_gz(genes, out_genes)
    _normalize_barcodes_to_tsv_gz(barcodes, out_barcodes)

    log.write(f"Normalized 10x triplet in-place under: {out_dir}\n")
    log.write(f"- {out_matrix.name}\n")
    log.write(f"- {out_genes.name}\n")
    log.write(f"- {out_barcodes.name}\n")


def run_count_matrix_10x(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    inp = study_spec.get("input", {})
    matrix = _require_path(inp.get("matrix"), "matrix")
    genes = _require_path(inp.get("genes"), "genes")
    barcodes = _require_path(inp.get("barcodes"), "barcodes")
    metadata_path = inp.get("metadata")

    log.write(f"[{study}] Reading MTX triplet\n")

    # Normalize in-place so downstream always sees a standard 10x directory.
    # This also avoids Scanpy's sensitivity to mixed compression and prefixes.
    _normalize_10x_triplet_inplace(matrix=matrix, genes=genes, barcodes=barcodes, log=log)
    adata = sc.read_10x_mtx(
        str(matrix.parent),
        var_names="gene_symbols",
        make_unique=True,
        # Keep default `compressed=True` and `gex_only=True` per Scanpy docs.
    )

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
