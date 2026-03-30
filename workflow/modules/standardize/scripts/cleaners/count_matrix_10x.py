from __future__ import annotations

import gzip
import shutil
from collections.abc import Callable
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

    if src.resolve() == dst.resolve():
        return

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

    if src.resolve() == dst.resolve():
        return

    if src.suffix == ".gz" and _is_gzip(src):
        shutil.copyfile(src, dst)
        return
    with src.open("rb") as fin, gzip.open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def _extract_barcodes_series(df: pd.DataFrame) -> pd.Series:
    if df.shape[1] == 1:
        return df.iloc[:, 0]
    key = None
    for candidate in ["barcode", "barcodes", "cell", "cells", "Cell", "Barcode"]:
        if candidate in df.columns:
            key = candidate
            break
    if key is None:
        key = df.columns[0]
    return df[key]


def _read_mtx_shape(path: Path) -> tuple[int, int]:
    """Return (n_genes, n_barcodes) from a MatrixMarket file path."""

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fin:
        for line in fin:
            s = line.strip()
            if not s or s.startswith("%"):
                continue
            parts = s.split()
            if len(parts) != 3:
                raise SystemExit(f"Invalid MatrixMarket shape line in {path}: {s!r}")
            n_genes, n_barcodes, _nnz = (int(parts[0]), int(parts[1]), int(parts[2]))
            return n_genes, n_barcodes
    raise SystemExit(f"Could not parse MatrixMarket shape from {path}")


def _select_by_expected_count(
    *,
    src: Path,
    label: str,
    expected_count: int,
    candidate_a: Callable[[], pd.Series],
    candidate_b: Callable[[], pd.Series],
) -> pd.Series:
    """Pick parsing mode that matches expected row count.

    This protects against CSV files that may or may not have headers.
    """

    series_a = candidate_a().astype(str)
    series_b = candidate_b().astype(str)

    len_a = len(series_a)
    len_b = len(series_b)
    a_matches = len_a == expected_count
    b_matches = len_b == expected_count

    if a_matches and not b_matches:
        return series_a
    if b_matches and not a_matches:
        return series_b
    if a_matches and b_matches:
        return series_a

    raise SystemExit(
        f"{label} count mismatch for {src}: expected {expected_count}, "
        f"got {len_a} (header=infer) and {len_b} (header=None)."
    )


def _normalize_genes_to_tsv_gz(src: Path, dst: Path, *, expected_count: int | None = None) -> None:
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

    def _extract(df: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
        if df.shape[1] == 1:
            ids = df.iloc[:, 0].astype(str)
            syms = df.iloc[:, 0].astype(str)
            return ids, syms

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
        return df[id_col].astype(str), df[symbol_col].astype(str)

    df_infer = pd.read_csv(src)
    ids_infer, syms_infer = _extract(df_infer)

    if expected_count is None:
        gene_ids, gene_symbols = ids_infer, syms_infer
    else:
        df_no_header = pd.read_csv(src, header=None)
        ids_no_header, syms_no_header = _extract(df_no_header)
        gene_ids = _select_by_expected_count(
            src=src,
            label="genes",
            expected_count=expected_count,
            candidate_a=lambda: ids_infer,
            candidate_b=lambda: ids_no_header,
        )
        gene_symbols = _select_by_expected_count(
            src=src,
            label="genes(symbols)",
            expected_count=expected_count,
            candidate_a=lambda: syms_infer,
            candidate_b=lambda: syms_no_header,
        )

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
    expected_genes, expected_barcodes = _read_mtx_shape(out_matrix)
    _normalize_genes_to_tsv_gz(genes, out_genes, expected_count=expected_genes)

    if barcodes.name.startswith("barcodes.tsv"):
        _ensure_gz_tsv(barcodes, out_barcodes)
    elif barcodes.suffix == ".csv":
        df_infer = pd.read_csv(barcodes)
        df_no_header = pd.read_csv(barcodes, header=None)
        barcodes_series = _select_by_expected_count(
            src=barcodes,
            label="barcodes",
            expected_count=expected_barcodes,
            candidate_a=lambda: _extract_barcodes_series(df_infer),
            candidate_b=lambda: _extract_barcodes_series(df_no_header),
        )
        with gzip.open(out_barcodes, "wt", encoding="utf-8") as fout:
            for val in barcodes_series.astype(str):
                fout.write(f"{val}\n")
    else:
        raise SystemExit(
            "Unsupported barcodes format. Expected barcodes.tsv(.gz) or a CSV. "
            f"Got: {barcodes}"
        )

    log.write(f"Normalized 10x triplet in-place under: {out_dir}\n")
    log.write(f"- {out_matrix.name}\n")
    log.write(f"- {out_genes.name}\n")
    log.write(f"- {out_barcodes.name}\n")


def _write_features_from_genes(genes_tsv_gz: Path, features_tsv_gz: Path) -> None:
    """Create a valid 10x-v3 features.tsv.gz (3 columns) from genes.tsv.gz."""

    with gzip.open(genes_tsv_gz, "rt", encoding="utf-8") as fin, gzip.open(
        features_tsv_gz, "wt", encoding="utf-8"
    ) as fout:
        for line in fin:
            row = line.rstrip("\n")
            if not row:
                continue
            parts = row.split("\t")
            if len(parts) == 1:
                gid = parts[0]
                gsym = parts[0]
            else:
                gid, gsym = parts[0], parts[1]
            fout.write(f"{gid}\t{gsym}\tGene Expression\n")


def _write_genes_from_features(features_tsv_gz: Path, genes_tsv_gz: Path) -> None:
    """Create a 10x-v2 genes.tsv.gz (2 columns) from features.tsv.gz."""

    with gzip.open(features_tsv_gz, "rt", encoding="utf-8") as fin, gzip.open(
        genes_tsv_gz, "wt", encoding="utf-8"
    ) as fout:
        for line in fin:
            row = line.rstrip("\n")
            if not row:
                continue
            parts = row.split("\t")
            if len(parts) == 1:
                gid = parts[0]
                gsym = parts[0]
            else:
                gid, gsym = parts[0], parts[1]
            fout.write(f"{gid}\t{gsym}\n")


def run_count_matrix_10x(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    inp = study_spec.get("input", {})
    matrix = _require_path(inp.get("matrix"), "matrix")
    genes = _require_path(inp.get("genes"), "genes")
    barcodes = _require_path(inp.get("barcodes"), "barcodes")
    metadata_path = inp.get("metadata")

    log.write(f"[{study}] Reading MTX triplet\n")
    log.flush()

    # Normalize in-place so downstream always sees a standard 10x directory.
    # This also avoids Scanpy's sensitivity to mixed compression and prefixes.
    _normalize_10x_triplet_inplace(matrix=matrix, genes=genes, barcodes=barcodes, log=log)

    resolved_dir = matrix.parent.resolve()
    log.write(f"[{study}] 10x resolved dir: {resolved_dir}\n")
    log.write(f"[{study}] 10x dir exists={resolved_dir.exists()} is_dir={resolved_dir.is_dir()}\n")
    try:
        entries = sorted(p.name for p in resolved_dir.iterdir())
    except Exception as exc:  # noqa: BLE001
        log.write(f"[{study}] ERROR listing 10x dir: {exc!r}\n")
        entries = []
    if entries:
        log.write(f"[{study}] 10x dir entries ({len(entries)}):\n")
        for name in entries:
            log.write(f"[{study}] - {name}\n")
    else:
        log.write(f"[{study}] 10x dir entries: <empty or unreadable>\n")
    log.flush()

    # Scanpy can look for either `genes.tsv.gz` (v2) or `features.tsv.gz` (v3)
    # depending on version/args. Ensure both exist to avoid FileNotFoundError.
    genes_tsv = resolved_dir / "genes.tsv.gz"
    features_tsv = resolved_dir / "features.tsv.gz"
    if genes_tsv.exists():
        try:
            _write_features_from_genes(genes_tsv, features_tsv)
            log.write(f"[{study}] Wrote 3-column features.tsv.gz from genes.tsv.gz\n")
        except Exception as exc:  # noqa: BLE001
            log.write(f"[{study}] WARNING: failed creating features.tsv.gz: {exc!r}\n")
    elif features_tsv.exists() and not genes_tsv.exists():
        try:
            _write_genes_from_features(features_tsv, genes_tsv)
            log.write(f"[{study}] Created 2-column genes.tsv.gz from features.tsv.gz\n")
        except Exception as exc:  # noqa: BLE001
            log.write(f"[{study}] WARNING: failed creating genes.tsv.gz: {exc!r}\n")
    log.flush()

    log.write(f"[{study}] scanpy.read_10x_mtx(dir={resolved_dir})\n")
    log.flush()
    adata = sc.read_10x_mtx(
        str(resolved_dir),
        var_names="gene_symbols",
        make_unique=True,
        # Keep default `compressed=True` and `gex_only=True` per Scanpy docs.
    )

    log.write(f"[{study}] Loaded matrix: n_obs={adata.n_obs} n_vars={adata.n_vars}\n")
    log.flush()

    if metadata_path:
        md = _require_path(metadata_path, "metadata")
        log.write(f"[{study}] Merging metadata: {md}\n")
        log.flush()
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
        log.flush()

    adata.uns["standardize"] = {"study": study, "process": "count_matrix_10x"}

    log.write(f"[{study}] Writing: {output_h5ad}\n")
    log.flush()
    adata.write_h5ad(output_h5ad)
