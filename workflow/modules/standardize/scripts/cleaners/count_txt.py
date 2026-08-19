from __future__ import annotations

from pathlib import Path
from typing import TextIO

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from rtlog import log_kv, log_line

# Reuse the shared count-table -> 10x-h5 helpers so both handlers stay in sync.
from cleaners.csv_handler import _infer_needs_transpose, _write_10x_h5_from_adata


def _require_path(value: str | None, name: str) -> Path:
    if not value:
        raise SystemExit(f"Missing required input '{name}'")
    path = Path(value)
    if not path.exists():
        raise SystemExit(f"Input '{name}' not found: {path}")
    return path


def _sample_id_from_name(filename: str, suffix: str) -> str:
    if suffix and filename.endswith(suffix):
        return filename[: -len(suffix)]
    return Path(filename).stem


def _pick_count_tables(*, raw_dir: Path, pattern: str, suffix: str) -> list[tuple[str, Path]]:
    """Collect one count table per sample, recursing per the glob pattern.

    Sample ids come from the file name with `suffix` stripped. Files can live in
    sub-directories (e.g. normal/ and pre_RCT/); a sample id collision across
    directories is an error so the caller disambiguates the filenames instead of
    silently dropping a sample.
    """

    grouped: dict[str, Path] = {}
    for p in sorted(raw_dir.glob(pattern)):
        if not p.is_file():
            continue
        sample_id = _sample_id_from_name(p.name, suffix)
        if sample_id in grouped:
            raise SystemExit(
                f"Duplicate sample_id '{sample_id}' from {p} and {grouped[sample_id]}. "
                "Disambiguate the filenames (or narrow input.pattern)."
            )
        grouped[sample_id] = p
    return [(sample_id, grouped[sample_id]) for sample_id in sorted(grouped)]


def _needs_transpose(df: pd.DataFrame, *, genes_in: str, sample_id: str) -> bool:
    """Decide whether rows are genes (True) or cells (False).

    `genes_in` pins the orientation explicitly ("rows"/"columns"); "auto" falls
    back to the shared shape/label heuristic, which assumes the larger axis is
    genes and only inspects labels on an exact tie. Pin it whenever the layout is
    known (e.g. genes-as-rows count tables) so a sample with more cells than
    genes cannot silently flip the axes.
    """

    orient = (genes_in or "auto").strip().lower()
    if orient == "rows":
        return True
    if orient in ("columns", "cols"):
        return False
    if orient == "auto":
        return _infer_needs_transpose(df)
    raise SystemExit(
        f"[{sample_id}] input.genes_in must be one of rows|columns|auto, got {genes_in!r}"
    )


def _table_to_adata(
    *, table_path: Path, sample_id: str, sep: str, genes_in: str, log: TextIO
) -> ad.AnnData:
    log_line(log, f"[{sample_id}] Reading count table: {table_path}")

    read_kwargs: dict[str, object] = {"index_col": 0}
    if sep in ("whitespace", r"\s+", "\\s+"):
        read_kwargs["sep"] = r"\s+"
        read_kwargs["engine"] = "python"
    else:
        read_kwargs["sep"] = sep
    df = pd.read_csv(table_path, **read_kwargs)
    if df.empty:
        raise SystemExit(f"Count table is empty: {table_path}")

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    if _needs_transpose(df, genes_in=genes_in, sample_id=sample_id):
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


def run_count_txt(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    """Standardize per-sample plain-text count tables (e.g. `*.count.txt`).

    Mirrors `csv_handler`: each table is converted to an intermediate 10x `.h5`
    and registered under `input.samples`, so a `count_txt -> tenx_h5_multi`
    pipeline picks them up for the shared QC/merge path.
    """

    inp = study_spec.setdefault("input", {})
    raw_dir = _require_path(inp.get("raw_dir"), "raw_dir")
    pattern = str(inp.get("pattern", "**/*.count.txt"))
    suffix = str(inp.get("suffix", ".count.txt"))
    sep = str(inp.get("sep", "\t"))
    genes_in = str(inp.get("genes_in", "auto"))

    tmp_h5_dir = output_h5ad.parent / ".tmp" / f"{study}_count_txt_h5"
    tmp_h5_dir.mkdir(parents=True, exist_ok=True)

    tables = _pick_count_tables(raw_dir=raw_dir, pattern=pattern, suffix=suffix)
    if not tables:
        raise SystemExit(f"No files matched '{pattern}' under {raw_dir}")

    log_kv(log, "count_txt.start", study=study, raw_dir=raw_dir, samples=len(tables))

    sample_entries: list[dict[str, str]] = []
    for sample_id, table_path in tables:
        adata = _table_to_adata(
            table_path=table_path, sample_id=sample_id, sep=sep, genes_in=genes_in, log=log
        )
        sample_h5 = tmp_h5_dir / f"{sample_id}.h5"
        _write_10x_h5_from_adata(adata=adata, out_h5=sample_h5, sample_id=sample_id, log=log)
        sample_entries.append({"sample_id": sample_id, "h5": str(sample_h5)})

    inp["samples"] = sample_entries

    placeholder = ad.AnnData(
        X=sparse.csr_matrix((0, 0), dtype=np.float32),
        obs=pd.DataFrame(index=pd.Index([], name="barcode")),
        var=pd.DataFrame(index=pd.Index([], name="gene_symbol")),
    )
    placeholder.uns["standardize"] = {"study": study, "process": "count_txt"}
    placeholder.write_h5ad(output_h5ad)
    log_kv(log, "count_txt.complete", study=study, output=output_h5ad, prepared_samples=len(sample_entries))
