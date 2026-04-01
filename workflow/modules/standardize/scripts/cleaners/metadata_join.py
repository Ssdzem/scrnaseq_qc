from __future__ import annotations

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


def run_metadata_join(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    inp = study_spec.get("input", {})
    h5ad_path = _require_path(inp.get("h5ad"), "h5ad")
    metadata_path = _require_path(inp.get("metadata"), "metadata")
    join_column = str(inp.get("join_column") or "").strip()
    if not join_column:
        raise SystemExit("Missing required input 'join_column' for metadata_join cleaner")

    log.write(f"[{study}] Reading h5ad: {h5ad_path}\n")
    log.flush()
    adata = sc.read_h5ad(str(h5ad_path))

    log.write(f"[{study}] Reading metadata: {metadata_path}\n")
    log.flush()
    md_df = pd.read_csv(metadata_path)
    if join_column not in md_df.columns:
        unnamed_cols = [col for col in md_df.columns if str(col).strip().lower().startswith("unnamed:")]
        if unnamed_cols:
            md_df = md_df.rename(columns={unnamed_cols[0]: join_column})
        else:
            md_with_index = pd.read_csv(metadata_path, index_col=0)
            md_with_index = md_with_index.copy()
            md_with_index[join_column] = md_with_index.index.astype(str)
            md_df = md_with_index.reset_index(drop=True)

    md_df = md_df.copy()
    md_df[join_column] = md_df[join_column].astype(str)
    md_df = md_df.set_index(join_column)

    obs_key = inp.get("obs_key")
    left = adata.obs.copy()
    if obs_key:
        obs_key = str(obs_key)
        if obs_key not in left.columns:
            raise SystemExit(f"obs_key '{obs_key}' not found in adata.obs columns")
    else:
        obs_key = join_column
        if obs_key not in left.columns:
            left[obs_key] = adata.obs_names.astype(str)

    left[obs_key] = left[obs_key].astype(str)
    adata.obs = left.merge(md_df, left_on=obs_key, right_index=True, how="left", sort=False)
    adata.obs.index = adata.obs_names

    adata.uns.setdefault("standardize", {})
    adata.uns["standardize"]["metadata_join"] = {"column": join_column, "obs_key": obs_key or "obs_names"}

    log.write(f"[{study}] Writing: {output_h5ad}\n")
    log.flush()
    adata.write_h5ad(output_h5ad)
