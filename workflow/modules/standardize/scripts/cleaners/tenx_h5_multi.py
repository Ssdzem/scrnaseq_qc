from __future__ import annotations

from pathlib import Path
from typing import TextIO

import anndata as ad
import pandas as pd
import scanpy as sc


def _require_list(value, name: str) -> list:
    if not value or not isinstance(value, list):
        raise SystemExit(f"Missing required list '{name}'")
    return value


def _require_path(value: str | None, name: str) -> Path:
    if not value:
        raise SystemExit(f"Missing required input '{name}'")
    path = Path(value)
    if not path.exists():
        raise SystemExit(f"Input '{name}' not found: {path}")
    return path


def _normalize_sample_table(samples: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(samples)
    if df.empty:
        raise SystemExit("'samples' list is empty")
    if "sample_id" not in df.columns:
        raise SystemExit("Each sample entry must contain 'sample_id'")
    if "h5" not in df.columns:
        raise SystemExit("Each sample entry must contain 'h5'")
    if df["sample_id"].isna().any():
        raise SystemExit("Some samples are missing 'sample_id'")
    if df["h5"].isna().any():
        raise SystemExit("Some samples are missing 'h5'")
    if df["sample_id"].duplicated().any():
        dups = sorted(df.loc[df["sample_id"].duplicated(), "sample_id"].astype(str).unique())
        raise SystemExit(f"Duplicate sample_id values in spec: {dups}")
    return df


def run_tenx_h5_multi(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    inp = study_spec.get("input", {})
    samples = _require_list(inp.get("samples"), "input.samples")
    samples_df = _normalize_sample_table(samples)

    adatas: list[ad.AnnData] = []
    for row in samples_df.itertuples(index=False):
        sample_id = str(row.sample_id)
        h5_path = _require_path(str(row.h5), f"samples[{sample_id}].h5")
        log.write(f"[{study}] Reading 10x h5: {h5_path} (sample_id={sample_id})\n")
        log.flush()
        a = sc.read_10x_h5(str(h5_path))
        a.var_names_make_unique()
        a.obs["sample_id"] = sample_id
        a.obs["study_accession"] = study
        log.write(f"[{study}] Loaded {sample_id}: n_obs={a.n_obs} n_vars={a.n_vars}\n")
        log.flush()
        adatas.append(a)

    log.write(f"[{study}] Concatenating {len(adatas)} samples\n")
    log.flush()
    merged = ad.concat(adatas, axis=0, join="outer", label=None, merge="same", index_unique="-")
    log.write(f"[{study}] Concatenated: n_obs={merged.n_obs} n_vars={merged.n_vars}\n")
    log.flush()

    metadata_path = inp.get("metadata")
    if metadata_path:
        md = _require_path(metadata_path, "metadata")
        log.write(f"[{study}] Merging metadata: {md}\n")
        log.flush()
        md_df = pd.read_csv(md)
        if "sample_id" in md_df.columns:
            md_df["sample_id"] = md_df["sample_id"].astype(str)
        key = None
        for candidate in ["cell", "barcode", "Cell", "Barcode"]:
            if candidate in md_df.columns:
                key = candidate
                break
        if key is not None:
            md_df = md_df.set_index(key)
            md_df.index = md_df.index.astype(str)
            merged.obs = merged.obs.join(md_df, how="left")
        else:
            log.write(f"[{study}] WARNING: metadata has no obvious barcode column; skipping join\n")
        log.flush()

    merged.uns["standardize"] = {"study": study, "process": "tenx_h5_multi"}
    log.write(f"[{study}] Writing: {output_h5ad}\n")
    log.flush()
    merged.write_h5ad(output_h5ad)
