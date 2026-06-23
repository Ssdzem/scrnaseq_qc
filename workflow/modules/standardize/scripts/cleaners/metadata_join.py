from __future__ import annotations

from pathlib import Path
from typing import TextIO

import pandas as pd
import scanpy as sc


CANONICAL_ID_COLUMNS: tuple[str, ...] = ("sample_id", "study_accession", "donor_id", "run_id")


def _require_path(value: str | None, name: str) -> Path:
    if not value:
        raise SystemExit(f"Missing required input '{name}'")
    path = Path(value)
    if not path.exists():
        raise SystemExit(f"Input '{name}' not found: {path}")
    return path


def _find_metadata_catalog(*, inp: dict, study: str, log: TextIO) -> Path | None:
    explicit = str(inp.get("metadata_catalog") or "").strip()
    if explicit:
        p = Path(explicit)
        if p.exists():
            log.write(f"[{study}] metadata catalog (explicit): {p}\n")
            return p
        log.write(f"[{study}] WARNING: explicit metadata_catalog not found: {p}\n")

    candidates = [
        Path("config/liver_scRNAseq_metadata_everything.csv"),
        Path("a1_p1_qc/config/liver_scRNAseq_metadata_everything.csv"),
    ]
    for p in candidates:
        if p.exists():
            log.write(f"[{study}] metadata catalog (auto): {p}\n")
            return p

    log.write(f"[{study}] metadata catalog not found; skipping value-based canonicalization\n")
    return None


def _first_n_unique_nonempty(values: pd.Series, *, n: int) -> list[str]:
    out: list[str] = []
    for v in values.astype(str):
        s = v.strip()
        if not s:
            continue
        if s not in out:
            out.append(s)
        if len(out) >= n:
            break
    return out


def _match_score(*, probes: list[str], candidate_values: pd.Series) -> tuple[int, int]:
    if not probes:
        return (0, 0)
    probe_set = set(probes)
    candidate_set = {str(v).strip() for v in candidate_values.astype(str) if str(v).strip()}
    hits = len(probe_set & candidate_set)
    return (hits, len(probe_set))


def _value_based_column_mapping(
    *,
    study: str,
    md_df: pd.DataFrame,
    inp: dict,
    log: TextIO,
) -> pd.DataFrame:
    catalog_path = _find_metadata_catalog(inp=inp, study=study, log=log)
    if catalog_path is None:
        return md_df

    try:
        catalog = pd.read_csv(catalog_path, dtype=str)
    except Exception as exc:  # noqa: BLE001
        log.write(f"[{study}] WARNING: failed reading metadata catalog {catalog_path}: {exc}\n")
        return md_df

    if "study_accession" not in catalog.columns:
        log.write(f"[{study}] WARNING: metadata catalog missing study_accession; skipping canonicalization\n")
        return md_df

    study_rows = catalog.loc[catalog["study_accession"].astype(str).str.strip() == study].copy()
    if study_rows.empty:
        log.write(f"[{study}] metadata catalog has no rows for study; skipping canonicalization\n")
        return md_df

    result = md_df.copy()
    already_used: set[str] = set()
    for canonical in CANONICAL_ID_COLUMNS:
        if canonical in result.columns:
            log.write(f"[{study}] canonicalization: {canonical} already present; keeping existing column\n")
            already_used.add(canonical)
            continue
        if canonical not in study_rows.columns:
            log.write(f"[{study}] canonicalization: {canonical} not in metadata catalog columns\n")
            continue

        probes = _first_n_unique_nonempty(study_rows[canonical], n=3)
        log.write(f"[{study}] canonicalization probes for {canonical}: {probes}\n")
        if not probes:
            continue

        best_col: str | None = None
        best_hits = -1
        best_total = len(probes)
        for col in result.columns:
            if col in already_used:
                continue
            hits, total = _match_score(probes=probes, candidate_values=result[col])
            log.write(f"[{study}] canonicalization score {canonical} <- {col}: {hits}/{total}\n")
            if hits > best_hits:
                best_hits = hits
                best_col = col

        if best_col is None or best_hits <= 0:
            log.write(f"[{study}] canonicalization: no value match found for {canonical}\n")
            continue

        if best_hits < min(2, best_total):
            log.write(
                f"[{study}] canonicalization: weak match for {canonical} from {best_col} ({best_hits}/{best_total}); skipping\n"
            )
            continue

        result = result.rename(columns={best_col: canonical})
        already_used.add(canonical)
        log.write(f"[{study}] canonicalization: mapped {best_col} -> {canonical} ({best_hits}/{best_total})\n")

    return result


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
    md_df = _value_based_column_mapping(study=study, md_df=md_df, inp=inp, log=log)
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
