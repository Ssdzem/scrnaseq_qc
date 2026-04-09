import argparse
import gzip
import json
import logging
import shutil
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd
import scanpy as sc
import yaml


DEFAULT_ZIP_URL = "https://www.livercellatlas.org/data_files/toDownload/rawData_human.zip"
DEFAULT_METADATA_URL = "https://www.livercellatlas.org/data_files/toDownload/annot_humanAll.csv"


def setup_logging(log_path: str | None):
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_path:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_path))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=handlers,
        force=True,
    )


def read_full_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def _download(url: str, target: Path):
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0:
        logging.info("Using cached file: %s", target)
        return
    logging.info("Downloading %s -> %s", url, target)
    with urllib.request.urlopen(url) as resp, target.open("wb") as out:
        shutil.copyfileobj(resp, out)


def _resolve_triplet_dir(extract_dir: Path, matrix_subdir: str) -> Path:
    candidate = extract_dir / matrix_subdir
    if candidate.exists():
        return candidate

    nested = list(extract_dir.glob(f"**/{matrix_subdir}"))
    if nested:
        return nested[0]

    raise FileNotFoundError(f"Could not find 10x matrix subdir '{matrix_subdir}' under {extract_dir}")


def _normalize_join_key(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.strip()
    # Common 10x mismatch: metadata may omit '-1' suffix present in barcodes.
    return s.str.replace(r"-1$", "", regex=True)


def _choose_var_names_mode(triplet_dir: Path) -> str:
    feature_candidates = (
        triplet_dir / "features.tsv.gz",
        triplet_dir / "features.tsv",
        triplet_dir / "genes.tsv.gz",
        triplet_dir / "genes.tsv",
    )
    feature_path = next((p for p in feature_candidates if p.exists()), None)
    if feature_path is None:
        logging.warning(
            "Could not locate features/genes table in %s. Defaulting var_names='gene_symbols'.",
            triplet_dir,
        )
        return "gene_symbols"

    preview = pd.read_csv(
        feature_path,
        sep="\t",
        header=None,
        nrows=10,
        compression="infer",
    )
    if preview.shape[1] >= 2:
        return "gene_symbols"

    logging.warning(
        "Detected single-column feature table (%s). "
        "Using var_names='gene_ids' for scanpy.read_10x_mtx.",
        feature_path,
    )
    return "gene_ids"


def _open_maybe_gzip(path: Path, mode: str):
    if ".gz" in path.suffixes:
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def _ensure_scanpy_compatible_features(triplet_dir: Path) -> None:
    """Patch single-column 10x feature tables to include expected columns.

    Some public datasets ship `features.tsv(.gz)` with only one column.
    Older/newer Scanpy builds still access column 1 internally, which raises
    `KeyError: 1` even when `var_names='gene_ids'`.
    """
    feature_candidates = (
        triplet_dir / "features.tsv.gz",
        triplet_dir / "features.tsv",
        triplet_dir / "genes.tsv.gz",
        triplet_dir / "genes.tsv",
    )
    feature_path = next((p for p in feature_candidates if p.exists()), None)
    if feature_path is None:
        return

    preview = pd.read_csv(
        feature_path,
        sep="\t",
        header=None,
        nrows=20,
        compression="infer",
    )
    if preview.shape[1] >= 2:
        return

    logging.warning(
        "Rewriting single-column feature table to Scanpy-compatible layout: %s",
        feature_path,
    )
    tmp_path = feature_path.with_suffix(feature_path.suffix + ".tmp")
    with _open_maybe_gzip(feature_path, "rt") as fin, _open_maybe_gzip(tmp_path, "wt") as fout:
        for raw in fin:
            gene = raw.rstrip("\n").split("\t")[0]
            if not gene:
                continue
            # 10x v3 compatible: gene_id, gene_symbol, feature_type
            fout.write(f"{gene}\t{gene}\tGene Expression\n")
    tmp_path.replace(feature_path)


def _prepare_symphony_reference(adata, full_cfg: dict):
    ann_cfg = full_cfg.get("annotation", {}) if isinstance(full_cfg, dict) else {}
    method = str(ann_cfg.get("method", "")).strip().lower()
    if method != "symphony":
        return adata

    import symphonypy as sp

    sym_cfg = ann_cfg.get("symphony", {}) if isinstance(ann_cfg, dict) else {}
    target_sum = float(sym_cfg.get("target_sum", 1e5))
    n_top_genes = int(sym_cfg.get("n_top_genes", 3000))
    n_comps = int(sym_cfg.get("n_comps", 30))
    batch_key_ref = str(sym_cfg.get("reference_batch_key", "")).strip()

    logging.info(
        "Preparing reference atlas for symphonypy (target_sum=%s, n_top_genes=%s, n_comps=%s, batch_key=%s)",
        target_sum,
        n_top_genes,
        n_comps,
        batch_key_ref or "<none>",
    )

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    hvg_kwargs = {"n_top_genes": min(max(100, n_top_genes), adata.n_vars)}
    if batch_key_ref and batch_key_ref in adata.obs.columns:
        hvg_kwargs["batch_key"] = batch_key_ref
    sc.pp.highly_variable_genes(adata, **hvg_kwargs)
    adata.raw = adata.copy()
    if "highly_variable" in adata.var.columns and bool(adata.var["highly_variable"].any()):
        adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    n_comps_eff = max(2, min(n_comps, adata.n_vars - 1)) if adata.n_vars > 2 else adata.n_vars
    sc.pp.pca(adata, n_comps=n_comps_eff, zero_center=False)

    if batch_key_ref and batch_key_ref in adata.obs.columns and adata.obs[batch_key_ref].nunique(dropna=True) > 1:
        sp.pp.harmony_integrate(adata, key=batch_key_ref)
    else:
        # Keep a Symphony-compatible slot even without Harmony correction.
        adata.obsm["X_pca_harmony"] = adata.obsm["X_pca"].copy()

    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    return adata


def build_reference_anndata(full_cfg: dict, output_h5ad: str) -> tuple[object, dict]:
    ref_cfg = full_cfg.get("reference_atlas", {})
    if not bool(ref_cfg.get("enabled", False)):
        raise ValueError("reference_atlas.enabled must be true to build the atlas reference")

    output_path = Path(output_h5ad)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cache_dir = Path(str(ref_cfg.get("cache_dir", "data/ref/cache")))
    extract_dir = Path(str(ref_cfg.get("extract_dir", "data/ref/rawData_human")))

    zip_url = str(ref_cfg.get("zip_url", DEFAULT_ZIP_URL)).strip()
    metadata_url = str(ref_cfg.get("metadata_url", DEFAULT_METADATA_URL)).strip()
    matrix_subdir = str(ref_cfg.get("matrix_subdir", "countTable_human")).strip()
    metadata_join_key = str(ref_cfg.get("metadata_join_key", "cell")).strip()
    metadata_label_col = str(ref_cfg.get("metadata_label_col", "annot")).strip()

    zip_path = cache_dir / Path(zip_url).name
    metadata_path = cache_dir / Path(metadata_url).name

    _download(zip_url, zip_path)
    _download(metadata_url, metadata_path)

    extract_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_dir)

    triplet_dir = _resolve_triplet_dir(extract_dir, matrix_subdir)
    logging.info("Reading 10x matrix from %s", triplet_dir)
    _ensure_scanpy_compatible_features(triplet_dir)
    var_names_mode = _choose_var_names_mode(triplet_dir)
    adata = sc.read_10x_mtx(str(triplet_dir), var_names=var_names_mode, make_unique=True)

    md = pd.read_csv(metadata_path)
    if metadata_join_key not in md.columns:
        raise ValueError(f"metadata_join_key not found in metadata file: {metadata_join_key}")
    if metadata_label_col not in md.columns:
        raise ValueError(f"metadata_label_col not found in metadata file: {metadata_label_col}")

    md = md.copy()
    md[metadata_join_key] = md[metadata_join_key].astype(str)
    md = md.drop_duplicates(subset=[metadata_join_key], keep="first").set_index(metadata_join_key)

    adata.obs["__barcode_raw"] = adata.obs_names.astype(str)
    adata.obs["__barcode_norm"] = _normalize_join_key(adata.obs["__barcode_raw"])
    md.index = _normalize_join_key(md.index.to_series())
    md = md[~md.index.duplicated(keep="first")]

    adata.obs = adata.obs.join(md, on="__barcode_norm", how="left")
    adata = _prepare_symphony_reference(adata, full_cfg)

    adata.uns.setdefault("integration_annotation", {})
    adata.uns["integration_annotation"]["reference_atlas"] = {
        "zip_url": zip_url,
        "metadata_url": metadata_url,
        "matrix_subdir": matrix_subdir,
        "metadata_join_key": metadata_join_key,
        "metadata_label_col": metadata_label_col,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_cells_with_label": int(adata.obs[metadata_label_col].notna().sum()),
    }

    adata.write_h5ad(str(output_path), compression="gzip")
    logging.info("Wrote reference atlas h5ad: %s", output_path)

    report = {
        "output_h5ad": str(output_path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "label_column": metadata_label_col,
        "n_cells_with_label": int(adata.obs[metadata_label_col].notna().sum()),
    }
    return adata, report


def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--log", default=None)
    args = p.parse_args(argv)

    setup_logging(args.log)
    cfg = read_full_config(args.config)
    _, report = build_reference_anndata(cfg, args.output)

    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    logging.info("Wrote reference build report: %s", report_path)


if __name__ == "__main__":
    main()
