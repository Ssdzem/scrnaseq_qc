import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml


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


def resolve_annotation_cfg(full_cfg: dict) -> dict:
    return full_cfg.get("annotation", {})


def _safe_obs_str(series: pd.Series, unknown_label: str) -> pd.Series:
    out = series.copy()
    out = out.where(~out.isna(), other=unknown_label)
    out = out.astype(str)
    out = out.replace({"nan": unknown_label, "None": unknown_label, "NA": unknown_label, "": unknown_label})
    return out


def _annotate_with_symphony(adata_query, cfg: dict, target_col: str, unknown_label: str):
    import symphonypy as sp

    reference_path = str(cfg.get("reference_h5ad", "")).strip()
    if not reference_path:
        raise ValueError("annotation.symphony.reference_h5ad must be set for method=symphony")

    reference_label_key = str(cfg.get("reference_label_key", "cell_type")).strip()
    reference_batch_key = str(cfg.get("reference_batch_key", "")).strip()
    query_batch_key = str(cfg.get("query_batch_key", "sample_id")).strip()
    symphony_k = int(cfg.get("k", 5))

    adata_ref = sc.read_h5ad(reference_path)
    if reference_label_key not in adata_ref.obs.columns:
        raise ValueError(f"Reference label key not found in atlas: {reference_label_key}")

    map_kwargs = {}
    if query_batch_key and query_batch_key in adata_query.obs.columns:
        map_kwargs["key"] = query_batch_key
    sp.tl.map_embedding(adata_query, adata_ref, **map_kwargs)
    if "X_pca_harmony" in adata_query.obsm:
        adata_query.obsm["X_symphony"] = adata_query.obsm["X_pca_harmony"]

    if reference_batch_key and reference_batch_key not in adata_ref.obs.columns:
        logging.warning(
            "annotation.symphony.reference_batch_key=%s is not present in reference obs; continuing without it",
            reference_batch_key,
        )

    sp.tl.transfer_labels_kNN(
        adata_query,
        adata_ref,
        reference_label_key,
        query_labels=target_col,
        k=symphony_k,
    )

    # Symphony mapping confidence metrics (lower distance means better map confidence).
    try:
        if hasattr(sp.tl, "per_cell_confidence"):
            sp.tl.per_cell_confidence(adata_query, adata_ref)
    except Exception as exc:
        logging.warning("Could not compute Symphony per-cell confidence: %s", exc)

    try:
        if hasattr(sp.tl, "per_cluster_confidence") and "leiden" in adata_query.obs.columns:
            sp.tl.per_cluster_confidence(adata_query, adata_ref, "leiden")
    except Exception as exc:
        logging.warning("Could not compute Symphony per-cluster confidence: %s", exc)

    adata_query.obs[target_col] = _safe_obs_str(adata_query.obs[target_col], unknown_label)
    adata_query.uns.setdefault("integration_annotation", {})
    per_cell_col = "symphony_per_cell_dist"
    if per_cell_col in adata_query.obs.columns:
        vals = pd.to_numeric(adata_query.obs[per_cell_col], errors="coerce").dropna()
        per_cell_summary = {
            "n_valid": int(vals.size),
            "mean": float(vals.mean()) if vals.size else None,
            "median": float(vals.median()) if vals.size else None,
            "p90": float(vals.quantile(0.90)) if vals.size else None,
            "p95": float(vals.quantile(0.95)) if vals.size else None,
            "min": float(vals.min()) if vals.size else None,
            "max": float(vals.max()) if vals.size else None,
            "interpretation": "lower_is_better",
        }
    else:
        per_cell_summary = {"available": False}

    adata_query.uns["integration_annotation"]["symphony"] = {
        "reference_h5ad": reference_path,
        "reference_label_key": reference_label_key,
        "reference_batch_key": reference_batch_key,
        "query_batch_key": query_batch_key,
        "k": symphony_k,
        "embedding_key": "X_symphony" if "X_symphony" in adata_query.obsm else "",
        "per_cell_score_key": per_cell_col if per_cell_col in adata_query.obs.columns else "",
        "per_cluster_score_key": "symphony_per_cluster_dist" if "symphony_per_cluster_dist" in adata_query.uns else "",
        "per_cell_summary": per_cell_summary,
    }

    return adata_query


def _annotate_with_scanvi(adata_query, cfg: dict, target_col: str, unknown_label: str):
    import scvi

    reference_path = str(cfg.get("reference_h5ad", "")).strip()
    if not reference_path:
        raise ValueError("annotation.scanvi.reference_h5ad must be set for method=scanvi")

    reference_label_key = str(cfg.get("reference_label_key", "cell_type")).strip()
    batch_key = str(cfg.get("batch_key", "sample_id")).strip()
    counts_layer = str(cfg.get("counts_layer", "counts")).strip()
    scvi_max_epochs = int(cfg.get("scvi_max_epochs", 100))
    scanvi_max_epochs = int(cfg.get("scanvi_max_epochs", 30))
    n_latent = int(cfg.get("n_latent", 30))
    seed = int(cfg.get("seed", 0))

    adata_ref = sc.read_h5ad(reference_path)
    if reference_label_key not in adata_ref.obs.columns:
        raise ValueError(f"Reference label key not found in atlas: {reference_label_key}")

    adata_ref = adata_ref.copy()
    adata_query = adata_query.copy()
    adata_ref.obs["__orig_obs_name"] = adata_ref.obs_names.astype(str)
    adata_query.obs["__orig_obs_name"] = adata_query.obs_names.astype(str)

    if batch_key and batch_key not in adata_ref.obs.columns:
        adata_ref.obs[batch_key] = "reference_batch"
    if batch_key and batch_key not in adata_query.obs.columns:
        adata_query.obs[batch_key] = "query_batch"

    train_label_key = "__scanvi_train_labels"
    source_key = "__scanvi_source"
    adata_ref.obs[train_label_key] = _safe_obs_str(adata_ref.obs[reference_label_key], unknown_label)
    adata_query.obs[train_label_key] = unknown_label

    adata_joint = sc.concat(
        [adata_ref, adata_query],
        label=source_key,
        keys=["reference", "query"],
        join="inner",
        index_unique="-",
        merge="same",
    )
    if adata_joint.n_vars == 0:
        raise ValueError("No overlapping genes between query and scanvi reference atlas")

    layer = counts_layer if counts_layer and counts_layer in adata_joint.layers else None
    scvi.settings.seed = seed
    scvi.model.SCVI.setup_anndata(adata_joint, batch_key=batch_key or None, labels_key=train_label_key, layer=layer)
    scvi_model = scvi.model.SCVI(adata_joint, n_latent=n_latent)
    scvi_model.train(max_epochs=scvi_max_epochs)

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key=train_label_key,
        unlabeled_category=unknown_label,
    )
    scanvi_model.train(max_epochs=scanvi_max_epochs)

    pred = scanvi_model.predict(adata_joint)
    pred_series = pd.Series(pred, index=adata_joint.obs_names)
    try:
        pred_soft = scanvi_model.predict(adata_joint, soft=True)
        if isinstance(pred_soft, pd.DataFrame):
            pred_soft_df = pred_soft.copy()
        else:
            pred_soft_df = pd.DataFrame(pred_soft, index=adata_joint.obs_names)
    except Exception as exc:
        logging.warning("Could not compute SCANVI soft predictions: %s", exc)
        pred_soft_df = pd.DataFrame(index=adata_joint.obs_names)

    try:
        latent_joint = scanvi_model.get_latent_representation(adata_joint)
    except Exception as exc:
        logging.warning("Could not extract SCANVI latent representation: %s", exc)
        latent_joint = None

    query_mask = adata_joint.obs[source_key].eq("query")
    pred_table = pd.DataFrame(
        {
            "orig_obs_name": adata_joint.obs.loc[query_mask, "__orig_obs_name"].astype(str).values,
            "pred_label": pred_series.loc[query_mask].astype(str).values,
        }
    )
    pred_map = pred_table.drop_duplicates(subset=["orig_obs_name"]).set_index("orig_obs_name")["pred_label"]

    adata_query.obs[target_col] = unknown_label
    adata_query.obs.loc[:, target_col] = adata_query.obs["__orig_obs_name"].map(pred_map).fillna(unknown_label).astype(str)
    adata_query.obs[target_col] = _safe_obs_str(adata_query.obs[target_col], unknown_label)

    if not pred_soft_df.empty:
        pred_soft_df = pred_soft_df.loc[query_mask]
        pred_soft_df.index = adata_joint.obs.loc[query_mask, "__orig_obs_name"].astype(str).values
        pred_soft_df = pred_soft_df.reindex(adata_query.obs["__orig_obs_name"].astype(str).values)
        probs = pred_soft_df.to_numpy(dtype=float, copy=True)
        probs = np.clip(probs, 1e-12, 1.0)
        row_sums = probs.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        probs = probs / row_sums
        conf = probs.max(axis=1)
        entropy = -(probs * np.log(probs)).sum(axis=1)
        if probs.shape[1] > 1:
            top2 = np.partition(probs, kth=-2, axis=1)[:, -2:]
            margin = top2[:, 1] - top2[:, 0]
        else:
            margin = np.ones(probs.shape[0], dtype=float)
        adata_query.obs["scanvi_pred_confidence"] = conf
        adata_query.obs["scanvi_pred_entropy"] = entropy
        adata_query.obs["scanvi_pred_margin"] = margin
        adata_query.obs["scanvi_pred_confidence"] = pd.to_numeric(adata_query.obs["scanvi_pred_confidence"], errors="coerce")
        adata_query.obs["scanvi_pred_entropy"] = pd.to_numeric(adata_query.obs["scanvi_pred_entropy"], errors="coerce")
        adata_query.obs["scanvi_pred_margin"] = pd.to_numeric(adata_query.obs["scanvi_pred_margin"], errors="coerce")
    if latent_joint is not None:
        latent_query = np.asarray(latent_joint)[query_mask.to_numpy(), :]
        adata_query.obsm["X_scanvi"] = latent_query

    adata_query.uns.setdefault("integration_annotation", {})
    conf_series = (
        pd.to_numeric(adata_query.obs["scanvi_pred_confidence"], errors="coerce")
        if "scanvi_pred_confidence" in adata_query.obs.columns
        else pd.Series(dtype=float)
    )
    entropy_series = (
        pd.to_numeric(adata_query.obs["scanvi_pred_entropy"], errors="coerce")
        if "scanvi_pred_entropy" in adata_query.obs.columns
        else pd.Series(dtype=float)
    )
    margin_series = (
        pd.to_numeric(adata_query.obs["scanvi_pred_margin"], errors="coerce")
        if "scanvi_pred_margin" in adata_query.obs.columns
        else pd.Series(dtype=float)
    )
    adata_query.uns["integration_annotation"]["scanvi"] = {
        "reference_h5ad": reference_path,
        "reference_label_key": reference_label_key,
        "batch_key": batch_key,
        "counts_layer": layer or "X",
        "scvi_max_epochs": scvi_max_epochs,
        "scanvi_max_epochs": scanvi_max_epochs,
        "n_latent": n_latent,
        "seed": seed,
        "n_reference_cells": int(adata_ref.n_obs),
        "n_query_cells": int(adata_query.n_obs),
        "n_shared_genes": int(adata_joint.n_vars),
        "embedding_key": "X_scanvi" if "X_scanvi" in adata_query.obsm else "",
        "mapping_score_analog": {
            "confidence_key": "scanvi_pred_confidence" if "scanvi_pred_confidence" in adata_query.obs.columns else "",
            "entropy_key": "scanvi_pred_entropy" if "scanvi_pred_entropy" in adata_query.obs.columns else "",
            "margin_key": "scanvi_pred_margin" if "scanvi_pred_margin" in adata_query.obs.columns else "",
            "confidence_mean": float(conf_series.mean()) if conf_series.notna().any() else None,
            "confidence_median": float(conf_series.median()) if conf_series.notna().any() else None,
            "entropy_mean": float(entropy_series.mean()) if entropy_series.notna().any() else None,
            "entropy_median": float(entropy_series.median()) if entropy_series.notna().any() else None,
            "margin_mean": float(margin_series.mean()) if margin_series.notna().any() else None,
            "margin_median": float(margin_series.median()) if margin_series.notna().any() else None,
            "interpretation": {
                "scanvi_pred_confidence": "higher_is_better",
                "scanvi_pred_entropy": "lower_is_better",
                "scanvi_pred_margin": "higher_is_better",
            },
        },
    }
    return adata_query


def annotate(adata, cfg: dict):
    method = str(cfg.get("method", "from_existing_column"))
    source_col = str(cfg.get("source_column", "cell_type"))
    target_col = str(cfg.get("target_column", "cell_type"))
    unknown_label = str(cfg.get("unknown_label", "Unknown"))

    if method == "from_existing_column" and source_col in adata.obs.columns:
        adata.obs[target_col] = _safe_obs_str(adata.obs[source_col], unknown_label)
    elif method == "from_leiden":
        if "leiden" in adata.obs.columns:
            adata.obs[target_col] = "Cluster_" + adata.obs["leiden"].astype(str)
        else:
            adata.obs[target_col] = unknown_label
    elif method == "symphony":
        adata = _annotate_with_symphony(adata, cfg.get("symphony", {}), target_col, unknown_label)
    elif method == "scanvi":
        adata = _annotate_with_scanvi(adata, cfg.get("scanvi", {}), target_col, unknown_label)
    else:
        adata.obs[target_col] = unknown_label

    adata.uns.setdefault("integration_annotation", {})
    adata.uns["integration_annotation"]["annotation"] = {
        "method": method,
        "source_column": source_col,
        "target_column": target_col,
        "unknown_label": unknown_label,
    }
    return adata, target_col, unknown_label


def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser()
    p.add_argument("--scope", choices=["all_studies", "study"], required=True)
    p.add_argument("--study", default="")
    p.add_argument("--input", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--log", default=None)
    args = p.parse_args(argv)

    setup_logging(args.log)
    cfg = read_full_config(args.config)
    ann_cfg = resolve_annotation_cfg(cfg)

    adata = sc.read_h5ad(args.input)
    adata, target_col, unknown_label = annotate(adata, ann_cfg)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(out), compression="gzip")
    logging.info("Wrote annotated h5ad: %s", out)

    counts = adata.obs[target_col].astype(str).value_counts(dropna=False).rename_axis("label").reset_index(name="n_cells")
    counts["is_unknown"] = counts["label"].eq(unknown_label)
    summary = Path(args.summary)
    summary.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(summary, index=False)
    logging.info("Wrote annotation summary: %s", summary)


if __name__ == "__main__":
    main()
