import argparse
import json
import logging
import os
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


def resolve_eval_cfg(full_cfg: dict) -> dict:
    return full_cfg.get("evaluation", {})


def _to_float_or_none(value):
    try:
        if value is None:
            return None
        if isinstance(value, (bool, str)):
            return None
        out = float(value)
        if np.isnan(out) or np.isinf(out):
            return None
        return out
    except Exception:
        return None


def _series_summary(series: pd.Series, direction: str) -> dict:
    vals = pd.to_numeric(series, errors="coerce").dropna()
    if vals.empty:
        return {"available": False}
    return {
        "available": True,
        "n_valid": int(vals.size),
        "mean": float(vals.mean()),
        "median": float(vals.median()),
        "p90": float(vals.quantile(0.9)),
        "p95": float(vals.quantile(0.95)),
        "min": float(vals.min()),
        "max": float(vals.max()),
        "interpretation": direction,
    }


def _coerce_scib_output_to_dict(result) -> dict:
    if isinstance(result, pd.DataFrame):
        if result.empty:
            return {}
        row = result.iloc[0].to_dict()
    elif isinstance(result, pd.Series):
        row = result.to_dict()
    elif isinstance(result, dict):
        row = dict(result)
    else:
        return {"raw": str(result)}
    out = {}
    for k, v in row.items():
        fv = _to_float_or_none(v)
        out[str(k)] = fv if fv is not None else v
    return out


def _compute_scib_for_embedding(adata, batch_key: str, label_key: str, embed_key: str, cfg: dict) -> dict:
    import scib

    metrics_kwargs = {
        "adata": adata,
        "adata_int": adata,
        "batch_key": batch_key,
        "label_key": label_key,
        "embed": embed_key,
        "cluster_key": f"cluster_{embed_key}",
        "ari_": bool(cfg.get("ari", True)),
        "nmi_": bool(cfg.get("nmi", True)),
        "silhouette_": bool(cfg.get("silhouette", True)),
        "pcr_": bool(cfg.get("pcr", True)),
        "cell_cycle_": bool(cfg.get("cell_cycle", True)),
        "hvg_score_": bool(cfg.get("hvg_overlap", True)),
        "isolated_labels_": bool(cfg.get("isolated_labels", True)),
        "isolated_labels_f1_": bool(cfg.get("isolated_labels_f1", True)),
        "isolated_labels_asw_": bool(cfg.get("isolated_labels_asw", True)),
        "graph_conn_": bool(cfg.get("graph_connectivity", True)),
        "trajectory_": bool(cfg.get("trajectory", True)),
        "kBET_": bool(cfg.get("kbet", True)),
        "lisi_graph_": bool(cfg.get("lisi_graph", True)),
        "ilisi_": bool(cfg.get("ilisi", True)),
        "clisi_": bool(cfg.get("clisi", True)),
        "subsample": float(cfg.get("lisi_subsample", 0.5)),
        "n_cores": int(os.environ.get("SCRNASEQ_QC_IA_EVAL_THREADS", cfg.get("n_cores", 1))),
        "type_": str(cfg.get("type", "embed")),
        "verbose": bool(cfg.get("verbose", False)),
    }
    result = scib.metrics.metrics(**metrics_kwargs)
    return _coerce_scib_output_to_dict(result)


def compute_scib_metrics(adata, cfg: dict, batch_key: str, label_key: str) -> dict:
    out = {
        "enabled": bool(cfg.get("enabled", True)),
        "results_by_embedding": {},
        "errors_by_embedding": {},
        "available_embeddings": sorted(list(adata.obsm.keys())),
    }
    if not out["enabled"]:
        return out

    if batch_key not in adata.obs.columns:
        out["error"] = f"batch_key missing in obs: {batch_key}"
        return out
    if label_key not in adata.obs.columns:
        out["error"] = f"label_key missing in obs: {label_key}"
        return out

    embed_candidates = cfg.get(
        "embeddings",
        ["X_integrated", "X_scVI", "X_scanvi", "X_symphony", "X_harmony", "X_pca_harmony", "X_pca"],
    )
    embed_keys = [str(x) for x in embed_candidates if str(x) in adata.obsm]
    out["evaluated_embeddings"] = embed_keys
    if not embed_keys:
        out["error"] = "no requested embeddings were found in adata.obsm"
        return out

    for embed in embed_keys:
        try:
            out["results_by_embedding"][embed] = _compute_scib_for_embedding(adata, batch_key, label_key, embed, cfg)
        except Exception as exc:
            out["errors_by_embedding"][embed] = str(exc)
    return out


def _resolve_label_key(adata, cfg: dict) -> str:
    configured = str(cfg.get("label_key", "")).strip()
    if configured:
        return configured
    ann_meta = adata.uns.get("integration_annotation", {}).get("annotation", {})
    target_col = str(ann_meta.get("target_column", "")).strip()
    if target_col:
        return target_col
    return "cell_type"


def compute_metrics(adata, cfg: dict) -> dict:
    batch_key = str(cfg.get("batch_key", "sample_id")).strip()
    label_key = _resolve_label_key(adata, cfg)
    unknown_label = str(cfg.get("unknown_label", "Unknown"))

    metrics: dict[str, object] = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "has_X_integrated": bool("X_integrated" in adata.obsm),
        "has_X_umap": bool("X_umap" in adata.obsm),
        "has_leiden": bool("leiden" in adata.obs.columns),
        "batch_key": batch_key,
        "label_key": label_key,
        "n_batches": int(adata.obs[batch_key].nunique()) if batch_key in adata.obs.columns else 0,
        "n_labels": int(adata.obs[label_key].nunique()) if label_key in adata.obs.columns else 0,
    }

    if label_key in adata.obs.columns:
        labels = adata.obs[label_key].astype(str)
        metrics["unknown_fraction"] = float(labels.eq(unknown_label).mean())
    else:
        metrics["unknown_fraction"] = 1.0

    if "leiden" in adata.obs.columns:
        metrics["n_leiden_clusters"] = int(adata.obs["leiden"].astype(str).nunique())
    else:
        metrics["n_leiden_clusters"] = 0

    if "integration_annotation" in adata.uns:
        metrics["integration_annotation"] = adata.uns["integration_annotation"]

    mapping = {}
    if "symphony_per_cell_dist" in adata.obs.columns:
        mapping["symphony_per_cell_dist"] = _series_summary(adata.obs["symphony_per_cell_dist"], "lower_is_better")
    if "scanvi_pred_confidence" in adata.obs.columns:
        mapping["scanvi_pred_confidence"] = _series_summary(adata.obs["scanvi_pred_confidence"], "higher_is_better")
    if "scanvi_pred_entropy" in adata.obs.columns:
        mapping["scanvi_pred_entropy"] = _series_summary(adata.obs["scanvi_pred_entropy"], "lower_is_better")
    if "scanvi_pred_margin" in adata.obs.columns:
        mapping["scanvi_pred_margin"] = _series_summary(adata.obs["scanvi_pred_margin"], "higher_is_better")
    metrics["mapping_scores"] = mapping

    metrics["scib"] = compute_scib_metrics(adata, cfg.get("scib", {}), batch_key=batch_key, label_key=label_key)

    return metrics


def _json_safe(value):
    """Let json.dumps write the numpy values that scanpy and scib return."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def flatten_metrics(metrics: dict) -> pd.DataFrame:
    flat = {}
    for k, v in metrics.items():
        if isinstance(v, (dict, list)):
            flat[k] = json.dumps(v, sort_keys=True, default=_json_safe)
        else:
            flat[k] = v
    return pd.DataFrame([flat])


def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser()
    p.add_argument("--scope", choices=["all_studies", "study"], required=True)
    p.add_argument("--study", default="")
    p.add_argument("--input", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--output-json", required=True)
    p.add_argument("--output-csv", required=True)
    p.add_argument("--log", default=None)
    args = p.parse_args(argv)

    setup_logging(args.log)
    cfg = read_full_config(args.config)
    eval_cfg = resolve_eval_cfg(cfg)

    adata = sc.read_h5ad(args.input)
    metrics = compute_metrics(adata, eval_cfg)

    out_json = Path(args.output_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(metrics, indent=2, default=_json_safe) + "\n", encoding="utf-8")

    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    flatten_metrics(metrics).to_csv(out_csv, index=False)

    logging.info("Wrote evaluation metrics: %s and %s", out_json, out_csv)


if __name__ == "__main__":
    main()
