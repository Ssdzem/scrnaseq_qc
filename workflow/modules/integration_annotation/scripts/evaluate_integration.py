import argparse
import json
import logging
from pathlib import Path

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


def compute_metrics(adata, cfg: dict) -> dict:
    batch_key = str(cfg.get("batch_key", "sample_id"))
    label_key = str(cfg.get("label_key", "cell_type"))
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

    return metrics


def flatten_metrics(metrics: dict) -> pd.DataFrame:
    flat = {}
    for k, v in metrics.items():
        if isinstance(v, (dict, list)):
            flat[k] = json.dumps(v, sort_keys=True)
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
    out_json.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    flatten_metrics(metrics).to_csv(out_csv, index=False)

    logging.info("Wrote evaluation metrics: %s and %s", out_json, out_csv)


if __name__ == "__main__":
    main()
