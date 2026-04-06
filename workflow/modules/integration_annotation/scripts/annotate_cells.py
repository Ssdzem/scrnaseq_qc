import argparse
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

    adata_query.obs[target_col] = _safe_obs_str(adata_query.obs[target_col], unknown_label)

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

    adata_query.uns.setdefault("integration_annotation", {})
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
