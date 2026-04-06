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


def annotate(adata, cfg: dict):
    method = str(cfg.get("method", "from_existing_column"))
    source_col = str(cfg.get("source_column", "cell_type"))
    target_col = str(cfg.get("target_column", "cell_type"))
    unknown_label = str(cfg.get("unknown_label", "Unknown"))

    if method == "from_existing_column" and source_col in adata.obs.columns:
        adata.obs[target_col] = adata.obs[source_col].astype(str).fillna(unknown_label)
    elif method == "from_leiden":
        if "leiden" in adata.obs.columns:
            adata.obs[target_col] = "Cluster_" + adata.obs["leiden"].astype(str)
        else:
            adata.obs[target_col] = unknown_label
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
