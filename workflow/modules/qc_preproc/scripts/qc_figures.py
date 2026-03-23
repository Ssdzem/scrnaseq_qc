import argparse
import logging
from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
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
    )


def read_qc_config(path: str) -> dict:
    with open(path, "r") as f:
        cfg = yaml.safe_load(f) or {}
    return cfg.get("qc", {})


def add_qc_metrics(adata: ad.AnnData, mito_prefix: str, mito_prefixes_extra: list[str]):
    prefixes = [mito_prefix] + list(mito_prefixes_extra or [])
    mito_mask = np.zeros(adata.n_vars, dtype=bool)
    for prefix in prefixes:
        mito_mask |= adata.var_names.str.upper().str.startswith(str(prefix).upper())
    adata.var["mt"] = mito_mask
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


def _save_hist(series: pd.Series, outpath: str, title: str, xlabel: str, bins: int = 60):
    out = Path(outpath)
    out.parent.mkdir(parents=True, exist_ok=True)
    values = pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    plt.figure(figsize=(6, 4))
    plt.hist(values.values, bins=bins, color="#4C72B0", alpha=0.9)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Cells")
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()


def main(argv: list[str] | None = None):
    matplotlib.use("Agg")

    p = argparse.ArgumentParser()
    p.add_argument("--study", required=True)
    p.add_argument("--input-h5ad", required=True)
    p.add_argument("--qc-config", required=True)
    p.add_argument("--output-metrics", required=True)
    p.add_argument("--output-genes-fig", required=True)
    p.add_argument("--output-counts-fig", required=True)
    p.add_argument("--output-mito-fig", required=True)
    p.add_argument("--log", default=None)
    args = p.parse_args(argv)

    setup_logging(args.log)
    qc = read_qc_config(args.qc_config)
    if not bool(qc.get("make_figures", True)):
        logging.info("qc.make_figures is false; writing metrics only")

    mito_prefix = str(qc.get("mito_prefix", "MT-"))
    mito_prefixes_extra = qc.get("mito_prefixes_extra", [])

    logging.info("Reading h5ad: %s", args.input_h5ad)
    adata = sc.read_h5ad(args.input_h5ad)

    add_qc_metrics(adata, mito_prefix=mito_prefix, mito_prefixes_extra=mito_prefixes_extra)

    obs = adata.obs
    metrics = pd.DataFrame(
        {
            "study": args.study,
            "n_cells": [adata.n_obs],
            "n_genes": [adata.n_vars],
            "median_genes_per_cell": [float(np.median(obs.get("n_genes_by_counts", pd.Series(dtype=float))))],
            "median_counts_per_cell": [float(np.median(obs.get("total_counts", pd.Series(dtype=float))))],
            "median_mito_pct": [float(np.median(obs.get("pct_counts_mt", pd.Series(dtype=float))))],
        }
    )
    out_metrics = Path(args.output_metrics)
    out_metrics.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(out_metrics, index=False)

    if not bool(qc.get("make_figures", True)):
        for path in (args.output_genes_fig, args.output_counts_fig, args.output_mito_fig):
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            Path(path).touch()
        return

    _save_hist(
        obs.get("n_genes_by_counts", pd.Series(dtype=float)),
        args.output_genes_fig,
        title=f"{args.study}: genes per cell (pre-filter)",
        xlabel="n_genes_by_counts",
    )
    _save_hist(
        obs.get("total_counts", pd.Series(dtype=float)),
        args.output_counts_fig,
        title=f"{args.study}: counts per cell (pre-filter)",
        xlabel="total_counts",
    )
    _save_hist(
        obs.get("pct_counts_mt", pd.Series(dtype=float)),
        args.output_mito_fig,
        title=f"{args.study}: mito % (pre-filter)",
        xlabel="pct_counts_mt",
    )


if __name__ == "__main__":
    main()
