import argparse
import logging
from pathlib import Path

import numpy as np
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


def read_qc_config(path: str) -> dict:
    with open(path, "r") as f:
        cfg = yaml.safe_load(f) or {}
    return cfg.get("qc", {})


def ensure_obs_key(adata, key: str, fallback_value: str):
    if key not in adata.obs.columns:
        adata.obs[key] = fallback_value


def prepare_before_integration(adata, study: str, cfg: dict):
    batch_key = str(cfg.get("hvg_batch_key", "sample_id"))
    scvi_batch_key = str(cfg.get("scvi_batch_key", "sample_id"))
    n_top_genes = int(cfg.get("hvg_n_top_genes", 3000))
    seed = int(cfg.get("seed", 0))
    np.random.seed(seed)
    logging.info(
        "Pre-integration config: batch_key=%s scvi_batch_key=%s n_top_genes=%d seed=%d",
        batch_key,
        scvi_batch_key,
        n_top_genes,
        seed,
    )

    ensure_obs_key(adata, "study", study)
    ensure_obs_key(adata, "sample_id", study)
    ensure_obs_key(adata, batch_key, study)
    ensure_obs_key(adata, scvi_batch_key, study)

    adata.layers["counts"] = adata.X.copy()
    logging.info("Created counts layer")

    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.layers["normalized"] = adata.X.copy()
    logging.info("Completed normalize_total and stored normalized layer")

    sc.pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()
    logging.info("Completed log1p and stored log1p layer")

    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            batch_key=batch_key,
            flavor="seurat_v3",
            layer="counts",
        )
        logging.info("Computed highly variable genes with flavor=seurat_v3")
    except (ImportError, ValueError) as exc:
        logging.warning(
            "HVG with flavor=seurat_v3 failed (%s). Falling back to flavor=seurat.",
            exc,
        )
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            batch_key=batch_key,
            flavor="seurat",
        )
        logging.info("Computed highly variable genes with fallback flavor=seurat")

    adata.uns.setdefault("qc_preproc", {})
    adata.uns["qc_preproc"]["pre_integration_config"] = {
        "batch_key": batch_key,
        "scvi_batch_key": scvi_batch_key,
        "n_top_genes": n_top_genes,
        "seed": seed,
    }

    adata.uns["qc_preproc"]["hvg_var_names"] = adata.var_names[adata.var["highly_variable"]].tolist()
    logging.info(
        "Prepared data before integration: n_cells=%d n_vars=%d n_hvg=%d",
        adata.n_obs,
        adata.n_vars,
        len(adata.uns["qc_preproc"]["hvg_var_names"]),
    )
    return adata


def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["pre", "noop"], default="pre")
    p.add_argument("--study", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--qc-config", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--log", default=None)
    args = p.parse_args(argv)

    setup_logging(args.log)
    logging.info("explore_scanpy args: %s", vars(args))
    cfg = read_qc_config(args.qc_config)
    logging.info("Loaded QC config from %s: %s", args.qc_config, cfg)

    logging.info("Reading h5ad: %s", args.input)
    adata = sc.read_h5ad(args.input)

    if args.mode == "pre":
        logging.info("Running pre-integration preprocessing for: %s", args.study)
        adata = prepare_before_integration(adata, study=args.study, cfg=cfg)
    else:
        adata.uns.setdefault("qc_preproc", {})
        adata.uns["qc_preproc"]["pre_integration_skipped"] = True
        logging.info("Skipping pre-integration for this run (mode=noop)")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Writing output: %s", args.output)
    adata.write_h5ad(str(out), compression="gzip")
    logging.info("Done writing output: %s", args.output)


if __name__ == "__main__":
    main()
