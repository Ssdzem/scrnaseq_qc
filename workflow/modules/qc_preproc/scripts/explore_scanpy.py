import argparse
import logging
from pathlib import Path

import numpy as np
import scanpy as sc
import scvi
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


def ensure_obs_key(adata, key: str, fallback_value: str):
    if key not in adata.obs.columns:
        adata.obs[key] = fallback_value


def run_exploration(adata, study: str, cfg: dict):
    batch_key = str(cfg.get("hvg_batch_key", "sample_id"))
    scvi_batch_key = str(cfg.get("scvi_batch_key", "sample_id"))
    n_top_genes = int(cfg.get("hvg_n_top_genes", 3000))
    n_pcs = int(cfg.get("n_pcs", 50))
    n_neighbors = int(cfg.get("n_neighbors", 15))
    leiden_resolution = float(cfg.get("leiden_resolution", 1.0))
    max_epochs = int(cfg.get("scvi_max_epochs", 100))
    seed = int(cfg.get("seed", 0))
    np.random.seed(seed)
    scvi.settings.seed = seed

    ensure_obs_key(adata, "study", study)
    ensure_obs_key(adata, "sample_id", study)
    ensure_obs_key(adata, batch_key, study)
    ensure_obs_key(adata, scvi_batch_key, study)

    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.layers["normalized"] = adata.X.copy()

    sc.pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        batch_key=batch_key,
        flavor="seurat_v3",
        layer="counts",
    )

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    adata_hvg.layers["log1p"] = adata_hvg.X.copy()

    sc.pp.scale(adata_hvg, max_value=10)
    adata_hvg.layers["scaled"] = adata_hvg.X.copy()

    sc.tl.pca(adata_hvg, n_comps=n_pcs, svd_solver="arpack")

    scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts", batch_key=scvi_batch_key)
    model = scvi.model.SCVI(adata_hvg)
    model.train(max_epochs=max_epochs)
    adata_hvg.obsm["X_scVI"] = model.get_latent_representation()

    sc.pp.neighbors(adata_hvg, use_rep="X_scVI", n_neighbors=n_neighbors)
    sc.tl.umap(adata_hvg)
    sc.tl.leiden(adata_hvg, key_added="leiden_scvi", resolution=leiden_resolution)

    adata.uns.setdefault("qc_preproc", {})
    adata.uns["qc_preproc"]["explore_config"] = {
        "batch_key": batch_key,
        "scvi_batch_key": scvi_batch_key,
        "n_top_genes": n_top_genes,
        "n_pcs": n_pcs,
        "n_neighbors": n_neighbors,
        "leiden_resolution": leiden_resolution,
        "scvi_max_epochs": max_epochs,
        "seed": seed,
    }

    adata.uns["qc_preproc"]["hvg_var_names"] = adata_hvg.var_names.tolist()
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.obsm["X_scVI"] = adata_hvg.obsm["X_scVI"]
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.obs["leiden_scvi"] = adata_hvg.obs["leiden_scvi"].astype(str).reindex(adata.obs_names)

    return adata


def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser()
    p.add_argument("--study", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--qc-config", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--log", default=None)
    args = p.parse_args(argv)

    setup_logging(args.log)
    cfg = read_qc_config(args.qc_config)

    logging.info("Reading h5ad: %s", args.input)
    adata = sc.read_h5ad(args.input)

    logging.info("Running exploratory preprocessing for: %s", args.study)
    adata = run_exploration(adata, study=args.study, cfg=cfg)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Writing output: %s", args.output)
    adata.write_h5ad(str(out), compression="gzip")


if __name__ == "__main__":
    main()
