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
    except ImportError as exc:
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


def run_integration(adata, cfg: dict, do_integration: bool = True):
    integration_method = str(cfg.get("integration_method", "scvi")).strip().lower()
    n_pcs = int(cfg.get("n_pcs", 50))
    n_neighbors = int(cfg.get("n_neighbors", 15))
    leiden_resolution = float(cfg.get("leiden_resolution", 1.0))
    max_epochs = int(cfg.get("scvi_max_epochs", 100))
    seed = int(cfg.get("seed", 0))
    np.random.seed(seed)

    pre_cfg = adata.uns.get("qc_preproc", {}).get("pre_integration_config", {})
    scvi_batch_key = str(pre_cfg.get("scvi_batch_key", cfg.get("scvi_batch_key", "sample_id")))
    hvg_var_names = list(adata.uns.get("qc_preproc", {}).get("hvg_var_names", []))
    if not hvg_var_names:
        raise ValueError("Missing qc_preproc.hvg_var_names; run pre-integration mode first.")

    adata_hvg = adata[:, hvg_var_names].copy()
    adata_hvg.layers["log1p"] = adata_hvg.X.copy()
    logging.info("Created HVG subset for integration: n_cells=%d n_hvg=%d", adata_hvg.n_obs, adata_hvg.n_vars)

    sc.pp.scale(adata_hvg, max_value=10)
    adata_hvg.layers["scaled"] = adata_hvg.X.copy()
    logging.info("Scaled HVG data")

    sc.tl.pca(adata_hvg, n_comps=n_pcs, svd_solver="arpack")
    logging.info("Computed PCA")

    if do_integration and integration_method == "scvi":
        import scvi

        scvi.settings.seed = seed
        scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts", batch_key=scvi_batch_key)
        model = scvi.model.SCVI(adata_hvg)
        model.train(
          max_epochs=400,
          validation_size=0.1,
          check_val_every_n_epoch=1,
          early_stopping=True,
          early_stopping_monitor="elbo_validation",
          early_stopping_patience=20,
        )        
        adata_hvg.obsm["X_scVI"] = model.get_latent_representation()
        integrated_rep = "X_scVI"
        logging.info("Trained scVI and saved latent representation")
    elif do_integration and integration_method == "harmony":
        import harmonypy as hm

        ho = hm.run_harmony(
            adata_hvg.obsm["X_pca"],
            adata_hvg.obs,
            [scvi_batch_key],
        )
        z_corr = ho.Z_corr
        if hasattr(z_corr, "detach"):
            z_corr = z_corr.detach().cpu().numpy()
        z_corr = np.asarray(z_corr)

        if z_corr.ndim != 2:
            raise ValueError(
                f"Harmony returned unexpected ndim={z_corr.ndim} with shape={z_corr.shape}."
            )

        if z_corr.shape == (adata_hvg.n_obs, n_pcs):
            pass
        elif z_corr.shape == (n_pcs, adata_hvg.n_obs):
            z_corr = z_corr.T
        else:
            raise ValueError(
                "Harmony returned unexpected shape="
                f"{z_corr.shape}; expected {(adata_hvg.n_obs, n_pcs)} or {(n_pcs, adata_hvg.n_obs)}."
            )

        adata_hvg.obsm["X_harmony"] = z_corr
        integrated_rep = "X_harmony"
        logging.info("Ran Harmony integration and saved X_harmony embedding with shape=%s", z_corr.shape)
    elif do_integration:
        raise ValueError(f"Unsupported integration_method={integration_method!r}; expected 'scvi' or 'harmony'.")
    else:
        integrated_rep = "X_pca"
        adata_hvg.obsm["X_scVI"] = adata_hvg.obsm["X_pca"]
        integration_method = "none"
        logging.info("Skipping integration backend; using PCA representation directly")

    sc.pp.neighbors(adata_hvg, use_rep=integrated_rep, n_neighbors=n_neighbors)
    sc.tl.umap(adata_hvg)
    sc.tl.leiden(adata_hvg, key_added="leiden_scvi", resolution=leiden_resolution)
    logging.info("Computed neighbors, UMAP, and Leiden clusters")

    adata.uns.setdefault("qc_preproc", {})
    adata.uns["qc_preproc"]["integration_config"] = {
        "integration_method": integration_method,
        "integration_batch_key": scvi_batch_key,
        "n_pcs": n_pcs,
        "n_neighbors": n_neighbors,
        "leiden_resolution": leiden_resolution,
        "scvi_max_epochs": max_epochs,
        "seed": seed,
    }

    adata.uns["qc_preproc"]["hvg_var_names"] = adata_hvg.var_names.tolist()
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.obsm["X_scVI"] = adata_hvg.obsm["X_scVI"]
    if "X_harmony" in adata_hvg.obsm:
        adata.obsm["X_harmony"] = adata_hvg.obsm["X_harmony"]
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.obs["leiden_scvi"] = adata_hvg.obs["leiden_scvi"].astype(str).reindex(adata.obs_names)
    logging.info("Transferred embeddings and clusters back to full AnnData")

    return adata


def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["full", "pre", "integration", "no_integration", "noop"], default="full")
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

    if args.mode in {"full", "pre"}:
        logging.info("Running pre-integration preprocessing for: %s", args.study)
        adata = prepare_before_integration(adata, study=args.study, cfg=cfg)
    if args.mode in {"full", "integration", "no_integration"}:
        if args.mode == "no_integration":
            logging.info("Running downstream embedding/clustering without integration for: %s", args.study)
        else:
            logging.info("Running integration step for: %s", args.study)
        adata = run_integration(adata, cfg=cfg, do_integration=(args.mode != "no_integration"))
    if args.mode == "noop":
        adata.uns.setdefault("qc_preproc", {})
        adata.uns["qc_preproc"]["study_explore_integration_skipped"] = True
        logging.info("Skipping integration for this run (mode=noop)")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Writing output: %s", args.output)
    adata.write_h5ad(str(out), compression="gzip")
    logging.info("Done writing output: %s", args.output)


if __name__ == "__main__":
    main()
