import argparse
import json
import logging
from pathlib import Path

import anndata as ad
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


def resolve_ia_cfg(full_cfg: dict) -> dict:
    return full_cfg.get("integration", {})


def ensure_obs_key(adata, key: str, fallback_value: str):
    if key not in adata.obs.columns:
        adata.obs[key] = fallback_value


def _as_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        value = [value]
    return [str(v).strip() for v in value if str(v).strip()]


def _coerce_device_value(value):
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.lower() == "auto":
            return "auto"
        try:
            return int(stripped)
        except ValueError:
            return stripped
    return value


def _require_obs_keys(adata, keys: list[str]):
    missing = [key for key in keys if key not in adata.obs.columns]
    if missing:
        raise ValueError(f"Required integration obs column(s) missing: {missing}")


def _sanitize_categorical_obs(adata, keys: list[str], unknown_label: str = "Unknown"):
    for key in keys:
        series = adata.obs[key]
        series = series.where(~series.isna(), other=unknown_label)
        series = series.astype(str).replace({"": unknown_label, "nan": unknown_label, "None": unknown_label})
        adata.obs[key] = pd.Categorical(series)


def run_integration(adata, cfg: dict, scope_label: str):
    integration_method = str(cfg.get("method", cfg.get("integration_method", "scvi"))).strip().lower()
    do_integration = bool(cfg.get("enabled", True))
    hvg_key = str(cfg.get("hvg_key", "qc_preproc.hvg_var_names"))
    n_pcs = int(cfg.get("n_pcs", 50))
    n_neighbors = int(cfg.get("n_neighbors", 15))
    leiden_resolution = float(cfg.get("leiden_resolution", 1.0))
    compute_pca = bool(cfg.get("compute_pca", integration_method != "scvi" or not do_integration))
    compute_umap = bool(cfg.get("compute_umap", False))
    compute_leiden = bool(cfg.get("compute_leiden", True))
    scvi_batch_key = str(cfg.get("batch_key", "sample_id"))
    categorical_covariate_keys = _as_list(cfg.get("categorical_covariate_keys", []))
    continuous_covariate_keys = _as_list(cfg.get("continuous_covariate_keys", []))
    batch_keys_for_reporting = [scvi_batch_key, *categorical_covariate_keys, *continuous_covariate_keys]
    scvi_max_epochs = int(cfg.get("scvi_max_epochs", 100))
    n_latent = int(cfg.get("n_latent", 30))
    batch_size = int(cfg.get("batch_size", 4096))
    validation_size = float(cfg.get("validation_size", 0.1))
    check_val_every_n_epoch = int(cfg.get("check_val_every_n_epoch", 1))
    early_stopping = bool(cfg.get("early_stopping", True))
    early_stopping_monitor = str(cfg.get("early_stopping_monitor", "elbo_validation"))
    early_stopping_patience = int(cfg.get("early_stopping_patience", 20))
    dl_num_workers = int(cfg.get("dl_num_workers", 4))
    num_threads = int(cfg.get("num_threads", 16))
    use_gpu = bool(cfg.get("use_gpu", False))
    require_gpu = bool(cfg.get("require_gpu", False))
    accelerator = str(cfg.get("accelerator", "gpu" if use_gpu else "cpu")).strip().lower()
    devices = _coerce_device_value(cfg.get("devices", 1 if accelerator != "auto" else "auto"))
    seed = int(cfg.get("seed", 0))
    np.random.seed(seed)

    ensure_obs_key(adata, "sample_id", scope_label)
    _require_obs_keys(adata, batch_keys_for_reporting)
    _sanitize_categorical_obs(adata, [scvi_batch_key, *categorical_covariate_keys])

    # Prefer the precomputed HVGs from qc_preproc; fallback to all vars.
    hvg_var_names = []
    if hvg_key == "qc_preproc.hvg_var_names":
        hvg_var_names = list(adata.uns.get("qc_preproc", {}).get("hvg_var_names", []))
    if not hvg_var_names:
        hvg_var_names = adata.var_names.tolist()

    adata_hvg = adata[:, hvg_var_names].copy()
    if "counts" in adata_hvg.layers:
        layer_for_scvi = "counts"
    else:
        layer_for_scvi = None
    if "log1p" not in adata_hvg.layers:
        adata_hvg.layers["log1p"] = adata_hvg.X.copy()

    if compute_pca or integration_method in {"harmony", "none"}:
        pca_source = adata_hvg.copy()
        sc.pp.scale(pca_source, max_value=10)
        sc.tl.pca(pca_source, n_comps=n_pcs, svd_solver="arpack")
        adata_hvg.obsm["X_pca"] = pca_source.obsm["X_pca"]
        adata_hvg.uns["pca"] = pca_source.uns["pca"]
        adata_hvg.varm["PCs"] = pca_source.varm["PCs"]

    if do_integration and integration_method == "scvi":
        import scvi
        import torch

        scvi.settings.seed = seed
        if hasattr(torch, "set_float32_matmul_precision"):
            torch.set_float32_matmul_precision(str(cfg.get("float32_matmul_precision", "high")))
        cuda_available = bool(torch.cuda.is_available())
        if require_gpu and accelerator == "gpu" and not cuda_available:
            raise RuntimeError("integration.require_gpu=true but torch.cuda.is_available() is false")

        scvi.model.SCVI.setup_anndata(
            adata_hvg,
            layer=layer_for_scvi,
            batch_key=scvi_batch_key,
            categorical_covariate_keys=categorical_covariate_keys or None,
            continuous_covariate_keys=continuous_covariate_keys or None,
        )
        model = scvi.model.SCVI(adata_hvg, n_latent=n_latent)
        scvi.settings.dl_num_workers = dl_num_workers
        scvi.settings.num_threads = num_threads
        scvi.settings.dl_persistent_workers = dl_num_workers > 0
        logging.info(
            "Training scVI: batch_key=%s categorical_covariates=%s layer=%s accelerator=%s devices=%s "
            "cuda_available=%s batch_size=%d max_epochs=%d",
            scvi_batch_key,
            categorical_covariate_keys,
            layer_for_scvi or "X",
            accelerator,
            devices,
            cuda_available,
            batch_size,
            scvi_max_epochs,
        )
        model.train(
            max_epochs=scvi_max_epochs,
            validation_size=validation_size,
            accelerator=accelerator,
            devices=devices,
            batch_size=batch_size,
            check_val_every_n_epoch=check_val_every_n_epoch,
            early_stopping=early_stopping,
            early_stopping_monitor=early_stopping_monitor,
            early_stopping_patience=early_stopping_patience,
        )
        adata_hvg.obsm["X_integrated"] = model.get_latent_representation()
        adata_hvg.obsm["X_scVI"] = adata_hvg.obsm["X_integrated"]
        cluster_key = "leiden_scvi"
    elif do_integration and integration_method == "harmony":
        import harmonypy as hm

        ho = hm.run_harmony(adata_hvg.obsm["X_pca"], adata_hvg.obs, [scvi_batch_key, *categorical_covariate_keys])
        z_corr = np.asarray(ho.Z_corr)
        if z_corr.shape == (n_pcs, adata_hvg.n_obs):
            z_corr = z_corr.T
        if z_corr.shape != (adata_hvg.n_obs, n_pcs):
            raise ValueError(f"Unexpected Harmony shape: {z_corr.shape}")
        adata_hvg.obsm["X_integrated"] = z_corr
        adata_hvg.obsm["X_harmony"] = z_corr
        cluster_key = "leiden_harmony"
    else:
        adata_hvg.obsm["X_integrated"] = adata_hvg.obsm["X_pca"]
        integration_method = "none"
        cluster_key = "leiden_pca"

    if compute_leiden or compute_umap:
        sc.pp.neighbors(adata_hvg, use_rep="X_integrated", n_neighbors=n_neighbors)
    if compute_umap:
        sc.tl.umap(adata_hvg)
    if compute_leiden:
        sc.tl.leiden(adata_hvg, key_added=cluster_key, resolution=leiden_resolution)

    if "X_pca" in adata_hvg.obsm:
        adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.obsm["X_integrated"] = adata_hvg.obsm["X_integrated"]
    if "X_scVI" in adata_hvg.obsm:
        adata.obsm["X_scVI"] = adata_hvg.obsm["X_scVI"]
    if "X_harmony" in adata_hvg.obsm:
        adata.obsm["X_harmony"] = adata_hvg.obsm["X_harmony"]
    if "X_umap" in adata_hvg.obsm:
        adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    if compute_leiden and cluster_key in adata_hvg.obs.columns:
        adata.obs["leiden"] = adata_hvg.obs[cluster_key].astype(str).reindex(adata.obs_names)
        adata.obs[cluster_key] = adata.obs["leiden"]

    adata.uns.setdefault("integration_annotation", {})
    adata.uns["integration_annotation"]["integration"] = {
        "method": integration_method,
        "enabled": do_integration,
        "n_pcs": n_pcs,
        "n_neighbors": n_neighbors,
        "leiden_resolution": leiden_resolution,
        "compute_pca": compute_pca,
        "compute_umap": compute_umap,
        "compute_leiden": compute_leiden,
        "batch_key": scvi_batch_key,
        "categorical_covariate_keys": categorical_covariate_keys,
        "continuous_covariate_keys": continuous_covariate_keys,
        "scvi_max_epochs": scvi_max_epochs,
        "n_latent": n_latent,
        "batch_size": batch_size,
        "use_gpu": use_gpu,
        "require_gpu": require_gpu,
        "accelerator": accelerator,
        "devices": devices,
        "layer": layer_for_scvi or "X",
        "seed": seed,
        "scope": scope_label,
        "n_hvg_used": len(hvg_var_names),
    }

    return adata, adata.uns["integration_annotation"]["integration"]


def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser()
    p.add_argument("--scope", choices=["all_studies", "study"], required=True)
    p.add_argument("--study", default="")
    p.add_argument("--input", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--log", default=None)
    args = p.parse_args(argv)

    setup_logging(args.log)
    logging.info("Runtime: scanpy=%s anndata=%s", sc.__version__, ad.__version__)
    logging.info("Runtime paths: scanpy=%s anndata=%s", sc.__file__, ad.__file__)
    cfg = read_full_config(args.config)
    ia_cfg = resolve_ia_cfg(cfg)

    scope_label = args.study if args.scope == "study" and args.study else "all_studies"

    logging.info("Reading input: %s", args.input)
    adata = sc.read_h5ad(args.input)
    adata, report = run_integration(adata, ia_cfg, scope_label=scope_label)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(out), compression="gzip")
    logging.info("Wrote integrated h5ad: %s", out)

    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    logging.info("Wrote integration report: %s", report_path)


if __name__ == "__main__":
    main()
