import argparse
import json
import logging
from pathlib import Path

import anndata as ad
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


def read_full_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def resolve_ia_cfg(full_cfg: dict) -> dict:
    return full_cfg.get("integration", {})


def ensure_obs_key(adata, key: str, fallback_value: str):
    if key not in adata.obs.columns:
        adata.obs[key] = fallback_value


def resolve_accelerator(cfg: dict) -> tuple[str, bool]:
    accelerator_raw = str(cfg.get("accelerator", "")).strip().lower()
    use_gpu_raw = cfg.get("use_gpu", None)

    if accelerator_raw in {"cpu", "gpu"}:
        return accelerator_raw, accelerator_raw == "gpu"

    if use_gpu_raw is None:
        # Default to CPU for portability.
        return "cpu", False

    if isinstance(use_gpu_raw, bool):
        return ("gpu", True) if use_gpu_raw else ("cpu", False)

    use_gpu = str(use_gpu_raw).strip().lower() in {"1", "true", "yes", "y", "gpu"}
    return ("gpu", True) if use_gpu else ("cpu", False)


def run_integration(adata, cfg: dict, scope_label: str):
    integration_method = str(cfg.get("method", cfg.get("integration_method", "scvi"))).strip().lower()
    do_integration = bool(cfg.get("enabled", True))
    hvg_key = str(cfg.get("hvg_key", "qc_preproc.hvg_var_names"))
    n_pcs = int(cfg.get("n_pcs", 50))
    n_neighbors = int(cfg.get("n_neighbors", 15))
    leiden_resolution = float(cfg.get("leiden_resolution", 1.0))
    scvi_batch_key = str(cfg.get("batch_key", "sample_id"))
    scvi_max_epochs = int(cfg.get("scvi_max_epochs", 100))
    scvi_cpu_num_workers = int(cfg.get("scvi_cpu_num_workers", 8))
    scvi_cpu_num_threads = int(cfg.get("scvi_cpu_num_threads", 16))
    scvi_cpu_persistent_workers = bool(cfg.get("scvi_cpu_persistent_workers", True))
    seed = int(cfg.get("seed", 0))
    np.random.seed(seed)

    ensure_obs_key(adata, "sample_id", scope_label)
    ensure_obs_key(adata, scvi_batch_key, scope_label)

    # Prefer the precomputed HVGs from qc_preproc; fallback to all vars.
    hvg_var_names = []
    if hvg_key == "qc_preproc.hvg_var_names":
        hvg_var_names = list(adata.uns.get("qc_preproc", {}).get("hvg_var_names", []))
    if not hvg_var_names:
        hvg_var_names = adata.var_names.tolist()

    adata_hvg = adata[:, hvg_var_names].copy()
    adata_hvg.layers["log1p"] = adata_hvg.X.copy()

    sc.pp.scale(adata_hvg, max_value=10)
    adata_hvg.layers["scaled"] = adata_hvg.X.copy()

    sc.tl.pca(adata_hvg, n_comps=n_pcs, svd_solver="arpack")

    if do_integration and integration_method == "scvi":
        import scvi
        import torch

        accelerator, use_gpu_requested = resolve_accelerator(cfg)
        if accelerator == "gpu" and not torch.cuda.is_available():
            logging.warning("Integration configured with GPU but CUDA is unavailable; falling back to CPU")
            accelerator = "cpu"

        scvi.settings.seed = seed
        if accelerator == "cpu":
            scvi.settings.dl_num_workers = scvi_cpu_num_workers
            scvi.settings.num_threads = scvi_cpu_num_threads
            scvi.settings.dl_persistent_workers = scvi_cpu_persistent_workers

        layer_for_scvi = "counts" if "counts" in adata_hvg.layers else None
        scvi.model.SCVI.setup_anndata(adata_hvg, layer=layer_for_scvi, batch_key=scvi_batch_key)
        model = scvi.model.SCVI(adata_hvg)
        model.train(max_epochs=scvi_max_epochs, accelerator=accelerator)
        adata_hvg.obsm["X_integrated"] = model.get_latent_representation()
        adata_hvg.obsm["X_scVI"] = adata_hvg.obsm["X_integrated"]
        cluster_key = "leiden_scvi"
    elif do_integration and integration_method == "harmony":
        import harmonypy as hm

        ho = hm.run_harmony(adata_hvg.obsm["X_pca"], adata_hvg.obs, [scvi_batch_key])
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

    sc.pp.neighbors(adata_hvg, use_rep="X_integrated", n_neighbors=n_neighbors)
    sc.tl.umap(adata_hvg)
    sc.tl.leiden(adata_hvg, key_added=cluster_key, resolution=leiden_resolution)

    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.obsm["X_integrated"] = adata_hvg.obsm["X_integrated"]
    if "X_scVI" in adata_hvg.obsm:
        adata.obsm["X_scVI"] = adata_hvg.obsm["X_scVI"]
    if "X_harmony" in adata_hvg.obsm:
        adata.obsm["X_harmony"] = adata_hvg.obsm["X_harmony"]
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.obs["leiden"] = adata_hvg.obs[cluster_key].astype(str).reindex(adata.obs_names)
    adata.obs[cluster_key] = adata.obs["leiden"]

    adata.uns.setdefault("integration_annotation", {})
    adata.uns["integration_annotation"]["integration"] = {
        "method": integration_method,
        "enabled": do_integration,
        "n_pcs": n_pcs,
        "n_neighbors": n_neighbors,
        "leiden_resolution": leiden_resolution,
        "batch_key": scvi_batch_key,
        "scvi_max_epochs": scvi_max_epochs,
        "accelerator": accelerator if do_integration and integration_method == "scvi" else "",
        "use_gpu_requested": bool(use_gpu_requested) if do_integration and integration_method == "scvi" else False,
        "scvi_cpu_num_workers": scvi_cpu_num_workers,
        "scvi_cpu_num_threads": scvi_cpu_num_threads,
        "scvi_cpu_persistent_workers": scvi_cpu_persistent_workers,
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
