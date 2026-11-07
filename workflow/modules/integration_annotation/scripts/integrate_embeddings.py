import argparse
import json
import logging
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml


VALID_INTEGRATION_METHODS = {"scvi", "sysvi", "harmony"}


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


def _as_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y", "gpu"}
    return bool(value)


def _resolve_backend(cfg: dict) -> str:
    backend = str(cfg.get("backend", cfg.get("preprocessing_backend", "cpu"))).strip().lower()
    if backend not in {"cpu", "gpu", "auto"}:
        raise ValueError(f"Unsupported integration.backend={backend!r}. Choose one of: cpu, gpu, auto")
    return backend


def _resolve_rapids_runtime(*, backend: str, require_gpu: bool):
    if backend == "cpu":
        return None, False

    try:
        import cupy as cp
        import rapids_singlecell as rsc
    except ImportError as exc:
        if backend == "gpu" or require_gpu:
            raise RuntimeError("integration GPU backend requested but RAPIDS SingleCell is not installed") from exc
        logging.warning("RAPIDS SingleCell unavailable; falling back to CPU graph/integration steps: %s", exc)
        return None, False

    try:
        gpu_count = int(cp.cuda.runtime.getDeviceCount())
    except Exception as exc:  # pragma: no cover - depends on host CUDA runtime
        if backend == "gpu" or require_gpu:
            raise RuntimeError("integration GPU backend requested but CUDA runtime is unavailable") from exc
        logging.warning("CUDA runtime unavailable; falling back to CPU graph/integration steps: %s", exc)
        return None, False

    if gpu_count < 1:
        if backend == "gpu" or require_gpu:
            raise RuntimeError("integration GPU backend requested but no CUDA devices were detected")
        logging.warning("No CUDA devices detected; falling back to CPU graph/integration steps")
        return None, False

    return rsc, True


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


def _ensure_composite_obs_key(adata, key: str, components: list[str], sep: str = "__", unknown_label: str = "Unknown"):
    if not components:
        return
    _require_obs_keys(adata, components)
    parts = []
    for component in components:
        series = adata.obs[component]
        series = series.where(~series.isna(), other=unknown_label)
        series = series.astype(str).replace({"": unknown_label, "nan": unknown_label, "None": unknown_label})
        parts.append(series)
    composite = parts[0]
    for series in parts[1:]:
        composite = composite + sep + series
    if key in adata.obs.columns and not adata.obs[key].astype(str).equals(composite.astype(str)):
        raise ValueError(
            f"Refusing to overwrite obs['{key}'] with a composite of {components}. "
            "Use a new column name for the key, or leave the components empty to use "
            f"obs['{key}'] as it is."
        )
    adata.obs[key] = pd.Categorical(composite)


def _materialize_obs_aliases(adata, alias_cfg, unknown_label: str = "Unknown"):
    if not isinstance(alias_cfg, dict):
        return

    for target_key, spec in alias_cfg.items():
        target_key = str(target_key).strip()
        if not target_key or target_key in adata.obs.columns:
            continue

        if isinstance(spec, dict):
            components = _as_list(spec.get("components", []))
            if components and all(component in adata.obs.columns for component in components):
                _ensure_composite_obs_key(
                    adata,
                    target_key,
                    components,
                    sep=str(spec.get("separator", "__")),
                    unknown_label=unknown_label,
                )
                continue
            candidates = _as_list(spec.get("candidates", spec.get("source", [])))
        else:
            candidates = _as_list(spec)

        for candidate in candidates:
            if candidate not in adata.obs.columns:
                continue
            series = adata.obs[candidate]
            series = series.where(~series.isna(), other=unknown_label)
            series = series.astype(str).replace({"": unknown_label, "nan": unknown_label, "None": unknown_label})
            adata.obs[target_key] = pd.Categorical(series)
            break


def _resolve_integration_method(cfg: dict, enabled: bool) -> str:
    method = str(cfg.get("method", cfg.get("integration_method", ""))).strip().lower()
    if not enabled:
        return "none"
    if not method:
        raise ValueError(
            "integration.method must be explicitly set before running integration_annotation. "
            "Choose one of: sysvi, scvi, harmony."
        )
    if method not in VALID_INTEGRATION_METHODS:
        raise ValueError(
            f"Unsupported integration.method={method!r}. Choose one of: {sorted(VALID_INTEGRATION_METHODS)}"
        )
    return method


def _prepare_hvg_view(adata, hvg_key: str):
    hvg_var_names = []
    if hvg_key == "qc_preproc.hvg_var_names":
        hvg_var_names = list(adata.uns.get("qc_preproc", {}).get("hvg_var_names", []))
    if not hvg_var_names:
        hvg_var_names = adata.var_names.tolist()

    hvg_var_names = [gene for gene in hvg_var_names if gene in adata.var_names]
    if not hvg_var_names:
        raise ValueError("No valid HVG genes found for integration")
    return adata[:, hvg_var_names].copy(), hvg_var_names


def _compute_pca(adata_hvg, n_pcs: int):
    n_comps_eff = min(n_pcs, max(1, adata_hvg.n_obs - 1), max(1, adata_hvg.n_vars - 1))
    if n_comps_eff < 2:
        raise ValueError(
            f"Cannot compute PCA with n_obs={adata_hvg.n_obs}, n_vars={adata_hvg.n_vars}; need at least 2 components"
        )
    pca_source = adata_hvg.copy()
    sc.pp.scale(pca_source, max_value=10)
    sc.tl.pca(pca_source, n_comps=n_comps_eff, svd_solver="arpack")
    adata_hvg.obsm["X_pca"] = pca_source.obsm["X_pca"]
    adata_hvg.uns["pca"] = pca_source.uns["pca"]
    adata_hvg.varm["PCs"] = pca_source.varm["PCs"]
    return n_comps_eff


def run_integration(adata, cfg: dict, scope_label: str):
    do_integration = _as_bool(cfg.get("enabled", True))
    integration_method = _resolve_integration_method(cfg, do_integration)
    backend = _resolve_backend(cfg)
    hvg_key = str(cfg.get("hvg_key", "qc_preproc.hvg_var_names")).strip()
    n_pcs = int(cfg.get("n_pcs", 50))
    n_neighbors = int(cfg.get("n_neighbors", 15))
    leiden_resolution = float(cfg.get("leiden_resolution", 1.0))
    scvi_batch_key = str(cfg.get("batch_key", "sample_id")).strip()
    categorical_covariate_keys = _as_list(cfg.get("categorical_covariate_keys", []))
    continuous_covariate_keys = _as_list(cfg.get("continuous_covariate_keys", []))
    system_key_components = _as_list(cfg.get("system_key_components", []))
    system_key_separator = str(cfg.get("system_key_separator", "__"))
    obs_key_aliases = cfg.get("obs_key_aliases", {}) or {}
    batch_keys_for_reporting = [
        *system_key_components,
        scvi_batch_key,
        *categorical_covariate_keys,
        *continuous_covariate_keys,
    ]
    compute_pca = _as_bool(cfg.get("compute_pca", integration_method in {"harmony", "none"}))
    compute_umap = _as_bool(cfg.get("compute_umap", True))
    compute_leiden = _as_bool(cfg.get("compute_leiden", True))
    scvi_max_epochs = int(cfg.get("scvi_max_epochs", 100))
    sysvi_max_epochs = int(cfg.get("sysvi_max_epochs", cfg.get("scvi_max_epochs", 100)))
    n_latent = int(cfg.get("n_latent", 30))
    batch_size = int(cfg.get("batch_size", 4096))
    validation_size = float(cfg.get("validation_size", 0.1))
    check_val_every_n_epoch = int(cfg.get("check_val_every_n_epoch", 1))
    early_stopping = _as_bool(cfg.get("early_stopping", True))
    early_stopping_monitor = str(cfg.get("early_stopping_monitor", "elbo_validation"))
    early_stopping_patience = int(cfg.get("early_stopping_patience", 20))
    dl_num_workers = int(cfg.get("dl_num_workers", 4))
    num_threads = int(os.environ.get("SCRNASEQ_QC_IA_NUM_THREADS", cfg.get("num_threads", 16)))
    use_gpu = _as_bool(cfg.get("use_gpu", False))
    require_gpu = _as_bool(cfg.get("require_gpu", False))
    accelerator = str(cfg.get("accelerator", "gpu" if use_gpu else "cpu")).strip().lower()
    devices = _coerce_device_value(cfg.get("devices", 1 if accelerator != "auto" else "auto"))
    seed = int(cfg.get("seed", 0))
    rsc, use_rapids = _resolve_rapids_runtime(backend=backend, require_gpu=require_gpu)
    np.random.seed(seed)

    ensure_obs_key(adata, "sample_id", scope_label)
    _materialize_obs_aliases(adata, obs_key_aliases)
    if system_key_components:
        _ensure_composite_obs_key(adata, scvi_batch_key, system_key_components, sep=system_key_separator)
    elif scvi_batch_key:
        ensure_obs_key(adata, scvi_batch_key, scope_label)
    _require_obs_keys(adata, batch_keys_for_reporting)
    _sanitize_categorical_obs(adata, [scvi_batch_key, *categorical_covariate_keys])
    if adata.obs[scvi_batch_key].nunique() < 2:
        logging.warning("obs['%s'] has one value, so integration corrects nothing", scvi_batch_key)

    adata_hvg, hvg_var_names = _prepare_hvg_view(adata, hvg_key)
    if "log1p" not in adata_hvg.layers:
        adata_hvg.layers["log1p"] = adata_hvg.X.copy()

    configured_layer = cfg.get("layer")
    if configured_layer is not None:
        configured_layer = str(configured_layer).strip() or None
    if integration_method == "sysvi":
        layer_for_model = configured_layer if configured_layer is not None else "log1p"
        if layer_for_model is not None and layer_for_model not in adata_hvg.layers:
            raise ValueError(f"Configured sysVI layer not found: {layer_for_model}")
    elif integration_method == "scvi":
        if configured_layer is not None:
            layer_for_model = configured_layer
            if layer_for_model not in adata_hvg.layers:
                raise ValueError(f"Configured scVI layer not found: {layer_for_model}")
        elif "counts" in adata_hvg.layers:
            layer_for_model = "counts"
        else:
            layer_for_model = None
    else:
        layer_for_model = None

    pca_n_comps = None
    if compute_pca or integration_method in {"harmony", "none"}:
        if use_rapids:
            rsc.get.anndata_to_GPU(adata_hvg, convert_all=True)
            rsc.pp.pca(adata_hvg, n_comps=n_pcs, random_state=seed)
            pca_n_comps = int(adata_hvg.obsm["X_pca"].shape[1])
        else:
            pca_n_comps = _compute_pca(adata_hvg, n_pcs)

    if use_rapids and integration_method in {"scvi", "sysvi"}:
        # scvi-tools models operate on CPU-hosted AnnData / torch tensors; move any
        # RAPIDS-computed PCA results back before model setup/training.
        rsc.get.anndata_to_CPU(adata_hvg, convert_all=True)

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
            layer=layer_for_model,
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
            layer_for_model or "X",
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
        adata_hvg.obsm["X_integrated"] = model.get_latent_representation(batch_size=batch_size)
        adata_hvg.obsm["X_scVI"] = adata_hvg.obsm["X_integrated"]
        cluster_key = "leiden_scvi"
    elif do_integration and integration_method == "sysvi":
        import scvi
        import torch
        from scvi.external import SysVI

        scvi.settings.seed = seed
        if hasattr(torch, "set_float32_matmul_precision"):
            torch.set_float32_matmul_precision(str(cfg.get("float32_matmul_precision", "high")))
        cuda_available = bool(torch.cuda.is_available())
        if require_gpu and accelerator == "gpu" and not cuda_available:
            raise RuntimeError("integration.require_gpu=true but torch.cuda.is_available() is false")

        sysvi_cfg = cfg.get("sysvi", {}) or {}
        prior = str(sysvi_cfg.get("prior", "vamp"))
        n_prior_components = int(sysvi_cfg.get("n_prior_components", 20))
        trainable_priors = bool(sysvi_cfg.get("trainable_priors", True))
        n_layers = int(sysvi_cfg.get("n_layers", 2))
        n_hidden = int(sysvi_cfg.get("n_hidden", 256))
        dropout_rate = float(sysvi_cfg.get("dropout_rate", 0.1))
        out_var_mode = str(sysvi_cfg.get("out_var_mode", "feature"))
        embed_categorical_covariates = bool(sysvi_cfg.get("embed_categorical_covariates", True))
        weight_batches = bool(sysvi_cfg.get("weight_batches", True))
        get_latent_mc_samples = int(sysvi_cfg.get("get_latent_mc_samples", 5000))
        plan_kwargs = dict(sysvi_cfg.get("plan_kwargs", {}) or {})

        SysVI.setup_anndata(
            adata_hvg,
            layer=layer_for_model,
            batch_key=scvi_batch_key,
            categorical_covariate_keys=categorical_covariate_keys or None,
            continuous_covariate_keys=continuous_covariate_keys or None,
            weight_batches=weight_batches,
        )
        model = SysVI(
            adata_hvg,
            prior=prior,
            n_prior_components=n_prior_components,
            trainable_priors=trainable_priors,
            n_latent=n_latent,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            out_var_mode=out_var_mode,
            embed_categorical_covariates=embed_categorical_covariates,
        )
        scvi.settings.dl_num_workers = dl_num_workers
        scvi.settings.num_threads = num_threads
        scvi.settings.dl_persistent_workers = dl_num_workers > 0
        logging.info(
            "Training sysVI: system_key=%s system_components=%s categorical_covariates=%s layer=%s "
            "accelerator=%s devices=%s cuda_available=%s batch_size=%d max_epochs=%d prior=%s n_prior_components=%d "
            "n_latent=%d n_hidden=%d n_layers=%d weight_batches=%s",
            scvi_batch_key,
            system_key_components,
            categorical_covariate_keys,
            layer_for_model or "X",
            accelerator,
            devices,
            cuda_available,
            batch_size,
            sysvi_max_epochs,
            prior,
            n_prior_components,
            n_latent,
            n_hidden,
            n_layers,
            weight_batches,
        )
        model.train(
            max_epochs=sysvi_max_epochs,
            validation_size=validation_size,
            accelerator=accelerator,
            devices=devices,
            batch_size=batch_size,
            check_val_every_n_epoch=check_val_every_n_epoch,
            early_stopping=early_stopping,
            early_stopping_monitor=early_stopping_monitor,
            early_stopping_patience=early_stopping_patience,
            plan_kwargs=plan_kwargs or None,
        )
        adata_hvg.obsm["X_integrated"] = model.get_latent_representation(
            batch_size=batch_size,
            mc_samples=get_latent_mc_samples,
        )
        adata_hvg.obsm["X_sysVI"] = adata_hvg.obsm["X_integrated"]
        cluster_key = "leiden_sysvi"
    elif do_integration and integration_method == "harmony":
        if use_rapids:
            rsc.pp.harmony_integrate(adata_hvg, scvi_batch_key, basis="X_pca", adjusted_basis="X_harmony")
            adata_hvg.obsm["X_integrated"] = adata_hvg.obsm["X_harmony"]
        else:
            import harmonypy as hm

            ho = hm.run_harmony(adata_hvg.obsm["X_pca"], adata_hvg.obs, [scvi_batch_key])
            z_corr = np.asarray(ho.Z_corr)
            expected_pcs = adata_hvg.obsm["X_pca"].shape[1]
            if z_corr.shape == (expected_pcs, adata_hvg.n_obs):
                z_corr = z_corr.T
            if z_corr.shape != (adata_hvg.n_obs, expected_pcs):
                raise ValueError(f"Unexpected Harmony shape: {z_corr.shape}")
            adata_hvg.obsm["X_integrated"] = z_corr
            adata_hvg.obsm["X_harmony"] = z_corr
        cluster_key = "leiden_harmony"
    else:
        adata_hvg.obsm["X_integrated"] = adata_hvg.obsm["X_pca"]
        integration_method = "none"
        cluster_key = "leiden_pca"

    if use_rapids:
        rsc.get.anndata_to_GPU(adata_hvg, convert_all=True)
        if "X_integrated" not in adata_hvg.obsm:
            raise ValueError("X_integrated missing before RAPIDS neighbors computation")
        rsc.pp.neighbors(adata_hvg, use_rep="X_integrated", n_neighbors=n_neighbors, random_state=seed)
    else:
        sc.pp.neighbors(adata_hvg, use_rep="X_integrated", n_neighbors=n_neighbors)
    if compute_umap:
        if use_rapids:
            rsc.tl.umap(adata_hvg, random_state=seed)
        else:
            sc.tl.umap(adata_hvg)
    if compute_leiden:
        if use_rapids:
            rsc.tl.leiden(adata_hvg, key_added=cluster_key, resolution=leiden_resolution, random_state=seed)
        else:
            sc.tl.leiden(adata_hvg, key_added=cluster_key, resolution=leiden_resolution)

    if use_rapids:
        rsc.get.anndata_to_CPU(adata_hvg, convert_all=True)

    if "X_pca" in adata_hvg.obsm:
        adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.obsm["X_integrated"] = adata_hvg.obsm["X_integrated"]
    if "X_scVI" in adata_hvg.obsm:
        adata.obsm["X_scVI"] = adata_hvg.obsm["X_scVI"]
    if "X_sysVI" in adata_hvg.obsm:
        adata.obsm["X_sysVI"] = adata_hvg.obsm["X_sysVI"]
    if "X_harmony" in adata_hvg.obsm:
        adata.obsm["X_harmony"] = adata_hvg.obsm["X_harmony"]
    if "X_umap" in adata_hvg.obsm:
        adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    if cluster_key in adata_hvg.obs:
        adata.obs["leiden"] = adata_hvg.obs[cluster_key].astype(str).reindex(adata.obs_names)
        adata.obs[cluster_key] = adata.obs["leiden"]

    adata.uns.setdefault("integration_annotation", {})
    adata.uns["integration_annotation"]["integration"] = {
        "method": integration_method,
        "enabled": do_integration,
        "n_pcs": n_pcs,
        "pca_n_comps": pca_n_comps,
        "n_neighbors": n_neighbors,
        "leiden_resolution": leiden_resolution,
        "compute_pca": compute_pca,
        "compute_umap": compute_umap,
        "compute_leiden": compute_leiden,
        "batch_key": scvi_batch_key,
        "categorical_covariate_keys": categorical_covariate_keys,
        "continuous_covariate_keys": continuous_covariate_keys,
        "system_key_components": system_key_components,
        "system_key_separator": system_key_separator,
        "obs_key_aliases": obs_key_aliases,
        "scvi_max_epochs": scvi_max_epochs,
        "sysvi_max_epochs": sysvi_max_epochs,
        "n_latent": n_latent,
        "batch_size": batch_size,
        "use_gpu": use_gpu,
        "require_gpu": require_gpu,
        "accelerator": accelerator,
        "devices": devices,
        "layer": layer_for_model or "X",
        "backend": backend,
        "used_rapids_backend": use_rapids,
        "sysvi": cfg.get("sysvi", {}) or {},
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
