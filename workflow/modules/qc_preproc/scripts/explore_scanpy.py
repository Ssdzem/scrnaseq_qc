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


def _as_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y", "gpu"}
    return bool(value)


def resolve_backend(cfg: dict) -> tuple[str, bool]:
    backend = str(cfg.get("backend", cfg.get("preprocessing_backend", "cpu"))).strip().lower()
    if backend not in {"cpu", "gpu", "auto"}:
        raise ValueError(f"Unsupported qc.backend={backend!r}; choose one of: cpu, gpu, auto")
    require_gpu = _as_bool(cfg.get("require_gpu", False))
    return backend, require_gpu


def ensure_obs_key(adata, key: str, fallback_value: str):
    if key not in adata.obs.columns:
        adata.obs[key] = fallback_value


def _resolve_rapids_runtime(*, backend: str, require_gpu: bool):
    if backend == "cpu":
        return None, False

    try:
        import cupy as cp
        import rapids_singlecell as rsc
    except ImportError as exc:
        if backend == "gpu" or require_gpu:
            raise RuntimeError("qc GPU preprocessing requested but RAPIDS SingleCell is not installed") from exc
        logging.warning("RAPIDS SingleCell not available; falling back to CPU preprocessing: %s", exc)
        return None, False

    try:
        gpu_count = int(cp.cuda.runtime.getDeviceCount())
    except Exception as exc:  # pragma: no cover - depends on local CUDA runtime
        if backend == "gpu" or require_gpu:
            raise RuntimeError("qc GPU preprocessing requested but CUDA runtime is unavailable") from exc
        logging.warning("CUDA runtime unavailable; falling back to CPU preprocessing: %s", exc)
        return None, False

    if gpu_count < 1:
        if backend == "gpu" or require_gpu:
            raise RuntimeError("qc GPU preprocessing requested but no CUDA devices were detected")
        logging.warning("No CUDA devices detected; falling back to CPU preprocessing")
        return None, False

    return rsc, True


def prepare_before_integration(adata, study: str, cfg: dict):
    batch_key = str(cfg.get("hvg_batch_key", "sample_id"))
    scvi_batch_key = str(cfg.get("scvi_batch_key", "sample_id"))
    n_top_genes = int(cfg.get("hvg_n_top_genes", 3000))
    seed = int(cfg.get("seed", 0))
    backend, require_gpu = resolve_backend(cfg)
    rsc, using_gpu = _resolve_rapids_runtime(backend=backend, require_gpu=require_gpu)
    np.random.seed(seed)
    logging.info(
        "Pre-integration config: batch_key=%s scvi_batch_key=%s n_top_genes=%d seed=%d backend=%s using_gpu=%s",
        batch_key,
        scvi_batch_key,
        n_top_genes,
        seed,
        backend,
        using_gpu,
    )

    ensure_obs_key(adata, "study", study)
    ensure_obs_key(adata, "sample_id", study)
    ensure_obs_key(adata, batch_key, study)
    ensure_obs_key(adata, scvi_batch_key, study)

    if using_gpu:
        rsc.get.anndata_to_GPU(adata)
        logging.info("Moved AnnData .X to GPU memory")

    adata.layers["counts"] = adata.X.copy()
    logging.info("Created counts layer")

    pp = rsc.pp if using_gpu else sc.pp
    pp.normalize_total(adata, target_sum=1e4)
    adata.layers["normalized"] = adata.X.copy()
    logging.info("Completed normalize_total and stored normalized layer")

    pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()
    logging.info("Completed log1p and stored log1p layer")

    try:
        pp.highly_variable_genes(
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
        pp.highly_variable_genes(
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
        "backend": backend,
        "used_gpu": using_gpu,
    }

    adata.uns["qc_preproc"]["hvg_var_names"] = adata.var_names[adata.var["highly_variable"]].tolist()
    logging.info(
        "Prepared data before integration: n_cells=%d n_vars=%d n_hvg=%d",
        adata.n_obs,
        adata.n_vars,
        len(adata.uns["qc_preproc"]["hvg_var_names"]),
    )
    if using_gpu:
        rsc.get.anndata_to_CPU(adata, convert_all=True)
        logging.info("Moved AnnData matrices and arrays back to CPU memory")
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
