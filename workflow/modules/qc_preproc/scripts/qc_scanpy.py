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
    )


def read_qc_config(path: str) -> dict:
    with open(path, "r") as f:
        cfg = yaml.safe_load(f)
    return cfg.get("qc", {})


def read_h5ad(path: str):
    logging.info("Reading h5ad: %s", path)
    return sc.read_h5ad(path)


def add_qc_metrics(adata, mito_prefix: str, mito_prefixes_extra: list[str]):
    prefixes = [mito_prefix] + list(mito_prefixes_extra or [])
    mito_mask = np.zeros(adata.n_vars, dtype=bool)
    for prefix in prefixes:
        mito_mask |= adata.var_names.str.upper().str.startswith(str(prefix).upper())
    adata.var["mt"] = mito_mask
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


def apply_filters(
    adata,
    min_genes: int,
    min_counts: int,
    max_mito_frac: float,
    min_cells_per_gene: int,
):
    before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    if "pct_counts_mt" in adata.obs.columns:
        adata = adata[adata.obs["pct_counts_mt"] <= (max_mito_frac * 100)].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    logging.info("Filtered cells: %s -> %s", before, adata.n_obs)
    return adata


def write_h5ad(adata, output_path: str):
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Writing output: %s", output_path)
    adata.write_h5ad(str(out), compression="gzip")

def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser()
    p.add_argument("--study", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--qc-config", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--log", default=None)
    args = p.parse_args(argv)

    setup_logging(args.log)
    qc = read_qc_config(args.qc_config)

    min_genes = int(qc.get("min_genes", 200))
    min_counts = int(qc.get("min_counts", 500))
    max_mito_frac = float(qc.get("max_mito_frac", 0.2))
    min_cells_per_gene = int(qc.get("min_cells_per_gene", 3))
    mito_prefix = str(qc.get("mito_prefix", "MT-"))
    mito_prefixes_extra = qc.get("mito_prefixes_extra", [])

    input_path = args.input
    suffix = Path(input_path).suffix.lower()
    if suffix != ".h5ad":
        raise SystemExit(f"qc_scanpy expects .h5ad input, got: {input_path}")

    adata = read_h5ad(input_path)

    adata.obs["study"] = args.study
    adata.obs["sample_id"] = args.sample

    add_qc_metrics(adata, mito_prefix=mito_prefix, mito_prefixes_extra=mito_prefixes_extra)
    adata = apply_filters(
        adata,
        min_genes=min_genes,
        min_counts=min_counts,
        max_mito_frac=max_mito_frac,
        min_cells_per_gene=min_cells_per_gene,
    )

    adata.uns.setdefault("qc_preproc", {})
    adata.uns["qc_preproc"]["config"] = qc
    adata.uns["qc_preproc"]["inputs"] = {
        "study": args.study,
        "sample": args.sample,
        "input": args.input,
    }

    write_h5ad(adata, args.output)


if __name__ == "__main__":
    main()
