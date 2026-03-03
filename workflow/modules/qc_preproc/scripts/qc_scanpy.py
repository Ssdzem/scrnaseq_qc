import argparse
import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
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


def read_10x_h5(path: str) -> ad.AnnData:
    logging.info("Reading 10x h5: %s", path)
    return sc.read_10x_h5(path)


def read_h5ad(path: str) -> ad.AnnData:
    logging.info("Reading h5ad: %s", path)
    return sc.read_h5ad(path)


def _read_table_1col(path: str) -> np.ndarray:
    p = Path(path)
    if p.suffix.lower() in {".tsv", ".txt"}:
        df = pd.read_csv(path, sep="\t", header=None)
    else:
        df = pd.read_csv(path, header=None)
    return df.iloc[:, 0].astype(str).values


def read_mtx(matrix_path: str, genes_path: str, barcodes_path: str) -> ad.AnnData:
    logging.info("Reading MTX: %s", matrix_path)
    if str(matrix_path).endswith(".gz"):
        import gzip

        with gzip.open(matrix_path, "rt") as f:
            mat = scipy.io.mmread(f)
    else:
        mat = scipy.io.mmread(matrix_path)

    if not sp.issparse(mat):
        mat = sp.csr_matrix(mat)
    mat = mat.tocsr()

    var_names = _read_table_1col(genes_path)
    obs_names = _read_table_1col(barcodes_path)

    if mat.shape[0] == len(var_names) and mat.shape[1] == len(obs_names):
        mat = mat.T

    adata = ad.AnnData(X=mat)
    adata.var_names = var_names
    adata.obs_names = obs_names
    return adata


def add_qc_metrics(adata: ad.AnnData, mito_prefix: str, mito_prefixes_extra: list[str]):
    prefixes = [mito_prefix] + list(mito_prefixes_extra or [])
    mito_mask = np.zeros(adata.n_vars, dtype=bool)
    for prefix in prefixes:
        mito_mask |= adata.var_names.str.upper().str.startswith(str(prefix).upper())
    adata.var["mt"] = mito_mask
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


def apply_filters(
    adata: ad.AnnData,
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


def write_as_h5(adata: ad.AnnData, output_path: str):
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Writing output: %s", output_path)
    adata.write_h5ad(str(out), compression="gzip")


def _json_list_arg(value: str) -> list[str]:
    if not value:
        return []
    try:
        parsed = json.loads(value)
        if parsed is None:
            return []
        if isinstance(parsed, list):
            return [str(x) for x in parsed if str(x).strip()]
        return [str(parsed)]
    except json.JSONDecodeError:
        return [value]


def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser()
    p.add_argument("--study", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--genes", default="[]", help="JSON list with 0/1 paths")
    p.add_argument("--barcodes", default="[]", help="JSON list with 0/1 paths")
    p.add_argument("--h5ad", default="[]", help="JSON list with 0/1 paths")
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

    h5ad_list = _json_list_arg(args.h5ad)
    genes_list = _json_list_arg(args.genes)
    barcodes_list = _json_list_arg(args.barcodes)

    h5ad_path = h5ad_list[0] if h5ad_list else ""
    genes_path = genes_list[0] if genes_list else ""
    barcodes_path = barcodes_list[0] if barcodes_list else ""

    if h5ad_path:
        adata = read_h5ad(h5ad_path)
    elif suffix == ".h5":
        adata = read_10x_h5(input_path)
    elif suffix == ".mtx" or input_path.lower().endswith(".mtx.gz"):
        if not genes_path or not barcodes_path:
            raise SystemExit("MTX input requires genes+barcodes paths")
        adata = read_mtx(input_path, genes_path, barcodes_path)
    else:
        raise SystemExit(f"Unsupported input format for qc_scanpy: {input_path}")

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
        "genes": genes_path,
        "barcodes": barcodes_path,
        "h5ad": h5ad_path,
    }

    write_as_h5(adata, args.output)


if __name__ == "__main__":
    main()
