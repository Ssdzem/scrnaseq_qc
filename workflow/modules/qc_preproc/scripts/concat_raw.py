import argparse
import logging
from pathlib import Path

import anndata as ad
import scanpy as sc


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


def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser()
    p.add_argument("--inputs", nargs="+", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--log", default=None)
    args = p.parse_args(argv)

    setup_logging(args.log)
    logging.info("concat_raw args: %s", vars(args))

    adatas = []
    for h5 in args.inputs:
        logging.info("Reading: %s", h5)
        adata = sc.read_h5ad(h5)
        logging.info("Loaded %s with n_cells=%d n_genes=%d", h5, adata.n_obs, adata.n_vars)
        adatas.append(adata)

    logging.info("Concatenating %d studies", len(adatas))
    merged = ad.concat(adatas, join="outer", merge="same", label="study_source", index_unique="-")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Writing: %s", args.output)
    merged.write_h5ad(str(out), compression="gzip")
    logging.info("Wrote merged h5ad: %s (n_cells=%d n_genes=%d)", args.output, merged.n_obs, merged.n_vars)


if __name__ == "__main__":
    main()
