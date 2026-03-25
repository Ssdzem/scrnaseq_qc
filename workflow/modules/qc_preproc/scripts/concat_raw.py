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
    )


def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser()
    p.add_argument("--inputs", nargs="+", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--log", default=None)
    args = p.parse_args(argv)

    setup_logging(args.log)

    adatas = []
    for h5 in args.inputs:
        logging.info("Reading: %s", h5)
        adata = sc.read_h5ad(h5)
        adatas.append(adata)

    logging.info("Concatenating %d studies", len(adatas))
    merged = ad.concat(adatas, join="outer", merge="same", label="study_source", index_unique="-")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Writing: %s", args.output)
    merged.write_h5ad(str(out), compression="gzip")


if __name__ == "__main__":
    main()
