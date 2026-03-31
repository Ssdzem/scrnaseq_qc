import argparse
import csv
import logging
from pathlib import Path

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
    with open(path, "r") as f:
        return yaml.safe_load(f) or {}


def read_metadata_csv(path: str) -> dict[str, dict[str, str]]:
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        rows_by_sample: dict[str, dict[str, str]] = {}
        duplicate_sample_ids: set[str] = set()
        for row in reader:
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id:
                continue
            if sample_id in rows_by_sample:
                duplicate_sample_ids.add(sample_id)
                continue
            rows_by_sample[sample_id] = {k: v for k, v in row.items() if k != "sample_id"}

    if duplicate_sample_ids:
        logging.warning(
            "Ignoring duplicate sample_id rows in metadata CSV for %d IDs (first row kept): %s",
            len(duplicate_sample_ids),
            ", ".join(sorted(duplicate_sample_ids)[:10]),
        )
    return rows_by_sample


def attach_metadata(adata, metadata_by_sample: dict[str, dict[str, str]], metadata_path: str):
    if "sample_id" not in adata.obs.columns:
        raise ValueError("AnnData obs is missing required 'sample_id' column for metadata merge.")

    sample_ids = adata.obs["sample_id"].astype(str)
    matched_rows = [metadata_by_sample.get(sample_id, {}) for sample_id in sample_ids]
    metadata_columns = sorted({key for row in matched_rows for key in row.keys()})

    if not metadata_columns:
        logging.warning("No metadata columns found to attach from %s", metadata_path)
        return adata

    existing_cols = set(adata.obs.columns)
    for col in metadata_columns:
        if col in adata.obs.columns:
            continue
        adata.obs[col] = [row.get(col, "") for row in matched_rows]

    matched_count = sum(1 for row in matched_rows if row)
    logging.info(
        "Attached metadata from %s by sample_id: matched_cells=%d/%d, new_columns=%d",
        metadata_path,
        matched_count,
        adata.n_obs,
        len(set(adata.obs.columns) - existing_cols),
    )
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
    logging.info("attach_metadata args: %s", vars(args))

    cfg = read_full_config(args.qc_config)
    metadata_csv = cfg.get("inputs", {}).get("metadata_csv", "config/liver_scRNAseq_metadata_everything.csv")
    metadata_by_sample = read_metadata_csv(metadata_csv)

    logging.info("Reading h5ad: %s", args.input)
    adata = sc.read_h5ad(args.input)
    adata = attach_metadata(adata, metadata_by_sample, metadata_path=metadata_csv)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Writing output: %s", args.output)
    adata.write_h5ad(str(out), compression="gzip")


if __name__ == "__main__":
    main()
