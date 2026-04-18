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


def read_metadata_rows(path: str) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        return [{k: v for k, v in row.items()} for row in reader]


def build_rows_by_key(metadata_rows: list[dict[str, str]], key_column: str) -> dict[str, dict[str, str]]:
    rows_by_key: dict[str, dict[str, str]] = {}
    duplicate_keys: set[str] = set()
    for row in metadata_rows:
        key_value = (row.get(key_column) or "").strip()
        if not key_value:
            continue
        if key_value in rows_by_key:
            duplicate_keys.add(key_value)
            continue
        rows_by_key[key_value] = {k: v for k, v in row.items() if k != key_column}

    if duplicate_keys:
        logging.warning(
            "Ignoring duplicate rows for key column '%s' for %d IDs (first row kept): %s",
            key_column,
            len(duplicate_keys),
            ", ".join(sorted(duplicate_keys)[:10]),
        )
    return rows_by_key


def detect_fallback_key_column(adata_sample_ids: list[str], metadata_rows: list[dict[str, str]]) -> str | None:
    probe_ids = [v for v in adata_sample_ids[:5] if v]
    if not probe_ids or not metadata_rows:
        return None

    column_scores: dict[str, int] = {}
    for row in metadata_rows:
        for col, val in row.items():
            if val is None:
                continue
            if str(val).strip() in probe_ids:
                column_scores[col] = column_scores.get(col, 0) + 1

    if not column_scores:
        return None

    # Deterministic best key: most probe hits, then alphabetical column name.
    return sorted(column_scores.items(), key=lambda item: (-item[1], item[0]))[0][0]


def attach_metadata(
    adata,
    metadata_rows: list[dict[str, str]],
    metadata_path: str,
    primary_key_column: str = "sample_id",
):
    if "sample_id" not in adata.obs.columns:
        raise ValueError("AnnData obs is missing required 'sample_id' column for metadata merge.")

    sample_ids = adata.obs["sample_id"].astype(str).tolist()
    join_column = primary_key_column
    metadata_by_key = build_rows_by_key(metadata_rows, key_column=join_column)
    matched_rows = [metadata_by_key.get(sample_id, {}) for sample_id in sample_ids]

    if sum(1 for row in matched_rows if row) == 0:
        fallback_column = detect_fallback_key_column(sample_ids, metadata_rows)
        if fallback_column and fallback_column != join_column:
            logging.info(
                "No metadata matched by '%s'. Found fallback key column '%s' using first 5 adata sample IDs.",
                join_column,
                fallback_column,
            )
            join_column = fallback_column
            metadata_by_key = build_rows_by_key(metadata_rows, key_column=join_column)
            matched_rows = [metadata_by_key.get(sample_id, {}) for sample_id in sample_ids]
        elif fallback_column == join_column:
            logging.info("No metadata matched and fallback search points to the same key column '%s'.", join_column)
        else:
            logging.info(
                "No metadata matched by '%s' and no fallback key column matched first 5 adata sample IDs.",
                join_column,
            )

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
        "Attached metadata from %s by key='%s': matched_cells=%d/%d, new_columns=%d",
        metadata_path,
        join_column,
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
    metadata_rows = read_metadata_rows(metadata_csv)

    logging.info("Reading h5ad: %s", args.input)
    adata = sc.read_h5ad(args.input)
    adata = attach_metadata(adata, metadata_rows, metadata_path=metadata_csv, primary_key_column="sample_id")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Writing output: %s", args.output)
    adata.write_h5ad(str(out), compression="gzip")


if __name__ == "__main__":
    main()
