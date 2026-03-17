"""Fetch helper for Snakemake.

Reads the project metadata CSV and downloads assets referenced by columns:

- `fastq_r1`, `fastq_r2` -> data_nature `fastqs`
- `raw_counts` -> data_nature `raw`
- `processed_counts` -> data_nature `qconly`

Files are written under:

`out/{study_accession}/{data_nature}/{donor_id}/{sample_id}/{run_id}/<files>`

This module intentionally does not attempt to unpack/reshape study-level
archives (e.g. GEO `*_RAW.tar`). It stores what it downloads deterministically;
downstream workflow steps can unpack and standardize as needed.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set


def _normalize_header(name: str) -> str:
    return re.sub(r"\s+", " ", name).strip()


def _open_csv(path: Path):
    return path.open(newline="", encoding="utf-8", errors="replace")


DATA_COLUMNS = {
    "fastqs": ("fastq_r1", "fastq_r2"),
    "raw": ("raw_counts",),
    "qconly": ("processed_counts",),
}


def _is_missing(value: str) -> bool:
    return (value or "").strip() == "" or (value or "").strip().lower() == "na"


def read_rows(metadata_csv: Path) -> List[Dict[str, str]]:
    with _open_csv(metadata_csv) as handle:
        reader = csv.reader(handle)
        try:
            raw_headers = next(reader)
        except StopIteration:
            raise ValueError(f"Empty CSV: {metadata_csv}")

        headers = [_normalize_header(h) for h in raw_headers]
        index = {h: i for i, h in enumerate(headers)}

        required = {"study_accession", "donor_id", "sample_id", "run_id"}
        missing = required - set(index)
        if missing:
            raise ValueError(
                f"Missing columns {sorted(missing)} in {metadata_csv}. Found: {headers}"
            )

        rows: List[Dict[str, str]] = []
        for raw_row in reader:
            if not raw_row:
                continue
            row: Dict[str, str] = {}
            for h, i in index.items():
                row[h] = raw_row[i].strip() if i < len(raw_row) else ""
            rows.append(row)
    return rows


def list_studies(metadata_csv: Path) -> List[str]:
    studies: Set[str] = set()
    for row in read_rows(metadata_csv):
        study = (row.get("study_accession") or "").strip()
        if study and not _is_missing(study):
            studies.add(study)
    return sorted(studies)

    with _open_csv(metadata_csv) as handle:
        reader = csv.reader(handle)
        try:
            raw_headers = next(reader)
        except StopIteration:
            raise ValueError(f"Empty CSV: {metadata_csv}")

        headers = [_normalize_header(h) for h in raw_headers]
        index = {h: i for i, h in enumerate(headers)}
        required = {"study_accession", "raw_uri", "processed_uri"}
        missing = required - set(index)
        if missing:
            raise ValueError(
                f"Missing columns {sorted(missing)} in {metadata_csv}. Found: {headers}"
            )

        studies: Dict[str, Dict[str, Set[str]]] = {}
        for row in reader:
            if not row:
                continue
            study = (row[index["study_accession"]] or "").strip()
            if not study:
                continue
            raw_uri = (row[index["raw_uri"]] or "").strip()
            processed_uri = (row[index["processed_uri"]] or "").strip()

            entry = studies.setdefault(study, {"raw": set(), "processed": set()})
            if raw_uri and raw_uri.lower() != "na":
                entry["raw"].add(raw_uri)
            if processed_uri and processed_uri.lower() != "na":
                entry["processed"].add(processed_uri)

    return studies
def _dest_dir(out_root: Path, row: Dict[str, str], data_nature: str) -> Path:
    return (
        out_root
        / row["study_accession"]
        / data_nature
        / row["donor_id"]
        / row["sample_id"]
        / row["run_id"]
    )


def _dest_from_url(url: str, outdir: Path) -> Path:
    name = url.rstrip("/").split("/")[-1]
    if not name:
        raise ValueError(f"Could not infer filename from url: {url}")
    return outdir / name


def _best_level_for_url(url: str, row: Dict[str, str]) -> str:
    """Return the most specific hierarchy level implied by `url`.

    Levels from most->least specific:
    run_id < sample_id < donor_id < study_accession

    If none of the tokens appear in the URL, returns "unknown".
    """

    url_lc = url.lower()
    for level in ("run_id", "sample_id", "donor_id", "study_accession"):
        token = (row.get(level) or "").strip()
        if token and not _is_missing(token) and token.lower() in url_lc:
            return level
    return "unknown"


def _dest_dir_with_level(out_root: Path, row: Dict[str, str], data_nature: str, level: str) -> Path:
    base = out_root / row["study_accession"] / data_nature
    if level == "study_accession":
        return base
    if level == "donor_id":
        return base / row["donor_id"]
    if level == "sample_id":
        return base / row["donor_id"] / row["sample_id"]
    if level == "run_id":
        return base / row["donor_id"] / row["sample_id"] / row["run_id"]
    raise ValueError(f"Unknown level: {level}")


def _curl_download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)

    # If a previous attempt created an empty file, remove it so curl can resume
    # cleanly.
    try:
        if dest.exists() and dest.stat().st_size == 0:
            dest.unlink()
    except OSError:
        pass

    cmd = [
        "curl",
        "-L",
        "--fail",
        "--retry",
        # Large GEO downloads occasionally fail mid-transfer; be more patient.
        "10",
        "--retry-delay",
        "2",
        "--retry-all-errors",
        "--connect-timeout",
        "30",
        "--max-time",
        "0",
        "-C",
        "-",
        "-o",
        str(dest),
        url,
    ]
    subprocess.run(cmd, check=True)


def cmd_download_study(metadata_csv: Path, out_root: Path, study: str) -> None:
    rows = [r for r in read_rows(metadata_csv) if r.get("study_accession") == study]
    if not rows:
        raise SystemExit(
            f"Study '{study}' not found in {metadata_csv}. Available: {list_studies(metadata_csv)}"
        )

    # Precompute a fallback "level" per URL based on how unique it is across
    # different hierarchy keys. This handles study-level links that do not
    # include tokens like sample_id/run_id in the URL.
    url_to_rows: Dict[str, List[Dict[str, str]]] = {}
    for row in rows:
        for _, cols in DATA_COLUMNS.items():
            for col in cols:
                url = (row.get(col) or "").strip()
                if _is_missing(url):
                    continue
                url_to_rows.setdefault(url, []).append(row)

    url_fallback_level: Dict[str, str] = {}
    for url, urows in url_to_rows.items():
        for level in ("run_id", "sample_id", "donor_id", "study_accession"):
            values = {(r.get(level) or "").strip() for r in urows}
            values = {v for v in values if v and not _is_missing(v)}
            if len(values) == 1:
                url_fallback_level[url] = level
                break
        else:
            url_fallback_level[url] = "study_accession"

    seen: Set[tuple[str, Path]] = set()
    for row in rows:
        for data_nature, cols in DATA_COLUMNS.items():
            for col in cols:
                url = (row.get(col) or "").strip()
                if _is_missing(url):
                    continue

                level = _best_level_for_url(url, row)
                if level == "unknown":
                    level = url_fallback_level.get(url, "study_accession")

                dest_dir = _dest_dir_with_level(out_root, row, data_nature, level)
                dest = _dest_from_url(url, dest_dir)
                key = (url, dest)
                if key in seen:
                    continue
                seen.add(key)
                _curl_download(url, dest)


def cmd_list_studies(metadata_csv: Path) -> None:
    for s in list_studies(metadata_csv):
        print(s)


def cmd_list_studies_json(metadata_csv: Path) -> None:
    import json

    print(json.dumps(list_studies(metadata_csv)))


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--metadata",
        default="~/phd_coding/a1_p1_qc/config/liver_scRNAseq_metadata.merged.csv",
        help="Metadata CSV with study_accession/raw_uri/processed_uri",
    )
    parser.add_argument(
        "--out-root",
        default="~/phd_coding/a1_p1_qc/data",
        help="Output root directory",
    )

    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list-studies")
    sub.add_parser("list-studies-json")
    dl = sub.add_parser("download-study")
    dl.add_argument("study")

    args = parser.parse_args(argv)
    metadata_csv = Path(os.path.expanduser(args.metadata))
    out_root = Path(os.path.expanduser(args.out_root))

    if args.cmd == "list-studies":
        cmd_list_studies(metadata_csv)
        return 0
    if args.cmd == "list-studies-json":
        cmd_list_studies_json(metadata_csv)
        return 0
    if args.cmd == "download-study":
        cmd_download_study(metadata_csv, out_root, args.study)
        return 0

    raise AssertionError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())
