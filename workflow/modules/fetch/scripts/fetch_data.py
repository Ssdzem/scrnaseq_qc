"""Fetch helper for Snakemake.

Reads the project metadata CSV, groups entries by `study_accession`, and
downloads study-level `raw_uri` and `processed_uri` files into the expected
directory layout.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set


def _normalize_header(name: str) -> str:
    return re.sub(r"\s+", " ", name).strip()


def _open_csv(path: Path):
    return path.open(newline="", encoding="utf-8", errors="replace")


def read_study_uris(metadata_csv: Path) -> Dict[str, Dict[str, Set[str]]]:
    """Return mapping: study -> {raw, processed} -> set(urls)."""

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


def _dest_from_url(url: str, outdir: Path) -> Path:
    name = url.rstrip("/").split("/")[-1]
    if not name:
        raise ValueError(f"Could not infer filename from url: {url}")
    return outdir / name


def _curl_download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "curl",
        "-L",
        "--fail",
        "--retry",
        "3",
        "--retry-delay",
        "2",
        "-C",
        "-",
        "-o",
        str(dest),
        url,
    ]
    subprocess.run(cmd, check=True)


def cmd_download_study(metadata_csv: Path, out_root: Path, study: str) -> None:
    studies = read_study_uris(metadata_csv)
    if study not in studies:
        raise SystemExit(
            f"Study '{study}' not found in {metadata_csv}. Available: {sorted(studies)}"
        )

    outdir = out_root / study
    for url in sorted(studies[study]["raw"]):
        _curl_download(url, _dest_from_url(url, outdir / "raw"))
    for url in sorted(studies[study]["processed"]):
        _curl_download(url, _dest_from_url(url, outdir / "processed"))


def cmd_list_studies(metadata_csv: Path) -> None:
    studies = read_study_uris(metadata_csv)
    for s in sorted(studies):
        print(s)


def cmd_list_studies_json(metadata_csv: Path) -> None:
    import json

    studies = read_study_uris(metadata_csv)
    print(json.dumps(sorted(studies)))


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
