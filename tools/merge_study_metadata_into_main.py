#!/usr/bin/env python3
"""Merge per-study metadata into the main liver_scRNAseq_metadata.csv.

This script reads `config/liver_scRNAseq_metadata.csv` and enriches it with
selected variables (age, sex/gender, BMI) coming from study-specific metadata
tables located in the same `config/` directory.

Assumptions (based on current directory structure):
- Each study has a `{study_accession}_metadata.csv` file.
- Some studies also have a `{study_accession}_sup_metadata.csv` file containing
  extra rows/columns.

The merge key differs per study and is derived from the `donor_id` column of the
main CSV.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


CONFIG_DIR = Path(__file__).resolve().parents[1] / "config"


def _read_csv_rows(path: Path) -> List[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def _normalize_sex(value: str) -> str:
    v = (value or "").strip()
    if not v:
        return ""
    low = v.lower()
    if low in {"m", "male"}:
        return "male"
    if low in {"f", "female"}:
        return "female"
    if low in {"unknown", "na", "n/a"}:
        return ""
    return v


def _parse_float(value: str) -> str:
    v = (value or "").strip()
    if not v:
        return ""
    try:
        return str(float(v))
    except ValueError:
        return ""


def _parse_int(value: str) -> str:
    v = (value or "").strip()
    if not v:
        return ""
    try:
        return str(int(float(v)))
    except ValueError:
        return ""


@dataclass(frozen=True)
class StudyMapping:
    key_from_main_donor_id: callable
    extract_vars: callable


def _gse202379_key(donor_id: str) -> str:
    d = (donor_id or "").strip()
    if d.startswith("P") and len(d) > 1:
        return d[1:]
    if d.startswith("PHL"):
        return d.removeprefix("P")
    return d


def _gse202379_extract(row: dict) -> dict:
    # Works for both *_metadata and *_sup_metadata (different column names).
    age = row.get("Age") or row.get("age")
    sex = row.get("Sex") or row.get("Gender") or row.get("sex") or row.get("gender")
    bmi = row.get("BMI") or row.get("BMI (kg/m2)") or row.get("bmi")
    return {
        "age_years": _parse_int(age or ""),
        "sex": _normalize_sex(sex or ""),
        "bmi": _parse_float(bmi or ""),
    }


def _gse244832_key(donor_id: str) -> str:
    # Main donor_id examples: JB_288 etc. Study file has 'Donor ID snRNA/ATAC' with same.
    return (donor_id or "").strip()


_AGE_GENDER_RE = re.compile(r"^\s*(?P<age>\d+)\s*,\s*(?P<sex>[MF])\s*$", re.I)


def _gse244832_extract(row: dict) -> dict:
    age = ""
    sex = ""
    age_gender = (row.get("Age, Gender") or "").strip()
    match = _AGE_GENDER_RE.match(age_gender)
    if match:
        age = match.group("age")
        sex = match.group("sex")
    bmi = row.get("BMI")
    return {
        "age_years": _parse_int(age),
        "sex": _normalize_sex(sex),
        "bmi": _parse_float(bmi or ""),
    }


def _gse256398_key(donor_id: str) -> str:
    # Main donor_id examples: S19. Study table has 'Sample ID' that looks like S19.
    return (donor_id or "").strip()


def _gse256398_extract(row: dict) -> dict:
    return {
        "age_years": _parse_int(row.get("Age") or ""),
        "sex": _normalize_sex(row.get("Gender") or ""),
        "bmi": "",  # Not present in this supplementary table
    }


def _hra007511_key(donor_id: str) -> str:
    # Main donor_id examples: CTRL-173, MASH-0114. Study uses 'Sample name'.
    return (donor_id or "").strip()


def _hra007511_extract(row: dict) -> dict:
    return {
        "age_years": _parse_int(row.get("Age") or ""),
        "sex": _normalize_sex(row.get("Sex") or ""),
        "bmi": _parse_float((row.get("BMI (kg/m2)") or "").replace(",", "")),
    }


STUDY_MAPPINGS: Dict[str, StudyMapping] = {
    "GSE202379": StudyMapping(_gse202379_key, _gse202379_extract),
    "GSE244832": StudyMapping(_gse244832_key, _gse244832_extract),
    "GSE256398": StudyMapping(_gse256398_key, _gse256398_extract),
    "HRA007511": StudyMapping(_hra007511_key, _hra007511_extract),
}


def _build_lookup(study: str, config_dir: Path) -> Dict[str, dict]:
    mapping = STUDY_MAPPINGS[study]

    meta_path = config_dir / f"{study}_metadata.csv"
    if not meta_path.exists():
        raise FileNotFoundError(f"Missing per-study metadata file: {meta_path}")

    # Some metadata files include descriptive rows before the header.
    rows = _read_csv_rows(meta_path)
    if study in {"GSE256398", "HRA007511"}:
        with meta_path.open("r", encoding="utf-8-sig", newline="") as handle:
            lines = list(handle)

        if study == "GSE256398":
            header_idx = next(
                (
                    i
                    for i, line in enumerate(lines)
                    if "Sample ID" in line and "Gender" in line and "Age" in line
                ),
                None,
            )
        else:  # HRA007511
            header_idx = next(
                (i for i, line in enumerate(lines) if line.startswith("Sample name,Age,Sex")),
                None,
            )

        if header_idx is not None:
            reader = csv.DictReader(lines[header_idx:])
            rows = list(reader)

    sup_path = config_dir / f"{study}_sup_metadata.csv"
    sup_rows: List[dict] = []
    if sup_path.exists():
        sup_rows = _read_csv_rows(sup_path)

    lookup: Dict[str, dict] = {}

    def add_rows(rows_iter: Iterable[dict], key_fields: Sequence[str]) -> None:
        for row in rows_iter:
            key = ""
            for field in key_fields:
                if row.get(field):
                    key = str(row[field]).strip()
                    break
            if not key:
                continue
            lookup[key] = mapping.extract_vars(row)

    if study == "GSE202379":
        add_rows(rows, ["Patient ID"])
        add_rows(sup_rows, ["Patient_ID"])
    elif study == "GSE244832":
        add_rows(rows, ["Donor ID snRNA/ATAC", "Donor ID"])
    elif study == "GSE256398":
        add_rows(rows, ["   Sample ID", "Sample ID"])
    elif study == "HRA007511":
        add_rows(rows, ["Sample name"])
    else:
        raise KeyError(f"No lookup builder for study {study}")

    return lookup


def merge_main(
    main_path: Path,
    config_dir: Path,
    studies: Optional[Sequence[str]],
    overwrite_existing: bool,
) -> Tuple[List[dict], Dict[str, int]]:
    main_rows = _read_csv_rows(main_path)
    if not main_rows:
        return [], {}

    chosen_studies = sorted(set(studies) if studies else {r["study_accession"] for r in main_rows})
    for s in chosen_studies:
        if s not in STUDY_MAPPINGS:
            raise ValueError(
                f"Study {s} not supported yet. Supported: {sorted(STUDY_MAPPINGS)}"
            )

    lookups = {s: _build_lookup(s, config_dir) for s in chosen_studies}

    updated = 0
    unmatched = 0
    for row in main_rows:
        study = row.get("study_accession", "").strip()
        if study not in lookups:
            continue
        mapping = STUDY_MAPPINGS[study]
        key = mapping.key_from_main_donor_id(row.get("donor_id", ""))
        vars_row = lookups[study].get(key)
        if not vars_row:
            unmatched += 1
            continue

        changed_here = False
        for field, value in vars_row.items():
            if not value:
                continue
            if (row.get(field) or "").strip() and not overwrite_existing:
                continue
            if row.get(field) != value:
                row[field] = value
                changed_here = True
        if changed_here:
            updated += 1

    return main_rows, {"updated_rows": updated, "unmatched_rows": unmatched}


def _write_csv(path: Path, rows: List[dict]) -> None:
    if not rows:
        raise ValueError("No rows to write")
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Enrich liver_scRNAseq_metadata.csv with per-study metadata"
    )
    parser.add_argument(
        "--main",
        type=Path,
        default=CONFIG_DIR / "liver_scRNAseq_metadata.csv",
        help="Path to liver_scRNAseq_metadata.csv",
    )
    parser.add_argument(
        "--config-dir",
        type=Path,
        default=CONFIG_DIR,
        help="Directory containing per-study metadata CSVs",
    )
    parser.add_argument(
        "--studies",
        nargs="*",
        default=None,
        help="Subset of study_accession values to merge",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Write merged CSV to this path (default: print stats only)",
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Overwrite the main CSV in place (ignored if --out is provided)",
    )
    parser.add_argument(
        "--overwrite-existing",
        action="store_true",
        help="Overwrite non-empty values already present in main CSV",
    )
    args = parser.parse_args(argv)

    merged_rows, stats = merge_main(
        main_path=args.main,
        config_dir=args.config_dir,
        studies=args.studies,
        overwrite_existing=args.overwrite_existing,
    )
    print("Merge stats:")
    for k, v in stats.items():
        print(f"- {k}: {v}")

    if args.out is not None:
        _write_csv(args.out, merged_rows)
        print(f"Wrote: {args.out}")
    elif args.in_place:
        _write_csv(args.main, merged_rows)
        print(f"Updated in place: {args.main}")
    else:
        print("No output written (use --out or --in-place).")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
