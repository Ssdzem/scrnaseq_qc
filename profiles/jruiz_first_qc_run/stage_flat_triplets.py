#!/usr/bin/env python3
"""Stage flat 10x MTX deposits into the per-run triplet layout the standardize
module expects.

`count_matrix_10x_multi_runs` resolves each run to a directory and then looks for
fixed filenames inside it (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`).
GSE197461 and S-BSST1035 instead ship every sample side by side in one flat
`counts/` directory with a per-sample filename prefix, so this script builds

    <triplet-root>/<run_id>/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}

from `<source-dir>/<file_prefix>{matrix.mtx.gz,...}`, driven by the same runs CSV
the study spec points at (the extra `file_prefix` column is ignored by the
cleaner's own `_load_runs`).

The matrix is symlinked (they run to hundreds of MB); features and barcodes are
*copied*. That asymmetry is deliberate: `count_matrix_10x` rewrites
`features.tsv.gz` in place to guarantee a 3-column v3 file, and writing through a
symlink would corrupt the read-only deposit under raw_data/.

S-BSST1035 ships no barcode list at all; nothing is staged for it and the cleaner
fabricates positional barcodes into the staging dir when the study spec sets
`input.synthesize_missing_barcodes: true`.

Idempotent: re-running refreshes links/copies in place.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
from pathlib import Path

TRIPLET_FILES = {
    "matrix": ("matrix.mtx.gz", "symlink"),
    "features": ("features.tsv.gz", "copy"),
    "barcodes": ("barcodes.tsv.gz", "copy"),
}


def _stage_symlink(src: Path, dst: Path) -> None:
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    dst.symlink_to(os.path.relpath(src, dst.parent))


def _stage_copy(src: Path, dst: Path) -> None:
    if dst.is_symlink():
        dst.unlink()
    shutil.copyfile(src, dst)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--runs-csv", required=True, help="runs CSV with run_id + file_prefix columns")
    parser.add_argument("--study", required=True, help="value to match in study_accession")
    parser.add_argument("--source-dir", required=True, help="flat deposit directory holding every sample")
    parser.add_argument("--triplet-root", required=True, help="destination root; one subdir per run_id")
    args = parser.parse_args()

    source_dir = Path(args.source_dir)
    if not source_dir.is_dir():
        print(f"ERROR: --source-dir not found: {source_dir}", file=sys.stderr)
        return 1

    triplet_root = Path(args.triplet_root)
    triplet_root.mkdir(parents=True, exist_ok=True)

    with open(args.runs_csv, newline="", encoding="utf-8") as handle:
        rows = [r for r in csv.DictReader(handle) if r.get("study_accession") == args.study]

    if not rows:
        print(f"ERROR: no rows for study_accession={args.study} in {args.runs_csv}", file=sys.stderr)
        return 1

    staged = 0
    for row in rows:
        run_id = row["run_id"]
        prefix = row.get("file_prefix", "")
        run_dir = triplet_root / run_id
        run_dir.mkdir(parents=True, exist_ok=True)

        for kind, (canonical, mode) in TRIPLET_FILES.items():
            src = source_dir / f"{prefix}{canonical}"
            dst = run_dir / canonical
            if not src.exists():
                if kind == "barcodes":
                    # Deposit ships none; the cleaner synthesizes positional ids.
                    continue
                print(f"ERROR: missing {kind} for run_id={run_id}: {src}", file=sys.stderr)
                return 1
            if mode == "symlink":
                _stage_symlink(src, dst)
            else:
                _stage_copy(src, dst)

        staged += 1
        print(f"staged {args.study}/{run_id} -> {run_dir}")

    print(f"staged {staged} runs under {triplet_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
