"""Extract a GEO ``*_RAW.tar`` archive into the layout a standardize cleaner
expects. Reproducible, idempotent prep step run as part of the fetch module
(rule ``fetch_extract``) — it replaces the old, manually-run
``tools/extract_validation_raw.py``.

Two archive shapes are supported, selected per study via the fetch config
``extract:`` map (``extract: {<STUDY>: triplet|csv}``):

  triplet  per-sample 10x triplets -> one dir per GSM, ready for
           ``count_matrix_10x_multi_runs`` via ``triplet_layout``:
             <out_root>/<study>/raw/triplets/<GSM>/{matrix.mtx.gz,
                                                    features.tsv.gz,
                                                    barcodes.tsv.gz}
           Handles both flat GEO members (``<GSM>_<label>_matrix.mtx.gz``) and
           nested per-sample archives (``<GSM>_*.tar.gz`` holding a triplet).

  csv      per-sample count CSVs -> renamed to the ``csv_handler`` pattern:
             <out_root>/<study>/raw/<GSM>_<label>-raw_counts.csv.gz

Usage (invoked by the fetch_extract rule):
  python extract_raw_archive.py --study STUDY_A --format triplet --out-root out
"""

from __future__ import annotations

import argparse
import gzip
import io
import re
import shutil
import tarfile
from pathlib import Path

_GSM = re.compile(r"^(GSM\d+)_(.+)$")


def _members(archive: Path):
    if not archive.exists():
        raise SystemExit(
            f"missing archive: {archive} "
            f"(run `fetch_data.py download-study {archive.parent.parent.name}` first)"
        )
    with tarfile.open(archive) as tf:
        for m in tf.getmembers():
            if m.isfile():
                yield tf, m


def _canonical(basename: str) -> str | None:
    """Map a triplet member basename (gzipped or plain) to the canonical name.

    Handles both flat GEO files (``<label>_matrix.mtx.gz``) and nested-archive
    files (``matrix.mtx``) via suffix matching.
    """
    for stem, canon in (("matrix.mtx", "matrix.mtx.gz"),
                        ("features.tsv", "features.tsv.gz"),
                        ("genes.tsv", "features.tsv.gz"),
                        ("barcodes.tsv", "barcodes.tsv.gz")):
        if basename.endswith(stem) or basename.endswith(stem + ".gz"):
            return canon
    return None


def _write_gz(data: bytes, dst: Path, already_gz: bool) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if already_gz:
        dst.write_bytes(data)
    else:
        with gzip.open(dst, "wb") as fh:
            fh.write(data)


def _extract_nested(tf: tarfile.TarFile, m: tarfile.TarInfo, out_dir: Path) -> int:
    """A per-sample <GSM>_*.tar.gz holding an uncompressed 10x triplet."""
    n = 0
    with tf.extractfile(m) as raw, tarfile.open(fileobj=io.BytesIO(raw.read())) as inner:
        for im in inner.getmembers():
            if not im.isfile():
                continue
            canon = _canonical(Path(im.name).name)
            if canon is None:
                continue
            data = inner.extractfile(im).read()
            _write_gz(data, out_dir / canon, already_gz=im.name.endswith(".gz"))
            n += 1
    return n


def _extract_triplets(archive: Path, raw_root: Path) -> int:
    dest_root = raw_root / "triplets"
    n = 0
    for tf, m in _members(archive):
        name = Path(m.name).name
        mo = _GSM.match(name)
        if not mo:
            continue
        gsm, rest = mo.group(1), mo.group(2)
        out_dir = dest_root / gsm
        if rest.endswith(".tar.gz"):  # nested per-sample archive
            n += _extract_nested(tf, m, out_dir)
            continue
        canon = _canonical(rest)
        if canon is None:
            continue
        out_dir.mkdir(parents=True, exist_ok=True)
        with tf.extractfile(m) as src, open(out_dir / canon, "wb") as dst:
            shutil.copyfileobj(src, dst)
        n += 1
    n_samples = len(list(dest_root.iterdir())) if dest_root.exists() else 0
    print(f"  wrote {n} triplet files across {n_samples} samples -> {dest_root}")
    return n


def _extract_csvs(archive: Path, raw_root: Path) -> int:
    raw_root.mkdir(parents=True, exist_ok=True)
    n = 0
    for tf, m in _members(archive):
        name = Path(m.name).name
        mo = _GSM.match(name)
        if not mo:
            continue
        gsm, rest = mo.group(1), mo.group(2)
        # GSM5172568_S1.count.csv.gz -> label "S1"
        label = re.sub(r"\.count\.csv\.gz$|\.csv\.gz$|\.csv$", "", rest)
        out_name = f"{gsm}_{label}-raw_counts.csv.gz"
        with tf.extractfile(m) as src, open(raw_root / out_name, "wb") as dst:
            shutil.copyfileobj(src, dst)
        n += 1
    print(f"  wrote {n} renamed count CSVs -> {raw_root}")
    return n


_EXTRACTORS = {
    "triplet": _extract_triplets,
    "csv": _extract_csvs,
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--study", required=True, help="study accession")
    ap.add_argument("--format", required=True, choices=sorted(_EXTRACTORS),
                    help="archive layout (drives which cleaner the output feeds)")
    ap.add_argument("--out-root", default="out",
                    help="fetch output root (per-study dirs live under here)")
    args = ap.parse_args()

    raw_root = Path(args.out_root) / args.study / "raw"
    archive = raw_root / f"{args.study}_RAW.tar"
    print(f"Extracting {archive} (format={args.format}):")
    _EXTRACTORS[args.format](archive, raw_root)
    print("done.")


if __name__ == "__main__":
    main()
