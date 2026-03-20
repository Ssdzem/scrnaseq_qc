from __future__ import annotations

import argparse
import gzip
import shutil
import tarfile
from pathlib import Path
from typing import Iterable, TextIO


def _log(msg: str, log: TextIO) -> None:
    log.write(msg.rstrip() + "\n")


def _is_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def _safe_extract_tar(*, tar_path: Path, dst_dir: Path, log: TextIO) -> None:
    dst_dir.mkdir(parents=True, exist_ok=True)

    def is_within_directory(directory: Path, target: Path) -> bool:
        directory_resolved = directory.resolve()
        target_resolved = target.resolve()
        return str(target_resolved).startswith(str(directory_resolved) + "/")

    with tarfile.open(tar_path, "r:*") as tf:
        members = tf.getmembers()
        for m in members:
            target = dst_dir / m.name
            if not is_within_directory(dst_dir, target):
                raise SystemExit(f"Refusing to extract outside destination: {m.name}")

        tf.extractall(dst_dir)

    _log(f"Extracted tar: {tar_path} -> {dst_dir}", log)


def _gunzip_file(*, gz_path: Path, dst_path: Path, log: TextIO) -> None:
    dst_path.parent.mkdir(parents=True, exist_ok=True)

    # Skip if already done and non-empty.
    if dst_path.exists() and dst_path.stat().st_size > 0:
        _log(f"Skip gunzip (exists): {dst_path}", log)
        return

    with gzip.open(gz_path, "rb") as fin, dst_path.open("wb") as fout:
        shutil.copyfileobj(fin, fout)

    _log(f"Gunzip: {gz_path} -> {dst_path}", log)


def _iter_archives(root: Path) -> Iterable[Path]:
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        if p.suffix in {".tar", ".gz", ".tgz"}:
            yield p
        elif p.name.endswith(".tar.gz"):
            yield p


def _unpack_in_dir(root: Path, *, log: TextIO) -> None:
    if not root.exists():
        _log(f"Skip missing dir: {root}", log)
        return

    for path in _iter_archives(root):
        name = path.name

        # Do not touch fastqs: caller should only pass raw/qconly.
        if path.parent.name == "fastqs" or "fastqs" in path.parts:
            continue

        if name.endswith(".tar") or name.endswith(".tar.gz") or name.endswith(".tgz"):
            # Extract next to the tarball, into the same directory.
            _safe_extract_tar(tar_path=path, dst_dir=path.parent, log=log)
            continue

        # Plain .gz: gunzip if it's actually gzipped.
        if path.suffix == ".gz" and _is_gzip(path):
            out_path = path.with_suffix("")
            _gunzip_file(gz_path=path, dst_path=out_path, log=log)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Unpack study inputs in out/{study}/raw and out/{study}/qconly"
    )
    parser.add_argument("--study", required=True)
    parser.add_argument("--outdir", default="out", help="Pipeline out/ directory")
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    out_root = Path(args.outdir)
    study_root = out_root / args.study
    raw_dir = study_root / "raw"
    qconly_dir = study_root / "qconly"
    log_path = Path(args.log)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("a", encoding="utf-8") as log:
        _log(f"Unpacking inputs for study={args.study}", log)
        _unpack_in_dir(raw_dir, log=log)
        _unpack_in_dir(qconly_dir, log=log)


if __name__ == "__main__":
    main()

