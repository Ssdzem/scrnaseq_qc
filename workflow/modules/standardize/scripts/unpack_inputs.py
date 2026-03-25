from __future__ import annotations

import argparse
import gzip
import shutil
import tarfile
from pathlib import Path
from typing import Iterable

from rtlog import log_kv, log_line


def _is_gzip(path: Path) -> bool:
    """Return True if `path` looks like gzip by magic bytes."""

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
        for member in members:
            target = dst_dir / member.name
            if not is_within_directory(dst_dir, target):
                raise SystemExit(f"Refusing to extract outside destination: {member.name}")
        tf.extractall(dst_dir)

    log_line(log, f"Extracted tar: {tar_path} -> {dst_dir}")


def _maybe_create_top_level_symlink(*, tar_path: Path, log: TextIO) -> None:
    """Create a convenience symlink for single-folder tarball extractions.

    Many GEO archives are like `<name>.tar.gz` containing `processed_files/...`.
    In that case extraction yields `<parent>/<name>/processed_files/...`.
    If so, we create `<parent>/processed_files -> <parent>/<name>/processed_files`.
    """

    name = tar_path.name
    archive_stem = name
    for suffix in (".tar.gz", ".tgz", ".tar"):
        if archive_stem.endswith(suffix):
            archive_stem = archive_stem[: -len(suffix)]
            break

    nested_root = tar_path.parent / archive_stem
    if not (nested_root.exists() and nested_root.is_dir()):
        return

    top_level_dirs = [p for p in nested_root.iterdir() if p.is_dir()]
    if len(top_level_dirs) != 1:
        return

    top = top_level_dirs[0]
    link_path = tar_path.parent / top.name
    if link_path.exists():
        return

    try:
        link_path.symlink_to(top, target_is_directory=True)
        log_line(log, f"Symlink: {link_path} -> {top}")
    except OSError as e:
        log_line(log, f"WARNING: Failed to create symlink {link_path}: {e}")


def _gunzip_file(*, gz_path: Path, dst_path: Path, log: TextIO) -> None:
    dst_path.parent.mkdir(parents=True, exist_ok=True)

    if dst_path.exists() and dst_path.stat().st_size > 0:
        log_line(log, f"Skip gunzip (exists): {dst_path}")
        return

    try:
        with gzip.open(gz_path, "rb") as fin, dst_path.open("wb") as fout:
            shutil.copyfileobj(fin, fout)
    except gzip.BadGzipFile:
        if dst_path.exists():
            dst_path.unlink()
        log_line(log, f"WARNING: Not a gzipped file, skipping gunzip: {gz_path}")
        return

    log_line(log, f"Gunzip: {gz_path} -> {dst_path}")


def _iter_archives(root: Path) -> Iterable[Path]:
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue

        if p.suffix in {".tar", ".gz", ".tgz"}:
            yield p
            continue

        if p.name.endswith(".tar.gz"):
            yield p


def _unpack_in_dir(root: Path, *, log: TextIO) -> None:
    if not root.exists():
        log_line(log, f"Skip missing dir: {root}")
        return

    for path in _iter_archives(root):
        name = path.name

        # Guardrail: never touch FASTQs.
        if path.parent.name == "fastqs" or "fastqs" in path.parts:
            continue

        if name.endswith((".tar", ".tar.gz", ".tgz")):
            _safe_extract_tar(tar_path=path, dst_dir=path.parent, log=log)
            _maybe_create_top_level_symlink(tar_path=path, log=log)
            continue

        if path.suffix == ".gz":
            # Some providers ship files ending in .gz that are not gzipped.
            if not _is_gzip(path):
                log_line(log, f"Skip non-gzip .gz file: {path}")
                continue
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
        log_kv(log, "unpack.start", study=args.study, outdir=args.outdir)
        _unpack_in_dir(raw_dir, log=log)
        _unpack_in_dir(qconly_dir, log=log)
        log_kv(log, "unpack.complete", study=args.study)


if __name__ == "__main__":
    main()
