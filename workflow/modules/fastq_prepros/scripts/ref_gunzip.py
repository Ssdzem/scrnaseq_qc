import argparse
import gzip
import shutil
from pathlib import Path


def gunzip(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and dst.stat().st_size > 0:
        return
    with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--fasta-gz", required=True)
    p.add_argument("--gtf-gz", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--gtf", required=True)
    args = p.parse_args()

    gunzip(Path(args.fasta_gz), Path(args.fasta))
    gunzip(Path(args.gtf_gz), Path(args.gtf))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

