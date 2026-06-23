import argparse
from pathlib import Path

import urllib.request


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--ref-dir", required=True)
    p.add_argument("--fasta-url", required=True)
    p.add_argument("--gtf-url", required=True)
    args = p.parse_args()

    ref = Path(args.ref_dir)
    download(args.fasta_url, ref / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
    download(args.gtf_url, ref / "Homo_sapiens.GRCh38.110.gtf.gz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

