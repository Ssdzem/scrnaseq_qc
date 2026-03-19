import argparse
import os
import subprocess
from pathlib import Path


def _rename_fastqs_in_place(fastqs: list[str]) -> None:
    """Normalize FASTQ names in-place to CeleScope-compatible *_1/_2 scheme.

    This repo stores FASTQs as *_f1.fastq.gz and *_r2.fastq.gz.
    CeleScope multi_rna expects *_1.fastq.gz and *_2.fastq.gz (or R1/R2 variants).
    """

    for fq in fastqs:
        fq_path = Path(fq)
        name = fq_path.name
        if name.endswith("_f1.fastq.gz"):
            target = fq_path.with_name(name.removesuffix("_f1.fastq.gz") + "_1.fastq.gz")
        elif name.endswith("_r2.fastq.gz"):
            target = fq_path.with_name(name.removesuffix("_r2.fastq.gz") + "_2.fastq.gz")
        else:
            continue

        if fq_path == target:
            continue
        if target.exists():
            # Assume already normalized.
            continue
        fq_path.rename(target)


def write_mapfile(mapfile_path: Path, sample: str, fastq_dir: Path) -> None:
    """Create CeleScope multi_rna mapfile.

    Format (tab-separated):
      fastq_prefix  fastq_dir  sample

    This project stores FASTQs as:
      out/<study>/fastqs/<donor>/<sample>/<run>/<run>_f1.fastq.gz
      out/<study>/fastqs/<donor>/<sample>/<run>/<run>_r2.fastq.gz

    We assume `fastq_prefix` is the common prefix without read suffix,
    e.g. HRR1782813 for HRR1782813_f1.fastq.gz / HRR1782813_r2.fastq.gz.
    """

    mapfile_path.write_text(f"{sample}\t{fastq_dir}\t{sample}\n")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument(
        "--fastqs",
        required=True,
        help="Comma-separated list of FASTQ paths (R1+R2).",
    )
    p.add_argument("--genome-dir", required=True, help="CeleScope mkref output directory")
    p.add_argument("--threads", type=int, required=True)
    p.add_argument("--mod", default="shell")
    p.add_argument("--solo-cell-filter", default="EmptyDrops_CR")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mapfile = outdir / "rna.mapfile"

    fastqs = [p for p in args.fastqs.split(",") if p.strip()]
    _rename_fastqs_in_place(fastqs)
    fastq_dir = Path(fastqs[0]).parent
    write_mapfile(mapfile, args.sample, fastq_dir)

    cmd = [
        "multi_rna",
        "--mapfile",
        str(mapfile),
        "--genomeDir",
        args.genome_dir,
        "--thread",
        str(args.threads),
        "--soloCellFilter",
        args.solo_cell_filter,
        "--mod",
        args.mod,
    ]

    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", str(args.threads))
    subprocess.run(cmd, check=True, env=env)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
