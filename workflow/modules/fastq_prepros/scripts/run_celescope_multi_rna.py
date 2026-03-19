import argparse
import os
import subprocess
from pathlib import Path


def write_mapfile(mapfile_path: Path, fastqs: list[str], sample: str) -> None:
    """Create CeleScope multi_rna mapfile.

    Format (tab-separated):
      fastq_prefix  fastq_dir  sample

    This project stores FASTQs as:
      out/<study>/fastqs/<donor>/<sample>/<run>/<run>_f1.fastq.gz
      out/<study>/fastqs/<donor>/<sample>/<run>/<run>_r2.fastq.gz

    We assume `fastq_prefix` is the common prefix without read suffix,
    e.g. HRR1782813 for HRR1782813_f1.fastq.gz / HRR1782813_r2.fastq.gz.
    """

    # Group by directory and derive a prefix compatible with CeleScope's
    # expected R1 patterns (e.g. <prefix>*_1.fastq.gz, <prefix>*R1_001.fastq.gz).
    per_dir: dict[Path, set[str]] = {}
    for fq in fastqs:
        fq_path = Path(fq)
        per_dir.setdefault(fq_path.parent, set())
        name = fq_path.name
        # Common forms in this repo:
        #   <run>_f1.fastq.gz / <run>_r2.fastq.gz
        # CeleScope expects prefixes that match *_1/*_2 or *R1/*R2.
        if name.endswith("_f1.fastq.gz"):
            prefix = name.removesuffix("_f1.fastq.gz") + "_"
        elif name.endswith("_r2.fastq.gz"):
            prefix = name.removesuffix("_r2.fastq.gz") + "_"
        else:
            prefix = name.split(".", 1)[0]
        per_dir[fq_path.parent].add(prefix)

    lines: list[str] = []
    for fq_dir in sorted(per_dir):
        for prefix in sorted(per_dir[fq_dir]):
            lines.append(f"{prefix}\t{fq_dir}\t{sample}")

    mapfile_path.write_text("\n".join(lines) + "\n")


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
    write_mapfile(mapfile, fastqs, args.sample)

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
