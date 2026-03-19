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


def _symlink_force(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    dest.symlink_to(src)


def stage_fastqs(outdir: Path, fastqs: list[str], sample: str) -> tuple[list[str], Path, str]:
    """Stage FASTQs into a CeleScope-friendly naming scheme.

    CeleScope multi expects reads named like *_1/*_2 or *R1/*R2. Many of our
    inputs are named *_f1/*_r2, which CeleScope does not recognize.

    This creates symlinks under <outdir>/fastq_links/<sample>/ as:
      <sample>_1.fastq.gz -> original R1
      <sample>_2.fastq.gz -> original R2

    Returns: (staged_fastqs, staged_dir, library_prefix)
    """

    r1 = next((Path(p) for p in fastqs if Path(p).name.endswith("_f1.fastq.gz")), None)
    r2 = next((Path(p) for p in fastqs if Path(p).name.endswith("_r2.fastq.gz")), None)
    if r1 is None or r2 is None:
        # Fallback: do not stage; keep original behavior.
        return fastqs, Path(fastqs[0]).parent, sample

    staged_dir = outdir / "fastq_links" / sample
    staged_r1 = staged_dir / f"{sample}_1.fastq.gz"
    staged_r2 = staged_dir / f"{sample}_2.fastq.gz"
    _symlink_force(r1.resolve(), staged_r1)
    _symlink_force(r2.resolve(), staged_r2)
    return [str(staged_r1), str(staged_r2)], staged_dir, sample


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
    fastqs, staged_dir, library_prefix = stage_fastqs(outdir, fastqs, args.sample)
    write_mapfile(mapfile, fastqs, library_prefix)

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
