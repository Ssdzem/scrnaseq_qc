import argparse
import subprocess


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True)
    p.add_argument("--genome-name", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--filtered-gtf", required=True)
    p.add_argument("--mt-gene-list", required=True)
    args = p.parse_args()

    cmd = [
        "celescope",
        "rna",
        "mkref",
        "--outdir",
        args.outdir,
        "--genome_name",
        args.genome_name,
        "--fasta",
        args.fasta,
        "--gtf",
        args.filtered_gtf,
        "--mt_gene_list",
        args.mt_gene_list,
    ]
    subprocess.run(cmd, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

