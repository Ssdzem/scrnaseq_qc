import os
from pathlib import Path


REF_CFG = config.get("fastq_prepros", {}).get("celescope", {}).get("reference", {})

REF_DIR = os.path.normpath(os.path.expanduser(REF_CFG.get("ref_dir", "data/aux/ref/hs_ensembl_110")))
GENOME_NAME = REF_CFG.get("genome_name", "Homo_sapiens_ensembl_110_filtered")

FASTA_URL = REF_CFG.get(
    "fasta_url",
    "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
)
GTF_URL = REF_CFG.get(
    "gtf_url",
    "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz",
)


def _ref_paths():
    ref = Path(REF_DIR)
    return {
        "dir": str(ref),
        "fasta_gz": str(ref / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
        "gtf_gz": str(ref / "Homo_sapiens.GRCh38.110.gtf.gz"),
        "fasta": str(ref / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
        "gtf": str(ref / "Homo_sapiens.GRCh38.110.gtf"),
        "filtered_gtf": str(ref / "Homo_sapiens.GRCh38.110.filtered.gtf"),
        "mt_gene_list": str(ref / "mt_gene_list.txt"),
        # CeleScope mkref output typically includes a STAR index at {outdir}/star.
        "star_index": str(ref / "star"),
        "done": str(ref / ".celescope_ref_done"),
    }


P = _ref_paths()


rule celescope_ref_all:
    input:
        P["done"]


rule download_ensembl_fasta_gtf:
    output:
        fasta_gz=P["fasta_gz"],
        gtf_gz=P["gtf_gz"],
    params:
        ref_dir=P["dir"],
        fasta_url=FASTA_URL,
        gtf_url=GTF_URL,
    log:
        "out/logs/fastq_prepros/ref/download_ensembl_fasta_gtf.log",
    conda:
        "../env/fastq_prepros.yaml"
    shell:
        "python workflow/modules/fastq_prepros/scripts/ref_download_ensembl.py "
        "--ref-dir {params.ref_dir} "
        "--fasta-url {params.fasta_url} "
        "--gtf-url {params.gtf_url} "
        "> {log} 2>&1"


rule gunzip_ensembl:
    input:
        fasta_gz=P["fasta_gz"],
        gtf_gz=P["gtf_gz"],
    output:
        fasta=P["fasta"],
        gtf=P["gtf"],
    log:
        "out/logs/fastq_prepros/ref/gunzip_ensembl.log",
    conda:
        "../env/fastq_prepros.yaml"
    shell:
        "python workflow/modules/fastq_prepros/scripts/ref_gunzip.py "
        "--fasta-gz {input.fasta_gz} "
        "--gtf-gz {input.gtf_gz} "
        "--fasta {output.fasta} "
        "--gtf {output.gtf} "
        "> {log} 2>&1"


rule ensure_mt_gene_list:
    output:
        mt_gene_list=P["mt_gene_list"],
    log:
        "out/logs/fastq_prepros/ref/ensure_mt_gene_list.log",
    conda:
        "../env/fastq_prepros.yaml"
    shell:
        "python workflow/modules/fastq_prepros/scripts/ref_write_mt_gene_list.py "
        "--output {output.mt_gene_list} "
        "> {log} 2>&1"


rule celescope_mkgtf:
    input:
        gtf=P["gtf"],
    output:
        filtered_gtf=P["filtered_gtf"],
    log:
        "out/logs/fastq_prepros/ref/celescope_mkgtf.log",
    conda:
        "../env/fastq_prepros.yaml"
    shell:
        "python workflow/modules/fastq_prepros/scripts/ref_mkgtf.py "
        "--gtf {input.gtf} "
        "--filtered-gtf {output.filtered_gtf} "
        "> {log} 2>&1"


rule celescope_mkref:
    input:
        fasta=P["fasta"],
        filtered_gtf=P["filtered_gtf"],
        mt_gene_list=P["mt_gene_list"],
    output:
        done=touch(P["done"]),
    params:
        ref_dir=P["dir"],
        genome_name=GENOME_NAME,
    log:
        "out/logs/fastq_prepros/ref/celescope_mkref.log",
    conda:
        "../env/fastq_prepros.yaml"
    shell:
        "python workflow/modules/fastq_prepros/scripts/ref_mkref.py "
        "--outdir {params.ref_dir} "
        "--genome-name {params.genome_name} "
        "--fasta {input.fasta} "
        "--filtered-gtf {input.filtered_gtf} "
        "--mt-gene-list {input.mt_gene_list} "
        "> {log} 2>&1"
