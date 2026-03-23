import os
from pathlib import Path


QC_CONFIG_PATH = str(Path("config/qc.yaml"))


def _studies_to_process():
    studies = config.get("inputs", {}).get("studies", [])
    if studies:
        return list(studies)

    meta = config.get("inputs", {}).get("metadata_csv") or "config/liver_scRNAseq_metadata_everything.csv"
    try:
        import csv
        import re

        def _is_missing(v: str) -> bool:
            return (v or "").strip() == "" or (v or "").strip().lower() == "na"

        with open(meta, newline="", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle)
            studies_set = set()
            for row in reader:
                study = (row.get("study_accession") or "").strip()
                if not _is_missing(study):
                    studies_set.add(re.sub(r"\s+", "", study))
        return sorted(studies_set)
    except FileNotFoundError:
        return []


def _study_raw_h5ad(study: str) -> str:
    return f"out/{study}/raw.h5ad"


STUDIES = _studies_to_process()


rule qc_all:
    input:
        expand("out/{study}/qc/qc.h5ad", study=STUDIES),
        expand("out/{study}/qc/pre_filter_metrics.csv", study=STUDIES),
        expand(
            "out/{study}/qc/figures/{fig}",
            study=STUDIES,
            fig=[
                "pre_filter_genes_per_cell.png",
                "pre_filter_counts_per_cell.png",
                "pre_filter_mito_pct.png",
            ],
        )


rule qc_figures:
    input:
        h5ad=lambda wc: _study_raw_h5ad(wc.study),
        qc_config=QC_CONFIG_PATH,
    output:
        metrics="out/{study}/qc/pre_filter_metrics.csv",
        genes="out/{study}/qc/figures/pre_filter_genes_per_cell.png",
        counts="out/{study}/qc/figures/pre_filter_counts_per_cell.png",
        mito="out/{study}/qc/figures/pre_filter_mito_pct.png",
    params:
        study=lambda wc: wc.study,
    log:
        "out/logs/qc_preproc/qc_figures/{study}.log"
    conda:
        "../env/qc_preproc_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        $CONDA_PREFIX/bin/python workflow/modules/qc_preproc/scripts/qc_figures.py \
          --study {params.study} \
          --input-h5ad {input.h5ad} \
          --qc-config {input.qc_config} \
          --output-metrics {output.metrics} \
          --output-genes-fig {output.genes} \
          --output-counts-fig {output.counts} \
          --output-mito-fig {output.mito} \
          --log {log}
        """


rule qc_scanpy:
    input:
        raw=lambda wc: _study_raw_h5ad(wc.study),
        qc_config=QC_CONFIG_PATH,
        preqc="out/{study}/qc/pre_filter_metrics.csv",
    output:
        qc_h5="out/{study}/qc/qc.h5ad",
    params:
        study=lambda wc: wc.study,
    log:
        "out/logs/qc_preproc/qc_scanpy/{study}.log"
    conda:
        "../env/qc_preproc_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        $CONDA_PREFIX/bin/python workflow/modules/qc_preproc/scripts/qc_scanpy.py \
          --study {params.study} \
          --sample {params.study} \
          --input {input.raw} \
          --qc-config {input.qc_config} \
          --output {output.qc_h5} \
          --log {log}
        """
