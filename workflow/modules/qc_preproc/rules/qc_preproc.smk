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


def _qc_log(study: str) -> str:
    return f"out/logs/qc/{study}.log"


STUDIES = _studies_to_process()


rule qc_all:
    input:
        expand("out/{study}/qc/.pretrim_pictures_done", study=STUDIES),
        expand("out/{study}/qc/.study_explore_done", study=STUDIES),
        "out/all_studies/.all_raw_done",
        "out/all_studies/.all_explore_done",
    output:
        touch(expand("out/{study}/qc/.qc_all_done", study=STUDIES))


rule qc_pretrim_pictures:
    input:
        h5ad=lambda wc: _study_raw_h5ad(wc.study),
        qc_config=QC_CONFIG_PATH,
    output:
        metrics="out/{study}/qc/pre_filter_metrics.csv",
        genes="out/{study}/qc/pictures/pre_filter_genes_per_cell.png",
        counts="out/{study}/qc/pictures/pre_filter_counts_per_cell.png",
        mito="out/{study}/qc/pictures/pre_filter_mito_pct.png",
        done=touch("out/{study}/qc/.pretrim_pictures_done"),
    params:
        study=lambda wc: wc.study,
    log:
        lambda wc: _qc_log(wc.study)
    conda:
        "../env/qc_preproc_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        : > {log}
        {{
          echo "[$(date -Is)] [qc_pretrim_pictures] study={params.study}"
          echo "[$(date -Is)] INPUT input_h5ad={input.h5ad}"
          echo "[$(date -Is)] INPUT qc_config={input.qc_config}"
          echo "[$(date -Is)] OUTPUT metrics={output.metrics}"
          echo "[$(date -Is)] OUTPUT genes={output.genes}"
          echo "[$(date -Is)] OUTPUT counts={output.counts}"
          echo "[$(date -Is)] OUTPUT mito={output.mito}"
          $CONDA_PREFIX/bin/python workflow/modules/qc_preproc/scripts/qc_figures.py \
          --study {params.study} \
          --input-h5ad {input.h5ad} \
          --qc-config {input.qc_config} \
          --output-metrics {output.metrics} \
          --output-genes-fig {output.genes} \
          --output-counts-fig {output.counts} \
          --output-mito-fig {output.mito}
        }} 2>&1 | tee -a {log}
        """


rule qc_study_explore:
    input:
        raw=lambda wc: _study_raw_h5ad(wc.study),
        qc_config=QC_CONFIG_PATH,
        pretrim_done="out/{study}/qc/.pretrim_pictures_done",
    output:
        explore_h5="out/{study}/qc/{study}_explore.h5ad",
        done=touch("out/{study}/qc/.study_explore_done"),
    params:
        study=lambda wc: wc.study,
    log:
        lambda wc: _qc_log(wc.study)
    conda:
        "../env/qc_preproc_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        {{
          echo "[$(date -Is)] [qc_study_explore] study={params.study}"
          echo "[$(date -Is)] INPUT raw={input.raw}"
          echo "[$(date -Is)] INPUT qc_config={input.qc_config}"
          echo "[$(date -Is)] INPUT pretrim_done={input.pretrim_done}"
          echo "[$(date -Is)] OUTPUT explore_h5={output.explore_h5}"
          $CONDA_PREFIX/bin/python workflow/modules/qc_preproc/scripts/explore_scanpy.py \
          --study {params.study} \
          --input {input.raw} \
          --qc-config {input.qc_config} \
          --output {output.explore_h5}
        }} 2>&1 | tee -a {log}
        """


rule qc_concat_all_raw:
    input:
        h5ads=expand("out/{study}/raw.h5ad", study=STUDIES),
    output:
        all_raw="out/all_studies/all_raw.h5ad",
        done=touch("out/all_studies/.all_raw_done"),
    log:
        lambda wc: _qc_log("all_studies")
    conda:
        "../env/qc_preproc_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        : > {log}
        {{
          echo "[$(date -Is)] [qc_concat_all_raw] study=all_studies"
          echo "[$(date -Is)] INPUT h5ads={input.h5ads}"
          echo "[$(date -Is)] OUTPUT all_raw={output.all_raw}"
          $CONDA_PREFIX/bin/python workflow/modules/qc_preproc/scripts/concat_raw.py \
          --inputs {input.h5ads} \
          --output {output.all_raw}
        }} 2>&1 | tee -a {log}
        """


rule qc_all_explore:
    input:
        all_raw="out/all_studies/all_raw.h5ad",
        qc_config=QC_CONFIG_PATH,
        all_raw_done="out/all_studies/.all_raw_done",
    output:
        all_explore="out/all_studies/all_explore.h5ad",
        done=touch("out/all_studies/.all_explore_done"),
    log:
        lambda wc: _qc_log("all_studies")
    conda:
        "../env/qc_preproc_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        {{
          echo "[$(date -Is)] [qc_all_explore] study=all_studies"
          echo "[$(date -Is)] INPUT all_raw={input.all_raw}"
          echo "[$(date -Is)] INPUT qc_config={input.qc_config}"
          echo "[$(date -Is)] INPUT all_raw_done={input.all_raw_done}"
          echo "[$(date -Is)] OUTPUT all_explore={output.all_explore}"
          $CONDA_PREFIX/bin/python workflow/modules/qc_preproc/scripts/explore_scanpy.py \
          --study all_studies \
          --input {input.all_raw} \
          --qc-config {input.qc_config} \
          --output {output.all_explore}
        }} 2>&1 | tee -a {log}
        """
