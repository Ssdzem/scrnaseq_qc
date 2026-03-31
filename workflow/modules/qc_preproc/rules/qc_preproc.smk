import os
from pathlib import Path
import yaml


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


def _study_raw_metadata_h5ad(study: str) -> str:
    return f"out/{study}/raw_metadata.h5ad"


def _qc_log(study: str) -> str:
    return f"out/logs/qc/{study}.log"


STUDIES = _studies_to_process()

def _study_explore_mode() -> str:
    qc_cfg = config.get("qc", {})
    if not qc_cfg:
        try:
            with open(QC_CONFIG_PATH, "r", encoding="utf-8") as handle:
                qc_cfg = (yaml.safe_load(handle) or {}).get("qc", {})
        except FileNotFoundError:
            qc_cfg = {}

    setting = qc_cfg.get("study_explore_with_integration", True)
    if isinstance(setting, str):
        setting = setting.strip().lower() in {"1", "true", "yes", "y", "on"}

    if bool(setting):
        return "integration"
    return "no_integration"


rule qc_all:
    input:
        qc_config=QC_CONFIG_PATH,
        raw_metadata_done=expand("out/{study}/qc/.raw_metadata_done", study=STUDIES),
        pretrim_done=expand("out/{study}/qc/.pretrim_pictures_done", study=STUDIES),
        study_explore_done=expand("out/{study}/qc/.study_explore_done", study=STUDIES),
        all_raw_done="out/all_studies/.all_raw_done",
        all_explore_done="out/all_studies/.all_explore_done",
    output:
        touch("out/all_studies/.qc_all_done")


rule qc_attach_metadata:
    input:
        raw=lambda wc: _study_raw_h5ad(wc.study),
        qc_config=QC_CONFIG_PATH,
    output:
        raw_metadata=lambda wc: _study_raw_metadata_h5ad(wc.study),
        done=touch("out/{study}/qc/.raw_metadata_done"),
    params:
        study=lambda wc: wc.study,
    log:
        _qc_log("{study}")
    conda:
        "../env/qc_preproc_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        {{
          echo "[$(date -Is)] [qc_attach_metadata] study={params.study}"
          echo "[$(date -Is)] INPUT raw={input.raw}"
          echo "[$(date -Is)] INPUT qc_config={input.qc_config}"
          echo "[$(date -Is)] OUTPUT raw_metadata={output.raw_metadata}"
          $CONDA_PREFIX/bin/python workflow/modules/qc_preproc/scripts/attach_metadata.py \
          --study {params.study} \
          --input {input.raw} \
          --qc-config {input.qc_config} \
          --output {output.raw_metadata}
        }} 2>&1 | tee -a {log}
        """


rule qc_pretrim_pictures:
    input:
        h5ad=lambda wc: _study_raw_metadata_h5ad(wc.study),
        qc_config=QC_CONFIG_PATH,
        raw_metadata_done="out/{study}/qc/.raw_metadata_done",
    output:
        metrics="out/{study}/qc/pre_filter_metrics.csv",
        genes="out/{study}/qc/pictures/pre_filter_genes_per_cell.png",
        counts="out/{study}/qc/pictures/pre_filter_counts_per_cell.png",
        mito="out/{study}/qc/pictures/pre_filter_mito_pct.png",
        done=touch("out/{study}/qc/.pretrim_pictures_done"),
    params:
        study=lambda wc: wc.study,
    log:
        _qc_log("{study}")
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
          echo "[$(date -Is)] INPUT raw_metadata_done={input.raw_metadata_done}"
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
        pre_integration_h5="out/{study}/qc/{study}_pre_integration.h5ad",
        qc_config=QC_CONFIG_PATH,
        pre_integration_done="out/{study}/qc/.study_pre_integration_done",
    output:
        explore_h5="out/{study}/qc/{study}_explore.h5ad",
        done=touch("out/{study}/qc/.study_explore_done"),
    params:
        study=lambda wc: wc.study,
        mode=lambda wc: _study_explore_mode(),
    log:
        _qc_log("{study}")
    conda:
        "../env/qc_preproc_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        {{
          echo "[$(date -Is)] [qc_study_explore] study={params.study}"
          echo "[$(date -Is)] INPUT pre_integration_h5={input.pre_integration_h5}"
          echo "[$(date -Is)] INPUT qc_config={input.qc_config}"
          echo "[$(date -Is)] INPUT pre_integration_done={input.pre_integration_done}"
          echo "[$(date -Is)] PARAM mode={params.mode}"
          echo "[$(date -Is)] OUTPUT explore_h5={output.explore_h5}"
          $CONDA_PREFIX/bin/python workflow/modules/qc_preproc/scripts/explore_scanpy.py \
          --mode {params.mode} \
          --study {params.study} \
          --input {input.pre_integration_h5} \
          --qc-config {input.qc_config} \
          --output {output.explore_h5}
        }} 2>&1 | tee -a {log}
        """


rule qc_study_pre_integration:
    input:
        raw=lambda wc: _study_raw_metadata_h5ad(wc.study),
        qc_config=QC_CONFIG_PATH,
        pretrim_done="out/{study}/qc/.pretrim_pictures_done",
    output:
        pre_integration_h5="out/{study}/qc/{study}_pre_integration.h5ad",
        done=touch("out/{study}/qc/.study_pre_integration_done"),
    params:
        study=lambda wc: wc.study,
    log:
        _qc_log("{study}")
    conda:
        "../env/qc_preproc_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        {{
          echo "[$(date -Is)] [qc_study_pre_integration] study={params.study}"
          echo "[$(date -Is)] INPUT raw={input.raw}"
          echo "[$(date -Is)] INPUT qc_config={input.qc_config}"
          echo "[$(date -Is)] INPUT pretrim_done={input.pretrim_done}"
          echo "[$(date -Is)] OUTPUT pre_integration_h5={output.pre_integration_h5}"
          $CONDA_PREFIX/bin/python workflow/modules/qc_preproc/scripts/explore_scanpy.py \
          --mode pre \
          --study {params.study} \
          --input {input.raw} \
          --qc-config {input.qc_config} \
          --output {output.pre_integration_h5}
        }} 2>&1 | tee -a {log}
        """


rule qc_concat_all_raw:
    input:
        qc_config=QC_CONFIG_PATH,
        h5ads=expand("out/{study}/raw_metadata.h5ad", study=STUDIES),
    output:
        all_raw="out/all_studies/all_raw.h5ad",
        done=touch("out/all_studies/.all_raw_done"),
    log:
        _qc_log("all_studies")
    conda:
        "../env/qc_preproc_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        : > {log}
        {{
          echo "[$(date -Is)] [qc_concat_all_raw] study=all_studies"
          echo "[$(date -Is)] INPUT qc_config={input.qc_config}"
          echo "[$(date -Is)] INPUT h5ads={input.h5ads}"
          echo "[$(date -Is)] OUTPUT all_raw={output.all_raw}"
          $CONDA_PREFIX/bin/python workflow/modules/qc_preproc/scripts/concat_raw.py \
          --inputs {input.h5ads} \
          --output {output.all_raw}
        }} 2>&1 | tee -a {log}
        """


rule qc_all_explore:
    input:
        all_pre_integration="out/all_studies/all_pre_integration.h5ad",
        qc_config=QC_CONFIG_PATH,
        all_pre_integration_done="out/all_studies/.all_pre_integration_done",
    output:
        all_explore="out/all_studies/all_explore.h5ad",
        done=touch("out/all_studies/.all_explore_done"),
    log:
        _qc_log("all_studies")
    conda:
        "../env/qc_preproc_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        {{
          echo "[$(date -Is)] [qc_all_explore] study=all_studies"
          echo "[$(date -Is)] INPUT all_pre_integration={input.all_pre_integration}"
          echo "[$(date -Is)] INPUT qc_config={input.qc_config}"
          echo "[$(date -Is)] INPUT all_pre_integration_done={input.all_pre_integration_done}"
          echo "[$(date -Is)] OUTPUT all_explore={output.all_explore}"
          $CONDA_PREFIX/bin/python workflow/modules/qc_preproc/scripts/explore_scanpy.py \
          --mode integration \
          --study all_studies \
          --input {input.all_pre_integration} \
          --qc-config {input.qc_config} \
          --output {output.all_explore}
        }} 2>&1 | tee -a {log}
        """


rule qc_all_pre_integration:
    input:
        all_raw="out/all_studies/all_raw.h5ad",
        qc_config=QC_CONFIG_PATH,
        all_raw_done="out/all_studies/.all_raw_done",
    output:
        all_pre_integration="out/all_studies/all_pre_integration.h5ad",
        done=touch("out/all_studies/.all_pre_integration_done"),
    log:
        _qc_log("all_studies")
    conda:
        "../env/qc_preproc_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        {{
          echo "[$(date -Is)] [qc_all_pre_integration] study=all_studies"
          echo "[$(date -Is)] INPUT all_raw={input.all_raw}"
          echo "[$(date -Is)] INPUT qc_config={input.qc_config}"
          echo "[$(date -Is)] INPUT all_raw_done={input.all_raw_done}"
          echo "[$(date -Is)] OUTPUT all_pre_integration={output.all_pre_integration}"
          $CONDA_PREFIX/bin/python workflow/modules/qc_preproc/scripts/explore_scanpy.py \
          --mode pre \
          --study all_studies \
          --input {input.all_raw} \
          --qc-config {input.qc_config} \
          --output {output.all_pre_integration}
        }} 2>&1 | tee -a {log}
        """
