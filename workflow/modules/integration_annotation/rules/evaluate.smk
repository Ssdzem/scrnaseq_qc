rule ia_evaluate_all:
    input:
        annotated=f"{IA_ALL_OUTPUT_DIR}/all_annotated.h5ad",
        ia_config=IA_CONFIG_PATH,
        annotated_done=f"{IA_ALL_OUTPUT_DIR}/.annotated_done",
    output:
        metrics_json=f"{IA_ALL_OUTPUT_DIR}/evaluation_metrics.json",
        metrics_csv=f"{IA_ALL_OUTPUT_DIR}/evaluation_metrics.csv",
        done=touch(f"{IA_ALL_OUTPUT_DIR}/.evaluated_done"),
    log:
        _ia_log("all_studies")
    threads:
        _ia_evaluate_threads
    resources:
        mem_mb=_ia_evaluate_mem_mb,
    conda:
        "../env/integration_annotation_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        A1_P1_IA_EVAL_THREADS={threads} OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
        PYTHONNOUSERSITE=1 $CONDA_PREFIX/bin/python -s workflow/modules/integration_annotation/scripts/evaluate_integration.py \
          --scope all_studies \
          --input {input.annotated} \
          --config {input.ia_config} \
          --output-json {output.metrics_json} \
          --output-csv {output.metrics_csv} \
          --log {log}
        """


rule ia_evaluate_study:
    wildcard_constraints:
        study=r"(?!all_studies$)[^/]+"
    input:
        annotated=f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/{{study}}_annotated.h5ad",
        ia_config=IA_CONFIG_PATH,
        annotated_done=f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/.{{study}}_annotated_done",
    output:
        metrics_json=f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/{{study}}_evaluation_metrics.json",
        metrics_csv=f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/{{study}}_evaluation_metrics.csv",
        done=touch(f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/.{{study}}_evaluated_done"),
    params:
        study=lambda wc: wc.study,
    log:
        _ia_log("{study}")
    threads:
        _ia_evaluate_threads
    resources:
        mem_mb=_ia_evaluate_mem_mb,
    conda:
        "../env/integration_annotation_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        A1_P1_IA_EVAL_THREADS={threads} OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
        PYTHONNOUSERSITE=1 $CONDA_PREFIX/bin/python -s workflow/modules/integration_annotation/scripts/evaluate_integration.py \
          --scope study \
          --study {params.study} \
          --input {input.annotated} \
          --config {input.ia_config} \
          --output-json {output.metrics_json} \
          --output-csv {output.metrics_csv} \
          --log {log}
        """
