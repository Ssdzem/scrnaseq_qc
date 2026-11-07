rule ia_integrate_all:
    input:
        pre_integration=IA_ALL_PRE_INTEGRATION_H5AD,
        ia_config=IA_CONFIG_PATH,
    output:
        integrated=f"{IA_ALL_OUTPUT_DIR}/all_integrated.h5ad",
        report=f"{IA_ALL_OUTPUT_DIR}/integration_report.json",
        done=touch(f"{IA_ALL_OUTPUT_DIR}/.integrated_done"),
    params:
        integration_method=_ia_integration_method,
    log:
        _ia_log("all_studies")
    threads:
        _ia_integrate_threads
    resources:
        mem_mb=_ia_integrate_mem_mb,
        gpu=_ia_integrate_gpu,
    conda:
        _ia_integrate_env()
    shell:
        r"""
        mkdir -p $(dirname {log})
        SCRNASEQ_QC_IA_NUM_THREADS={threads} OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
        PYTHONNOUSERSITE=1 $CONDA_PREFIX/bin/python -s workflow/modules/integration_annotation/scripts/integrate_embeddings.py \
          --scope all_studies \
          --input {input.pre_integration} \
          --config {input.ia_config} \
          --output {output.integrated} \
          --report {output.report} \
          --log {log}
        """


rule ia_integrate_study:
    wildcard_constraints:
        study=r"(?!all_studies$)[^/]+"
    input:
        pre_integration=IA_STUDY_PRE_INTEGRATION_H5AD_TEMPLATE,
        ia_config=IA_CONFIG_PATH,
    output:
        integrated=f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/{{study}}_integrated.h5ad",
        report=f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/{{study}}_integration_report.json",
        done=touch(f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/.{{study}}_integrated_done"),
    params:
        study=lambda wc: wc.study,
        integration_method=_ia_integration_method,
    log:
        _ia_log("{study}")
    threads:
        _ia_integrate_threads
    resources:
        mem_mb=_ia_integrate_mem_mb,
        gpu=_ia_integrate_gpu,
    conda:
        _ia_integrate_env()
    shell:
        r"""
        mkdir -p $(dirname {log})
        SCRNASEQ_QC_IA_NUM_THREADS={threads} OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
        PYTHONNOUSERSITE=1 $CONDA_PREFIX/bin/python -s workflow/modules/integration_annotation/scripts/integrate_embeddings.py \
          --scope study \
          --study {params.study} \
          --input {input.pre_integration} \
          --config {input.ia_config} \
          --output {output.integrated} \
          --report {output.report} \
          --log {log}
        """
