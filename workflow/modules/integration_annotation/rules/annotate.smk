rule ia_annotate_all:
    input:
        integrated=f"{IA_ALL_OUTPUT_DIR}/all_integrated.h5ad",
        ia_config=IA_CONFIG_PATH,
        integrated_done=f"{IA_ALL_OUTPUT_DIR}/.integrated_done",
        reference=_ia_reference_targets,
    output:
        annotated=f"{IA_ALL_OUTPUT_DIR}/all_annotated.h5ad",
        summary=f"{IA_ALL_OUTPUT_DIR}/annotation_summary.csv",
        done=touch(f"{IA_ALL_OUTPUT_DIR}/.annotated_done"),
    params:
        reference_h5ad=IA_REFERENCE_H5AD,
    log:
        _ia_log("all_studies")
    threads:
        _ia_annotate_threads
    resources:
        mem_mb=_ia_annotate_mem_mb,
        gpu=_ia_annotate_gpu,
    conda:
        _ia_annotate_env()
    shell:
        r"""
        mkdir -p $(dirname {log})
        OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
        PYTHONNOUSERSITE=1 $CONDA_PREFIX/bin/python -s workflow/modules/integration_annotation/scripts/annotate_cells.py \
          --scope all_studies \
          --input {input.integrated} \
          --config {input.ia_config} \
          --reference-h5ad {params.reference_h5ad} \
          --output {output.annotated} \
          --summary {output.summary} \
          --log {log}
        """


rule ia_annotate_study:
    wildcard_constraints:
        study=r"(?!all_studies$)[^/]+"
    input:
        integrated=f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/{{study}}_integrated.h5ad",
        ia_config=IA_CONFIG_PATH,
        integrated_done=f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/.{{study}}_integrated_done",
        reference=_ia_reference_targets,
    output:
        annotated=f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/{{study}}_annotated.h5ad",
        summary=f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/{{study}}_annotation_summary.csv",
        done=touch(f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/.{{study}}_annotated_done"),
    params:
        study=lambda wc: wc.study,
        reference_h5ad=IA_REFERENCE_H5AD,
    log:
        _ia_log("{study}")
    threads:
        _ia_annotate_threads
    resources:
        mem_mb=_ia_annotate_mem_mb,
        gpu=_ia_annotate_gpu,
    conda:
        _ia_annotate_env()
    shell:
        r"""
        mkdir -p $(dirname {log})
        OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
        PYTHONNOUSERSITE=1 $CONDA_PREFIX/bin/python -s workflow/modules/integration_annotation/scripts/annotate_cells.py \
          --scope study \
          --study {params.study} \
          --input {input.integrated} \
          --config {input.ia_config} \
          --reference-h5ad {params.reference_h5ad} \
          --output {output.annotated} \
          --summary {output.summary} \
          --log {log}
        """
