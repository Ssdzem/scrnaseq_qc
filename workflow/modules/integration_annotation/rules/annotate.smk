rule ia_annotate_all:
    input:
        integrated="out/integration_annotation/all_studies/all_integrated.h5ad",
        ia_config=IA_CONFIG_PATH,
        integrated_done="out/integration_annotation/all_studies/.integrated_done",
        reference=_ia_reference_targets,
    output:
        annotated="out/integration_annotation/all_studies/all_annotated.h5ad",
        summary="out/integration_annotation/all_studies/annotation_summary.csv",
        done=touch("out/integration_annotation/all_studies/.annotated_done"),
    params:
        reference_h5ad=IA_REFERENCE_H5AD,
    log:
        _ia_log("all_studies")
    conda:
        "../env/integration_annotation_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
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
    input:
        integrated="out/integration_annotation/{study}/{study}_integrated.h5ad",
        ia_config=IA_CONFIG_PATH,
        integrated_done="out/integration_annotation/{study}/.integrated_done",
        reference=_ia_reference_targets,
    output:
        annotated="out/integration_annotation/{study}/{study}_annotated.h5ad",
        summary="out/integration_annotation/{study}/annotation_summary.csv",
        done=touch("out/integration_annotation/{study}/.annotated_done"),
    params:
        study=lambda wc: wc.study,
        reference_h5ad=IA_REFERENCE_H5AD,
    log:
        _ia_log("{study}")
    conda:
        "../env/integration_annotation_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
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
