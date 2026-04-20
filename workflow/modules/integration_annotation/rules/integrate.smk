rule ia_integrate_all:
    input:
        pre_integration="out/all_studies/all_pre_integration.h5ad",
        ia_config=IA_CONFIG_PATH,
    output:
        integrated="out/integration_annotation/all_studies/all_integrated.h5ad",
        report="out/integration_annotation/all_studies/integration_report.json",
        done=touch("out/integration_annotation/all_studies/.integrated_done"),
    params:
        integration_method=_ia_integration_method,
    log:
        _ia_log("all_studies")
    conda:
        "../env/integration_annotation_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        PYTHONNOUSERSITE=1 $CONDA_PREFIX/bin/python -s workflow/modules/integration_annotation/scripts/integrate_embeddings.py \
          --scope all_studies \
          --input {input.pre_integration} \
          --config {input.ia_config} \
          --output {output.integrated} \
          --report {output.report} \
          --log {log}
        """


rule ia_integrate_study:
    input:
        pre_integration="out/{study}/qc/{study}_pre_integration.h5ad",
        ia_config=IA_CONFIG_PATH,
    output:
        integrated="out/integration_annotation/{study}/{study}_integrated.h5ad",
        report="out/integration_annotation/{study}/integration_report.json",
        done=touch("out/integration_annotation/{study}/.integrated_done"),
    params:
        study=lambda wc: wc.study,
        integration_method=_ia_integration_method,
    log:
        _ia_log("{study}")
    conda:
        "../env/integration_annotation_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        PYTHONNOUSERSITE=1 $CONDA_PREFIX/bin/python -s workflow/modules/integration_annotation/scripts/integrate_embeddings.py \
          --scope study \
          --study {params.study} \
          --input {input.pre_integration} \
          --config {input.ia_config} \
          --output {output.integrated} \
          --report {output.report} \
          --log {log}
        """
