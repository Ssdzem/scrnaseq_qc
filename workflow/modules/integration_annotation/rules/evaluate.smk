rule ia_evaluate_all:
    input:
        annotated="out/integration_annotation/all_studies/all_annotated.h5ad",
        ia_config=IA_CONFIG_PATH,
        annotated_done="out/integration_annotation/all_studies/.annotated_done",
    output:
        metrics_json="out/integration_annotation/all_studies/evaluation_metrics.json",
        metrics_csv="out/integration_annotation/all_studies/evaluation_metrics.csv",
        done=touch("out/integration_annotation/all_studies/.evaluated_done"),
    log:
        _ia_log("all_studies")
    conda:
        "../env/integration_annotation_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        PYTHONNOUSERSITE=1 $CONDA_PREFIX/bin/python -s workflow/modules/integration_annotation/scripts/evaluate_integration.py \
          --scope all_studies \
          --input {input.annotated} \
          --config {input.ia_config} \
          --output-json {output.metrics_json} \
          --output-csv {output.metrics_csv} \
          --log {log}
        """


rule ia_evaluate_study:
    input:
        annotated="out/integration_annotation/{study}/{study}_annotated.h5ad",
        ia_config=IA_CONFIG_PATH,
        annotated_done="out/integration_annotation/{study}/.annotated_done",
    output:
        metrics_json="out/integration_annotation/{study}/evaluation_metrics.json",
        metrics_csv="out/integration_annotation/{study}/evaluation_metrics.csv",
        done=touch("out/integration_annotation/{study}/.evaluated_done"),
    params:
        study=lambda wc: wc.study,
    log:
        _ia_log("{study}")
    conda:
        "../env/integration_annotation_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        PYTHONNOUSERSITE=1 $CONDA_PREFIX/bin/python -s workflow/modules/integration_annotation/scripts/evaluate_integration.py \
          --scope study \
          --study {params.study} \
          --input {input.annotated} \
          --config {input.ia_config} \
          --output-json {output.metrics_json} \
          --output-csv {output.metrics_csv} \
          --log {log}
        """
