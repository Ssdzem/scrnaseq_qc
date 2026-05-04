rule ia_prepare_reference_atlas:
    input:
        ia_config=IA_CONFIG_PATH,
    output:
        ref_h5ad=IA_REFERENCE_H5AD,
        report=IA_REFERENCE_REPORT,
        done=touch(IA_REFERENCE_READY_MARKER),
    log:
        _ia_log("reference_atlas")
    conda:
        "../env/integration_annotation_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        PYTHONNOUSERSITE=1 $CONDA_PREFIX/bin/python -s workflow/modules/integration_annotation/scripts/fetch_reference_atlas.py \
          --config {input.ia_config} \
          --output {output.ref_h5ad} \
          --report {output.report} \
          --log {log}
        """
