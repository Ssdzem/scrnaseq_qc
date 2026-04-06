rule integration_annotation_all:
    input:
        ia_config=IA_CONFIG_PATH,
        evaluated_done="out/integration_annotation/all_studies/.evaluated_done",
    output:
        done=touch("out/integration_annotation/all_studies/.integration_annotation_all_done")


rule integration_annotation_study:
    input:
        ia_config=IA_CONFIG_PATH,
        evaluated_done="out/integration_annotation/{study}/.evaluated_done",
    output:
        done=touch("out/integration_annotation/{study}/.integration_annotation_done")
