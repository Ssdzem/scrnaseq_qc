rule integration_annotation_all:
    input:
        ia_config=IA_CONFIG_PATH,
        evaluated_done="out/integration_annotation/all_studies/.evaluated_done",
    output:
        done=touch("out/integration_annotation/all_studies/.integration_annotation_all_done")


rule integration_annotation_study:
    wildcard_constraints:
        study=r"(?!all_studies$)[^/]+"
    input:
        ia_config=IA_CONFIG_PATH,
        evaluated_done="out/integration_annotation/{study}/.{study}_evaluated_done",
    output:
        done=touch("out/integration_annotation/{study}/.{study}_integration_annotation_done")
