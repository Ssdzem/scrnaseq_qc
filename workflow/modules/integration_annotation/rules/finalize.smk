rule integration_annotation_all:
    input:
        ia_config=IA_CONFIG_PATH,
        evaluated_done=f"{IA_ALL_OUTPUT_DIR}/.evaluated_done",
    output:
        done=touch(f"{IA_ALL_OUTPUT_DIR}/.integration_annotation_all_done")


rule integration_annotation_study:
    wildcard_constraints:
        study=r"(?!all_studies$)[^/]+"
    input:
        ia_config=IA_CONFIG_PATH,
        evaluated_done=f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/.{{study}}_evaluated_done",
    output:
        done=touch(f"{IA_STUDY_OUTPUT_DIR_TEMPLATE}/.{{study}}_integration_annotation_done")
