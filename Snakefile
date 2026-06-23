# Project entrypoint Snakefile.

from pathlib import Path


def _profile_module_config(filename: str, explicit_key: str) -> str:
    explicit = str(config.get(explicit_key, "")).strip()
    if explicit:
        return explicit

    profile_dir = str(config.get("profile_config_dir", "")).strip()
    if profile_dir:
        profile_cfg = Path(profile_dir) / filename
        if profile_cfg.exists():
            return str(profile_cfg)

    return str(Path(workflow.basedir) / "config" / filename)


QC_CONFIG_PATH = _profile_module_config("qc.yaml", "qc_config")

configfile: QC_CONFIG_PATH

include: "workflow/modules/fetch/Snakefile"
include: "workflow/modules/fastq_prepros/Snakefile"
include: "workflow/modules/standardize/Snakefile"
include: "workflow/modules/qc_preproc/Snakefile"
include: "workflow/modules/integration_annotation/Snakefile"


def _study_from_config() -> str:
    study = (config.get("study") or "").strip()
    if not study:
        raise WorkflowError(
            "Missing study for qc_preintegration_study. "
            "Run with: snakemake qc_preintegration_study --config study=<STUDY_ACCESSION>"
        )
    return study


# Choose what runs by default.
# - `snakemake`            -> runs both fetch + qc
# - `snakemake fetch_all`  -> runs only fetch module
# - `snakemake qc_all`     -> runs only qc module
# - `snakemake qc_preintegration_study --config study=<STUDY>`
#                         -> runs fetch -> standardize -> study-level qc pre-integration for one study
rule all:
    input:
        rules.fetch_all.input,
        rules.fastq_prepros_all.input,
        rules.qc_all.input,


rule qc_preintegration_study:
    input:
        lambda wc: f"out/{_study_from_config()}/qc/.study_pre_integration_done"


rule integration_annotation_study_from_config:
    input:
        lambda wc: f"out/integration_annotation/{_study_from_config()}/.{_study_from_config()}_integration_annotation_done"
