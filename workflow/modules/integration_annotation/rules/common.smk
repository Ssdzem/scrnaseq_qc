from pathlib import Path

import yaml


IA_CONFIG_PATH = str(Path("config/integration_annotation.yaml"))


def _load_ia_config(path: str) -> dict:
    try:
        with open(path, "r", encoding="utf-8") as handle:
            return yaml.safe_load(handle) or {}
    except FileNotFoundError:
        return {}


IA_FULL_CFG = config if isinstance(config, dict) and config else _load_ia_config(IA_CONFIG_PATH)
IA_INTEGRATION_CFG = IA_FULL_CFG.get("integration", {}) if isinstance(IA_FULL_CFG, dict) else {}
IA_ANNOTATION_CFG = IA_FULL_CFG.get("annotation", {}) if isinstance(IA_FULL_CFG, dict) else {}
IA_REFERENCE_CFG = IA_FULL_CFG.get("reference_atlas", {}) if isinstance(IA_FULL_CFG, dict) else {}


def _ia_reference_h5ad() -> str:
    # The integration_annotation module intentionally uses the configured Liver Cell Atlas build target.
    # Annotation scripts receive this resolved path from Snakemake instead of accepting arbitrary atlas paths.
    return str(IA_REFERENCE_CFG.get("output_h5ad", "data/ref/lca_human_reference.h5ad")).strip()


def _ia_reference_report() -> str:
    return str(IA_REFERENCE_CFG.get("report", "data/ref/reference_build_report.json")).strip()


def _ia_reference_ready_marker() -> str:
    return str(IA_REFERENCE_CFG.get("ready_marker", "data/ref/.reference_ready")).strip()


IA_REFERENCE_H5AD = _ia_reference_h5ad()
IA_REFERENCE_REPORT = _ia_reference_report()
IA_REFERENCE_READY_MARKER = _ia_reference_ready_marker()


def _ia_reference_targets(_wildcards=None):
    enabled = bool(IA_REFERENCE_CFG.get("enabled", False))
    method = str(IA_ANNOTATION_CFG.get("method", "")).strip().lower()
    if enabled and method in {"symphony", "scanvi"}:
        return [IA_REFERENCE_H5AD, IA_REFERENCE_READY_MARKER]
    return []


def _studies_to_process():
    studies = config.get("inputs", {}).get("studies", [])
    if studies:
        return list(studies)

    meta = config.get("inputs", {}).get("metadata_csv") or "config/liver_scRNAseq_metadata_everything.csv"
    try:
        import csv
        import re

        def _is_missing(v: str) -> bool:
            return (v or "").strip() == "" or (v or "").strip().lower() == "na"

        with open(meta, newline="", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle)
            studies_set = set()
            for row in reader:
                study = (row.get("study_accession") or "").strip()
                if not _is_missing(study):
                    studies_set.add(re.sub(r"\s+", "", study))
        return sorted(studies_set)
    except FileNotFoundError:
        return []


def _ia_integration_method(_wildcards=None) -> str:
    integration_cfg = IA_INTEGRATION_CFG
    if not bool(integration_cfg.get("enabled", True)):
        return "none"
    method = str(integration_cfg.get("method", integration_cfg.get("integration_method", ""))).strip().lower()
    valid = {"sysvi", "scvi", "harmony"}
    if method not in valid:
        raise ValueError(
            "config/integration_annotation.yaml requires integration.method to be set before running "
            f"integration_annotation. Choose one of {sorted(valid)}; found {method!r}."
        )
    return method


def _ia_rapids_backend() -> str:
    raw = str(IA_INTEGRATION_CFG.get("backend", IA_INTEGRATION_CFG.get("preprocessing_backend", "cpu"))).strip().lower()
    if raw not in {"cpu", "gpu", "auto"}:
        raise ValueError(
            "config/integration_annotation.yaml requires integration.backend to be one of "
            "['cpu', 'gpu', 'auto']; found "
            f"{raw!r}."
        )
    return raw


def _ia_integrate_env() -> str:
    accelerator = str(IA_INTEGRATION_CFG.get("accelerator", "")).strip().lower()
    use_gpu = IA_INTEGRATION_CFG.get("use_gpu", None)
    if isinstance(use_gpu, str):
        use_gpu = use_gpu.strip().lower() in {"1", "true", "yes", "y", "gpu"}
    wants_scvi_gpu = accelerator == "gpu" or bool(use_gpu)
    wants_rapids_gpu = _ia_rapids_backend() in {"gpu", "auto"}
    return "../env/integration_annotation_gpu_py.yaml" if (wants_scvi_gpu or wants_rapids_gpu) else "../env/integration_annotation_py.yaml"


def _ia_log(study: str) -> str:
    return f"out/logs/integration_annotation/{study}.log"


STUDIES = _studies_to_process()
