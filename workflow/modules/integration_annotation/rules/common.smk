from pathlib import Path

import yaml


IA_CONFIG_PATH = str(Path("config/integration_annotation.yaml"))


def _load_ia_config(path: str) -> dict:
    try:
        with open(path, "r", encoding="utf-8") as handle:
            return yaml.safe_load(handle) or {}
    except FileNotFoundError:
        return {}


IA_FULL_CFG = _load_ia_config(IA_CONFIG_PATH)
IA_ANNOTATION_CFG = IA_FULL_CFG.get("annotation", {}) if isinstance(IA_FULL_CFG, dict) else {}
IA_REFERENCE_CFG = IA_FULL_CFG.get("reference_atlas", {}) if isinstance(IA_FULL_CFG, dict) else {}


def _ia_reference_h5ad() -> str:
    method = str(IA_ANNOTATION_CFG.get("method", "")).strip().lower()
    if method == "symphony":
        return str(IA_ANNOTATION_CFG.get("symphony", {}).get("reference_h5ad", "")).strip()
    if method == "scanvi":
        return str(IA_ANNOTATION_CFG.get("scanvi", {}).get("reference_h5ad", "")).strip()
    return ""


IA_REFERENCE_H5AD = _ia_reference_h5ad() or "out/integration_annotation/reference/lca_human_reference.h5ad"


def _ia_reference_targets():
    enabled = bool(IA_REFERENCE_CFG.get("enabled", False))
    method = str(IA_ANNOTATION_CFG.get("method", "")).strip().lower()
    if enabled and method in {"symphony", "scanvi"}:
        return [IA_REFERENCE_H5AD, "out/integration_annotation/reference/.reference_ready"]
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


def _ia_log(study: str) -> str:
    return f"out/logs/integration_annotation/{study}.log"


STUDIES = _studies_to_process()
