import os
import sys
from pathlib import Path

import yaml

# Resolve the module-local scripts directory robustly regardless of whether
# this file is included from the project root Snakefile or run standalone.
_THIS_FILE = Path(__file__).resolve()
_MODULE_DIR = _THIS_FILE.parents[1]  # .../workflow/modules/standardize
_SCRIPTS_DIR = _MODULE_DIR / "scripts"
sys.path.insert(0, str(_SCRIPTS_DIR))
from cleaner_registry import CLEANER_REGISTRY, normalize_cleaner_name


def _spec_path_for_study(study: str) -> Path:
    return Path("config") / "standardize" / f"{study}.yaml"


def _studies_to_process() -> list[str]:
    root = Path("config") / "standardize"
    if not root.exists():
        return []
    return sorted([p.stem for p in root.glob("*.yaml")])


def _study_out_h5ad(study: str) -> str:
    # Canonical output location is fixed by convention.
    return f"out/{study}/raw.h5ad"


def _load_spec(study: str) -> dict:
    spec_path = _spec_path_for_study(study)
    return yaml.safe_load(open(spec_path, "r", encoding="utf-8")) or {}


def _pipeline_cleaners(study: str) -> list[str]:
    spec = _load_spec(study)
    pipeline = spec.get("pipeline")
    if pipeline is None:
        pipeline = spec.get("Pipeline")
    if pipeline is None:
        raise WorkflowError(f"Missing required 'pipeline' in spec for study {study}")

    if isinstance(pipeline, list):
        cleaners = [normalize_cleaner_name(v) for v in pipeline if str(v).strip()]
    elif isinstance(pipeline, dict):
        cleaners = [normalize_cleaner_name(v) for _k, v in pipeline.items() if str(v).strip()]
    else:
        raise WorkflowError(f"Invalid pipeline format in spec for study {study}")

    unknown = [c for c in cleaners if c not in CLEANER_REGISTRY]
    if unknown:
        raise WorkflowError(f"Unknown cleaners in pipeline for {study}: {unknown}")
    return cleaners


def _required_spec_input_paths(study: str) -> list[str]:
    spec = _load_spec(study)
    inp = spec.get("input") or {}

    paths: set[str] = set()
    for cleaner in _pipeline_cleaners(study):
        reg = CLEANER_REGISTRY[cleaner]
        for key in reg.path_keys:
            if key.endswith("[]"):
                continue
            if key == "input.samples[].h5":
                for sample in inp.get("samples", []) if isinstance(inp, dict) else []:
                    if isinstance(sample, dict):
                        h5 = sample.get("h5")
                        if isinstance(h5, str) and h5.strip():
                            paths.add(h5)
                continue
            if key == "input.samples[].h5ad":
                for sample in inp.get("samples", []) if isinstance(inp, dict) else []:
                    if isinstance(sample, dict):
                        h5ad = sample.get("h5ad")
                        if isinstance(h5ad, str) and h5ad.strip() and not h5ad.startswith("in-memory://"):
                            paths.add(h5ad)
                continue

            cur = spec
            for part in key.split("."):
                if part == "input":
                    cur = inp
                elif isinstance(cur, dict):
                    cur = cur.get(part)
                else:
                    cur = None
                    break
            if isinstance(cur, str) and cur.strip() and not cur.startswith("in-memory://"):
                paths.add(cur)

    return sorted(paths)


def _pipeline_env(study: str) -> str:
    cleaners = _pipeline_cleaners(study)
    needs_r = any(CLEANER_REGISTRY[c].env == "r" for c in cleaners)
    return "../env/standardize_r.yaml" if needs_r else "../env/standardize_py.yaml"


STUDIES = _studies_to_process()


rule standardize_all:
    input:
        expand("{h5ad}", h5ad=[_study_out_h5ad(s) for s in STUDIES])


rule standardize_unpack_study:
    input:
        spec=lambda wc: str(_spec_path_for_study(wc.study)),
        fetch_done="out/{study}/.fetch_done",
    output:
        done=touch("out/{study}/.unpack_done"),
    params:
        study=lambda wc: wc.study,
    log:
        "out/logs/standardize/{study}.unpack.log"
    conda:
        "../env/standardize_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        $CONDA_PREFIX/bin/python -u workflow/modules/standardize/scripts/unpack_inputs.py \
          --study {params.study} \
          --outdir out \
          --log {log}
        """


rule standardize_validate_inputs:
    input:
        spec=lambda wc: str(_spec_path_for_study(wc.study)),
        unpack_done="out/{study}/.unpack_done",
    output:
        done=touch("out/{study}/.standardize_inputs_validated_done"),
    run:
        required = _required_spec_input_paths(wildcards.study)
        missing = [p for p in required if not Path(p).exists()]
        if missing:
            raise WorkflowError(
                f"Missing standardize input files for study {wildcards.study}: "
                + ", ".join(missing)
            )


rule standardize_study:
    input:
        spec=lambda wc: str(_spec_path_for_study(wc.study)),
        validated_done="out/{study}/.standardize_inputs_validated_done",
    output:
        h5ad=_study_out_h5ad("{study}"),
    params:
        study=lambda wc: wc.study,
    log:
        "out/logs/standardize/{study}.log"
    conda:
        lambda wc: _pipeline_env(wc.study)
    shell:
        r"""
        mkdir -p $(dirname {log})
        $CONDA_PREFIX/bin/python -u workflow/modules/standardize/scripts/standardize.py \
          --spec {input.spec} \
          --study {params.study} \
          --output {output.h5ad} \
          --log {log}
        """
