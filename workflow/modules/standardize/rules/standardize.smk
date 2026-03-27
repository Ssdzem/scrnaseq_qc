import os
from pathlib import Path

import yaml

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


STUDIES = _studies_to_process()


rule standardize_all:
    input:
        expand("{h5ad}", h5ad=[_study_out_h5ad(s) for s in STUDIES])


def _required_spec_input_paths(study: str) -> list[str]:
    spec_path = _spec_path_for_study(study)
    spec = yaml.safe_load(open(spec_path, "r", encoding="utf-8")) or {}
    inp = spec.get("input") or {}

    paths: list[str] = []

    # Flat key/value paths (e.g., seurat_rds / count_matrix_10x specs).
    for v in inp.values():
        if isinstance(v, str):
            paths.append(v)

    # Nested sample list paths (e.g., tenx_h5_multi specs).
    for sample in inp.get("samples", []) if isinstance(inp, dict) else []:
        if isinstance(sample, dict):
            h5 = sample.get("h5")
            if isinstance(h5, str):
                paths.append(h5)

    # Keep deterministic ordering and remove duplicates.
    return sorted(set(paths))


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
        lambda wc: "../env/standardize_r.yaml"
        if (yaml.safe_load(open(_spec_path_for_study(wc.study))).get("process") == "seurat_rds")
        else "../env/standardize_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        $CONDA_PREFIX/bin/python -u workflow/modules/standardize/scripts/standardize.py \
          --spec {input.spec} \
          --study {params.study} \
          --output {output.h5ad} \
          --log {log}
        """
