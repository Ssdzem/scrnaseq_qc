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


rule standardize_study:
    input:
        spec=lambda wc: str(_spec_path_for_study(wc.study)),
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
        $CONDA_PREFIX/bin/python workflow/modules/standardize/scripts/standardize.py \
          --spec {input.spec} \
          --study {params.study} \
          --output {output.h5ad} \
          --log {log}
        """
