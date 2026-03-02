import os
import sys
from pathlib import Path


# `workflow.basedir` is the directory of the top-level Snakefile.
# In this project layout, that is the project root (`a1_p1_qc/`).
configfile: os.path.join(workflow.basedir, "config", "fetch.yaml")


def _norm(path: str) -> str:
    return os.path.normpath(os.path.expanduser(path))


PROJECT_ROOT = _norm(config.get("paths", {}).get("project_root", workflow.basedir))


def _rooted(key: str) -> str:
    value = config["paths"][key]
    if os.path.isabs(os.path.expanduser(value)):
        return _norm(value)
    return _norm(os.path.join(PROJECT_ROOT, value))


MODULE_DIR = _rooted("module_dir")
sys.path.insert(0, MODULE_DIR)

from scripts.fetch_data import list_studies  # noqa: E402


METADATA_CSV = _rooted("metadata_csv")
OUT_ROOT = _rooted("out_root")


if Path(METADATA_CSV).exists():
    STUDIES = list_studies(Path(METADATA_CSV))
else:
    # Allow including this module even if fetch inputs
    # are not configured/present yet.
    STUDIES = []


rule fetch_all:
    input:
        expand(f"{OUT_ROOT}/{{study}}/.fetch_done", study=STUDIES)


rule fetch_study:
    input:
        metadata=METADATA_CSV,
        script=os.path.join(MODULE_DIR, "scripts", "fetch_data.py"),
    output:
        touch(f"{OUT_ROOT}/{{study}}/.fetch_done"),
    params:
        out_root=OUT_ROOT,
    conda:
        os.path.join(MODULE_DIR, "env", "fetch.yaml")
    shell:
        "python {input.script} --metadata {input.metadata} --out-root {params.out_root} download-study {wildcards.study}"
