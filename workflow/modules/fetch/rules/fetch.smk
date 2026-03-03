import os
import sys
from pathlib import Path


# Prefer the project-level config when included from the main workflow.
# When running this module standalone (via `workflow/modules/fetch/Snakefile`),
# fall back to the module-local config path.
_project_cfg = os.path.join(workflow.basedir, "config", "fetch.yaml")
_module_cfg = os.path.join(os.path.dirname(__file__), "..", "config", "fetch.yaml")
_cfg = _project_cfg if os.path.exists(_project_cfg) else _module_cfg
configfile: _cfg


def _norm(path: str) -> str:
    return os.path.normpath(os.path.expanduser(path))


_cfg_root = os.path.dirname(_cfg)
_project_root_cfg = config.get("paths", {}).get("project_root", ".")
if os.path.isabs(os.path.expanduser(_project_root_cfg)):
    PROJECT_ROOT = _norm(_project_root_cfg)
else:
    # Interpret project_root relative to the configfile location.
    PROJECT_ROOT = _norm(os.path.join(_cfg_root, _project_root_cfg))


def _rooted(key: str) -> str:
    value = config["paths"][key]
    expanded = os.path.expanduser(value)
    if os.path.isabs(expanded):
        return _norm(expanded)
    return _norm(os.path.join(PROJECT_ROOT, value))


MODULE_DIR = _rooted("module_dir")
sys.path.insert(0, MODULE_DIR)

from scripts.fetch_data import list_studies  # noqa: E402


METADATA_CSV = _rooted("metadata_csv")
OUT_ROOT = _rooted("out_root")


STUDIES = list_studies(Path(METADATA_CSV))
assert STUDIES, f"No studies parsed from metadata: {METADATA_CSV}"


rule fetch_all:
    input:
        expand(f"{OUT_ROOT}/{{study}}/.fetch_done", study=STUDIES)


rule fetch_done:
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


rule fetch_study:
    input:
        f"{OUT_ROOT}/{{study}}/.fetch_done"
