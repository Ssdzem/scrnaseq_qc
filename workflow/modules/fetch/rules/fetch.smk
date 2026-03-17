import os
import sys
from pathlib import Path


# Prefer the project-level config when included from the main workflow.
# When running this module standalone (via `workflow/modules/fetch/Snakefile`),
# fall back to the module-local config path.
_project_cfg = os.path.join(workflow.basedir, "config", "fetch.yaml")

# `workflow.basedir` can be surprising when running a module Snakefile directly,
# so also try resolving the project config relative to CWD.
_cwd_project_cfg = os.path.join(os.getcwd(), "config", "fetch.yaml")

_candidates = [_project_cfg, _cwd_project_cfg]
_cfg = next((p for p in _candidates if os.path.exists(p)), None)
if _cfg is None:
    raise WorkflowError(
        "Could not find fetch config. Expected one of: "
        + ", ".join(_candidates)
        + ". If you run this module standalone, pass --configfile config/fetch.yaml."
    )

configfile: _cfg

# The module is meant to be run from the `a1_p1_qc/` directory.
# If invoked via `-s workflow/modules/fetch/Snakefile` from elsewhere,
# Snakemake's path resolution can become ambiguous.
if os.path.basename(os.getcwd()) != "a1_p1_qc":
    raise WorkflowError(
        "Run the fetch module from the a1_p1_qc directory, e.g. `cd a1_p1_qc` "
        "and then `snakemake -s workflow/modules/fetch/Snakefile ...`."
    )


def _norm(path: str) -> str:
    return os.path.normpath(os.path.expanduser(path))


_cfg_root = os.path.dirname(_cfg)
_project_root_cfg = config.get("paths", {}).get("project_root", ".")

# `project_root` is treated as repo-root relative for this project.
# In practice this means running from `a1_p1_qc/` should keep relative paths
# like `config/...` resolving as expected.
PROJECT_ROOT = _norm(os.getcwd())


def _rooted(key: str) -> str:
    value = config["paths"][key]
    expanded = os.path.expanduser(value)
    if os.path.isabs(expanded):
        return _norm(expanded)
    # Interpret relative paths as relative to project root.
    return _norm(os.path.join(PROJECT_ROOT, expanded))


MODULE_DIR = _rooted("module_dir")
sys.path.insert(0, MODULE_DIR)

from scripts.fetch_data import list_studies  # noqa: E402


METADATA_CSV = _rooted("metadata_csv")
OUT_ROOT = _rooted("out_root")

ASPERA_CFG = config.get("aspera", {})


def _aspera_enabled() -> bool:
    val = (ASPERA_CFG.get("enabled") or "").strip()
    if val == "":
        return False
    return val.lower() in {"1", "true", "yes", "y", "on"}


def _aspera_val(key: str, env_key: str, default: str = "") -> str:
    v = (ASPERA_CFG.get(key) or "").strip()
    if v:
        return v
    return os.environ.get(env_key, default)


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
        (
            "python {input.script} --metadata {input.metadata} --out-root {params.out_root} download-study {wildcards.study} "
            + ("--aspera " if _aspera_enabled() else "")
            + ("--ascp-path '" + _aspera_val("ascp_path", "ASCP_PATH") + "' " if _aspera_enabled() else "")
            + ("--ascp-key '" + _aspera_val("key_path", "ASCP_KEY") + "' " if _aspera_enabled() else "")
            + ("--aspera-host '" + _aspera_val("host", "ASPERA_HOST") + "' " if _aspera_enabled() else "")
            + ("--aspera-remote-root '" + _aspera_val("remote_root", "ASPERA_REMOTE_ROOT") + "' " if _aspera_enabled() else "")
            + ("--aspera-port " + str(ASPERA_CFG.get("port", 33001)) + " " if _aspera_enabled() else "")
            + ("--aspera-rate '" + str(ASPERA_CFG.get("rate_limit", "100m")) + "'" if _aspera_enabled() else "")
        )


rule fetch_study:
    input:
        f"{OUT_ROOT}/{{study}}/.fetch_done"


rule fetch_run:
    input:
        metadata=METADATA_CSV,
        script=os.path.join(MODULE_DIR, "scripts", "fetch_data.py"),
    output:
        touch(f"{OUT_ROOT}/{{study}}/fastqs/{{run_id}}/.fetch_done"),
    params:
        out_root=OUT_ROOT,
    conda:
        os.path.join(MODULE_DIR, "env", "fetch.yaml")
    shell:
        (
            "python {input.script} --metadata {input.metadata} --out-root {params.out_root} download-run {wildcards.study} {wildcards.run_id} "
            + ("--aspera " if _aspera_enabled() else "")
            + ("--ascp-path '" + _aspera_val("ascp_path", "ASCP_PATH") + "' " if _aspera_enabled() else "")
            + ("--ascp-key '" + _aspera_val("key_path", "ASCP_KEY") + "' " if _aspera_enabled() else "")
            + ("--aspera-host '" + _aspera_val("host", "ASPERA_HOST") + "' " if _aspera_enabled() else "")
            + ("--aspera-remote-root '" + _aspera_val("remote_root", "ASPERA_REMOTE_ROOT") + "' " if _aspera_enabled() else "")
            + ("--aspera-port " + str(ASPERA_CFG.get("port", 33001)) + " " if _aspera_enabled() else "")
            + ("--aspera-rate '" + str(ASPERA_CFG.get("rate_limit", "100m")) + "'" if _aspera_enabled() else "")
        )


# Allow targeting `data/<study>/.fetch_done` directly when OUT_ROOT is something
# else (e.g. absolute path). This keeps CLI usage intuitive.
rule fetch_done_data:
    input:
        metadata=METADATA_CSV,
        script=os.path.join(MODULE_DIR, "scripts", "fetch_data.py"),
    output:
        touch("out/{study}/.fetch_done"),
    params:
        out_root="out",
    conda:
        os.path.join(MODULE_DIR, "env", "fetch.yaml")
    shell:
        (
            "python {input.script} --metadata {input.metadata} --out-root {params.out_root} download-study {wildcards.study} "
            + ("--aspera " if _aspera_enabled() else "")
            + ("--ascp-path '" + _aspera_val("ascp_path", "ASCP_PATH") + "' " if _aspera_enabled() else "")
            + ("--ascp-key '" + _aspera_val("key_path", "ASCP_KEY") + "' " if _aspera_enabled() else "")
            + ("--aspera-host '" + _aspera_val("host", "ASPERA_HOST") + "' " if _aspera_enabled() else "")
            + ("--aspera-remote-root '" + _aspera_val("remote_root", "ASPERA_REMOTE_ROOT") + "' " if _aspera_enabled() else "")
            + ("--aspera-port " + str(ASPERA_CFG.get("port", 33001)) + " " if _aspera_enabled() else "")
            + ("--aspera-rate '" + str(ASPERA_CFG.get("rate_limit", "100m")) + "'" if _aspera_enabled() else "")
        )
