import os
import shutil
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
# Import helper functions from `scripts/fetch_data.py`.
#
# `scripts` is not necessarily a Python package (no __init__.py), so importing
# it as `scripts.fetch_data` can fail depending on how Snakemake sets up the
# module search path. Add the scripts directory to sys.path and import the
# module directly instead.
SCRIPTS_DIR = os.path.join(MODULE_DIR, "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from fetch_data import list_studies, read_rows  # noqa: E402


METADATA_CSV = _rooted("metadata_csv")
OUT_ROOT = _rooted("out_root")

ASPERA_CFG = config.get("aspera", {})


def _aspera_enabled() -> bool:
    val = ASPERA_CFG.get("enabled")
    if val is None:
        return False
    if isinstance(val, bool):
        return val
    sval = str(val).strip()
    if sval == "":
        return False
    return sval.lower() in {"1", "true", "yes", "y", "on"}


def _aspera_val(key: str, env_key: str, default: str = "") -> str:
    v = (ASPERA_CFG.get(key) or "").strip()
    if v:
        return v
    return os.environ.get(env_key, default)


def _aspera_ascp_path() -> str:
    """Resolve the ascp binary path.

    Precedence:
      1) env var ASCP_PATH
      2) config.aspera.ascp_path
      3) `shutil.which("ascp")`
    """

    return (
        os.environ.get("ASCP_PATH")
        or (ASPERA_CFG.get("ascp_path") or "").strip()
        or (shutil.which("ascp") or "")
    )


def _aspera_key_path() -> str:
    """Resolve the Aspera SSH key path.

    Precedence:
      1) env var ASCP_KEY
      2) config.aspera.key_path
    """

    return os.environ.get("ASCP_KEY") or (ASPERA_CFG.get("key_path") or "").strip()


ROWS = read_rows(Path(METADATA_CSV))
assert ROWS, f"No rows parsed from metadata: {METADATA_CSV}"

STUDIES = sorted({r.get("study_accession", "").strip() for r in ROWS if r.get("study_accession")})
assert STUDIES, f"No studies parsed from metadata: {METADATA_CSV}"

RUNS_BY_STUDY = {}
TRIPLES_BY_STUDY = {}
for r in ROWS:
    study = (r.get("study_accession") or "").strip()
    run_id = (r.get("run_id") or "").strip()
    donor_id = (r.get("donor_id") or "").strip()
    sample_id = (r.get("sample_id") or "").strip()
    if not study or not run_id:
        continue
    RUNS_BY_STUDY.setdefault(study, set()).add(run_id)
    if donor_id and sample_id:
        TRIPLES_BY_STUDY.setdefault(study, set()).add((donor_id, sample_id, run_id))



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
            + ("--ascp-path '" + _aspera_ascp_path() + "' " if _aspera_enabled() else "")
            + ("--ascp-key '" + _aspera_key_path() + "' " if _aspera_enabled() else "")
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
            + ("--ascp-path '" + _aspera_ascp_path() + "' " if _aspera_enabled() else "")
            + ("--ascp-key '" + _aspera_key_path() + "' " if _aspera_enabled() else "")
            + ("--aspera-host '" + _aspera_val("host", "ASPERA_HOST") + "' " if _aspera_enabled() else "")
            + ("--aspera-remote-root '" + _aspera_val("remote_root", "ASPERA_REMOTE_ROOT") + "' " if _aspera_enabled() else "")
            + ("--aspera-port " + str(ASPERA_CFG.get("port", 33001)) + " " if _aspera_enabled() else "")
            + ("--aspera-rate '" + str(ASPERA_CFG.get("rate_limit", "100m")) + "'" if _aspera_enabled() else "")
        )


# Allow targeting `data/<study>/.fetch_done` directly when OUT_ROOT is something
# else (e.g. absolute path). This keeps CLI usage intuitive.
rule fetch_done_data:
    input:
        lambda wc: [
            f"out/{wc.study}/fastqs/{donor_id}/{sample_id}/{run_id}/.fetch_done"
            for donor_id, sample_id, run_id in sorted(TRIPLES_BY_STUDY.get(wc.study, set()))
        ]
    output:
        touch("out/{study,[^/]+}/.fetch_done"),
    conda:
        os.path.join(MODULE_DIR, "env", "fetch.yaml")
    shell:
        "true"



# Convenience: allow targeting a single run via
#   out/<study>/fastqs/<donor_id>/<sample_id>/<run_id>/.fetch_done
rule fetch_run_data:
    input:
        metadata=METADATA_CSV,
        script=os.path.join(MODULE_DIR, "scripts", "fetch_data.py"),
    output:
        touch("out/{study,[^/]+}/fastqs/{donor_id}/{sample_id}/{run_id}/.fetch_done"),
    params:
        out_root="out",
    conda:
        os.path.join(MODULE_DIR, "env", "fetch.yaml")
    shell:
        (
            "python {input.script} --metadata {input.metadata} --out-root {params.out_root} download-run {wildcards.study} {wildcards.run_id} "
            + ("--aspera " if _aspera_enabled() else "")
            + ("--ascp-path '" + _aspera_ascp_path() + "' " if _aspera_enabled() else "")
            + ("--ascp-key '" + _aspera_key_path() + "' " if _aspera_enabled() else "")
            + ("--aspera-host '" + _aspera_val("host", "ASPERA_HOST") + "' " if _aspera_enabled() else "")
            + ("--aspera-remote-root '" + _aspera_val("remote_root", "ASPERA_REMOTE_ROOT") + "' " if _aspera_enabled() else "")
            + ("--aspera-port " + str(ASPERA_CFG.get("port", 33001)) + " " if _aspera_enabled() else "")
            + ("--aspera-rate '" + str(ASPERA_CFG.get("rate_limit", "100m")) + "'" if _aspera_enabled() else "")
        )
