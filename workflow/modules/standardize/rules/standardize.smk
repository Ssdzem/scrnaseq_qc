import os
import sys
from pathlib import Path

import yaml

# Resolve scripts dir with multiple candidates because Snakemake can evaluate
# included files with different base-path semantics depending on invocation.
_candidate_scripts_dirs = [
    # Included from repo root Snakefile.
    Path(workflow.basedir) / "workflow" / "modules" / "standardize" / "scripts",
    # Running module Snakefile directly (workflow.basedir == module dir).
    Path(workflow.basedir) / "scripts",
    # Fallback from current working directory.
    Path(os.getcwd()) / "workflow" / "modules" / "standardize" / "scripts",
]
_SCRIPTS_DIR = next((p.resolve() for p in _candidate_scripts_dirs if p.exists()), None)
if _SCRIPTS_DIR is None:
    raise WorkflowError(
        "Could not locate standardize scripts directory for cleaner imports. Tried: "
        + ", ".join(str(p.resolve()) for p in _candidate_scripts_dirs)
    )
sys.path.insert(0, str(_SCRIPTS_DIR))
from cleaner_registry import CLEANER_REGISTRY, normalize_cleaner_name
from standardize_dependencies import unpack_requires_fetch_done


def _standardize_config_root() -> Path:
    explicit = str(config.get("standardize_config_dir", "")).strip()
    if explicit:
        return Path(explicit)

    profile_dir = str(config.get("profile_config_dir", "")).strip()
    if profile_dir:
        profile_root = Path(profile_dir) / "standardize"
        if profile_root.exists():
            return profile_root

    return Path("config") / "standardize"


def _spec_path_for_study(study: str) -> Path:
    return _standardize_config_root() / f"{study}.yaml"


def _studies_to_process() -> list[str]:
    root = _standardize_config_root()
    if not root.exists():
        return []
    return sorted([p.stem for p in root.glob("*.yaml")])


def _study_out_h5ad(study: str) -> str:
    # Canonical output location is fixed by convention.
    return f"out/{study}/raw.h5ad"


def _study_trim_table(study: str) -> str:
    return f"out/{study}/trimming_tables.csv"


def _study_qc_raw_metrics(study: str) -> str:
    return f"out/{study}/qc_raw/raw_qc_metrics.csv"


def _study_qc_raw_figs(study: str) -> list[str]:
    return [
        f"out/{study}/qc_raw/pictures/raw_log1p_total_counts.png",
        f"out/{study}/qc_raw/pictures/raw_log1p_n_genes_by_counts.png",
        f"out/{study}/qc_raw/pictures/raw_pct_counts_mt.png",
        f"out/{study}/qc_raw/pictures/raw_counts_vs_genes_by_mito.png",
    ]


def _study_mito_gene_prefix(study: str) -> str:
    spec = _load_spec(study)
    inp = spec.get("input") or {}
    return str(inp.get("mito_gene_prefix") or (inp.get("qc") or {}).get("mito_gene_prefix") or "MT-")


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


def _standardize_unpack_fetch_input(study: str) -> str | list[str]:
    # Studies whose raw inputs come from a GEO `*_RAW.tar` are declared in the
    # fetch `extract:` config map; their cleaner inputs only exist after the
    # fetch module's `fetch_extract` rule runs, so depend on its `.extract_done`
    # sentinel (which itself pulls `.fetch_done` -> the download).
    extract_cfg = config.get("extract", {}) or {}
    if str(study) in {str(k) for k in extract_cfg}:
        return f"out/{study}/raw/.extract_done"

    required_paths = _required_spec_input_paths(study)
    if unpack_requires_fetch_done(study=study, required_paths=required_paths):
        return f"out/{study}/.fetch_done"
    return []


def _pipeline_env(study: str) -> str:
    cleaners = _pipeline_cleaners(study)
    needs_r = any(CLEANER_REGISTRY[c].env == "r" for c in cleaners)
    return "../env/standardize_r.yaml" if needs_r else "../env/standardize_py.yaml"


def _workers_for_study(study: str) -> int:
    """Resolve worker count for a study.

    Precedence:
      1) env var SCRNASEQ_QC_STANDARDIZE_WORKERS
      2) spec input.workers
      3) default 1
    """
    env_workers = os.getenv("SCRNASEQ_QC_STANDARDIZE_WORKERS")
    if env_workers is not None and str(env_workers).strip():
        try:
            workers = int(env_workers)
        except ValueError as exc:
            raise WorkflowError("SCRNASEQ_QC_STANDARDIZE_WORKERS must be an integer >= 1") from exc
        if workers < 1:
            raise WorkflowError("SCRNASEQ_QC_STANDARDIZE_WORKERS must be an integer >= 1")
        return workers

    spec = _load_spec(study)
    inp = spec.get("input") or {}
    spec_workers = inp.get("workers")
    if spec_workers is not None:
        try:
            workers = int(spec_workers)
        except (TypeError, ValueError) as exc:
            raise WorkflowError(f"input.workers in config/standardize/{study}.yaml must be an integer >= 1") from exc
        if workers < 1:
            raise WorkflowError(f"input.workers in config/standardize/{study}.yaml must be an integer >= 1")
        return workers

    return 1


def _mem_mb_for_study(study: str) -> int:
    env_mem = os.getenv("SCRNASEQ_QC_STANDARDIZE_MEM_MB")
    if env_mem is not None and str(env_mem).strip():
        try:
            mem_mb = int(env_mem)
        except ValueError as exc:
            raise WorkflowError("SCRNASEQ_QC_STANDARDIZE_MEM_MB must be an integer >= 1") from exc
        if mem_mb < 1:
            raise WorkflowError("SCRNASEQ_QC_STANDARDIZE_MEM_MB must be an integer >= 1")
        return mem_mb

    spec = _load_spec(study)
    inp = spec.get("input") or {}
    spec_mem = inp.get("mem_mb")
    if spec_mem is not None:
        try:
            mem_mb = int(spec_mem)
        except (TypeError, ValueError) as exc:
            raise WorkflowError(f"input.mem_mb in config/standardize/{study}.yaml must be an integer >= 1") from exc
        if mem_mb < 1:
            raise WorkflowError(f"input.mem_mb in config/standardize/{study}.yaml must be an integer >= 1")
        return mem_mb

    return 32000


STUDIES = _studies_to_process()
PY_STUDIES = [study for study in STUDIES if _pipeline_env(study).endswith("standardize_py.yaml")]
R_STUDIES = [study for study in STUDIES if _pipeline_env(study).endswith("standardize_r.yaml")]


def _study_constraint(studies: list[str]) -> str:
    return "|".join(studies) if studies else r"$^"


rule standardize_all:
    input:
        [_study_out_h5ad(s) for s in STUDIES]
        + [_study_trim_table(s) for s in STUDIES]
        + [_study_qc_raw_metrics(s) for s in STUDIES]
        + [fig for s in STUDIES for fig in _study_qc_raw_figs(s)]


rule standardize_unpack_study:
    input:
        spec=lambda wc: str(_spec_path_for_study(wc.study)),
        fetch_done=lambda wc: _standardize_unpack_fetch_input(wc.study),
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


rule standardize_study_py:
    input:
        spec=lambda wc: str(_spec_path_for_study(wc.study)),
        validated_done="out/{study}/.standardize_inputs_validated_done",
    output:
        h5ad=_study_out_h5ad("{study}"),
        trim_csv=_study_trim_table("{study}"),
    wildcard_constraints:
        study=_study_constraint(PY_STUDIES)
    params:
        study=lambda wc: wc.study,
        workers=lambda wc: _workers_for_study(wc.study),
    log:
        "out/logs/standardize/{study}.log"
    threads:
        lambda wc: _workers_for_study(wc.study)
    resources:
        mem_mb=lambda wc: _mem_mb_for_study(wc.study)
    conda:
        "../env/standardize_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        $CONDA_PREFIX/bin/python -u workflow/modules/standardize/scripts/standardize.py \
          --spec {input.spec} \
          --study {params.study} \
          --output {output.h5ad} \
          --workers {params.workers} \
          --log {log}
        """


rule standardize_study_r:
    input:
        spec=lambda wc: str(_spec_path_for_study(wc.study)),
        validated_done="out/{study}/.standardize_inputs_validated_done",
    output:
        h5ad=_study_out_h5ad("{study}"),
        trim_csv=_study_trim_table("{study}"),
    wildcard_constraints:
        study=_study_constraint(R_STUDIES)
    params:
        study=lambda wc: wc.study,
        workers=lambda wc: _workers_for_study(wc.study),
    log:
        "out/logs/standardize/{study}.log"
    threads:
        lambda wc: _workers_for_study(wc.study)
    resources:
        mem_mb=lambda wc: _mem_mb_for_study(wc.study)
    conda:
        "../env/standardize_r.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        $CONDA_PREFIX/bin/python -u workflow/modules/standardize/scripts/standardize.py \
          --spec {input.spec} \
          --study {params.study} \
          --output {output.h5ad} \
          --workers {params.workers} \
          --log {log}
        """


rule standardize_qc_distributions:
    input:
        h5ad=_study_out_h5ad("{study}"),
    output:
        metrics=_study_qc_raw_metrics("{study}"),
        counts_fig=_study_qc_raw_figs("{study}")[0],
        genes_fig=_study_qc_raw_figs("{study}")[1],
        mito_fig=_study_qc_raw_figs("{study}")[2],
        scatter_fig=_study_qc_raw_figs("{study}")[3],
    params:
        study=lambda wc: wc.study,
        mito_gene_prefix=lambda wc: _study_mito_gene_prefix(wc.study),
    log:
        "out/logs/standardize/{study}.qc_distributions.log"
    conda:
        "../env/standardize_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        $CONDA_PREFIX/bin/python -u workflow/modules/standardize/scripts/qc_distribution_plots.py \
          --study {params.study} \
          --input {input.h5ad} \
          --mito-gene-prefix {params.mito_gene_prefix} \
          --output-metrics {output.metrics} \
          --output-counts-fig {output.counts_fig} \
          --output-genes-fig {output.genes_fig} \
          --output-mito-fig {output.mito_fig} \
          --output-scatter-fig {output.scatter_fig} \
          > {log} 2>&1
        """
