import os
import re
import sys
from pathlib import Path

import pandas as pd


_project_cfg = os.path.join(workflow.basedir, "config", "fastq_prepros.yaml")
_cwd_project_cfg = os.path.join(os.getcwd(), "config", "fastq_prepros.yaml")
_explicit_cfg = str(config.get("fastq_prepros_config", "")).strip()
_profile_dir = str(config.get("profile_config_dir", "")).strip()
_profile_cfg = os.path.join(_profile_dir, "fastq_prepros.yaml") if _profile_dir else ""

_candidates = [p for p in [_explicit_cfg, _profile_cfg, _project_cfg, _cwd_project_cfg] if p]
_cfg = next((p for p in _candidates if os.path.exists(p)), None)
if _cfg is None:
    raise WorkflowError(
        "Could not find fastq_prepros config. Expected one of: "
        + ", ".join(_candidates)
        + ". If you run this module standalone, pass --configfile config/fastq_prepros.yaml."
    )

configfile: _cfg

if os.path.basename(os.getcwd()) != "a1_p1_qc":
    raise WorkflowError(
        "Run the fastq_prepros module from the a1_p1_qc directory, e.g. `cd a1_p1_qc` "
        "and then `snakemake -s workflow/modules/fastq_prepros/Snakefile ...`."
    )


def _norm(path: str) -> str:
    return os.path.normpath(os.path.expanduser(path))


PROJECT_ROOT = _norm(os.getcwd())


def _rooted(key: str) -> str:
    value = config["paths"][key]
    expanded = os.path.expanduser(value)
    if os.path.isabs(expanded):
        return _norm(expanded)
    return _norm(os.path.join(PROJECT_ROOT, expanded))


def _maybe_relative(key: str) -> str:
    """Return config path normalized, preserving relativeness.

    For output paths (e.g. out_root), keeping a relative path makes
    Snakemake targets portable across machines.
    """
    value = config["paths"][key]
    expanded = os.path.expanduser(value)
    if os.path.isabs(expanded):
        return _norm(expanded)
    return _norm(expanded)


MODULE_DIR = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, MODULE_DIR)

METADATA_CSV = _rooted("metadata_csv")
OUT_ROOT = _maybe_relative("out_root")


FASTQ_PREPROS_CFG = config.get("fastq_prepros", {})
STUDIES = list(FASTQ_PREPROS_CFG.get("studies", []))
if not STUDIES:
    raise WorkflowError("config.fastq_prepros.studies is required (e.g. ['HRA007511'])")

TARGET_LEVEL = str(FASTQ_PREPROS_CFG.get("target_level", "run")).strip().lower()
if TARGET_LEVEL not in {"run", "sample"}:
    raise WorkflowError("config.fastq_prepros.target_level must be one of: run, sample")

ONLY_EXISTING_INPUTS = bool(FASTQ_PREPROS_CFG.get("only_existing_inputs", False))


CELESCOPE_CFG = FASTQ_PREPROS_CFG.get("celescope", {})
CELESCOPE_THREADS = int(CELESCOPE_CFG.get("threads", 16))
CELESCOPE_MEM_GB = int(CELESCOPE_CFG.get("mem_gb", 64))
CELESCOPE_STAR_INDEX = CELESCOPE_CFG.get("star_index")

if not CELESCOPE_STAR_INDEX:
    raise WorkflowError("config.fastq_prepros.celescope.star_index is required")


def _load_metadata() -> pd.DataFrame:
    df = pd.read_csv(METADATA_CSV)
    required = {"include", "study_accession", "platform", "donor_id", "sample_id", "run_id"}
    missing = required - set(df.columns)
    if missing:
        raise WorkflowError(f"Metadata CSV missing columns: {sorted(missing)}")

    df = df[df["include"].astype(str).str.lower().isin({"true", "1", "yes", "y"})]
    df = df[df["study_accession"].isin(STUDIES)]
    return df


def _fastq_dir(study: str, donor_id: str, sample_id: str, run_id: str) -> Path:
    return Path(OUT_ROOT) / study / "fastqs" / donor_id / sample_id / run_id


def _fastq_key_and_read(path: Path) -> tuple[str, int] | None:
    """Extract a pairing key and read number from FASTQ path."""
    name = path.name
    patterns = [
        (r"^(?P<key>.+)_f1\.fastq\.gz$", 1),
        (r"^(?P<key>.+)_r2\.fastq\.gz$", 2),
        (r"^(?P<key>.+)_R1(?:_[0-9]+)?\.fastq\.gz$", 1),
        (r"^(?P<key>.+)_R2(?:_[0-9]+)?\.fastq\.gz$", 2),
        (r"^(?P<key>.+)_1\.fastq\.gz$", 1),
        (r"^(?P<key>.+)_2\.fastq\.gz$", 2),
        (r"^(?P<key>.+)_1\.fq\.gz$", 1),
        (r"^(?P<key>.+)_2\.fq\.gz$", 2),
    ]
    for pattern, read in patterns:
        m = re.match(pattern, name)
        if m:
            return m.group("key"), read
    return None


def _discover_paired_fastqs(fq_dir: Path) -> list[str]:
    """Find paired FASTQs recursively in a run folder.

    Supports files directly under run_id and nested run_id/run_id layouts.
    """
    candidates = sorted(fq_dir.rglob("*.fastq.gz")) + sorted(fq_dir.rglob("*.fq.gz"))
    if not candidates:
        return []

    pairs: dict[str, dict[int, Path]] = {}
    for p in sorted(set(candidates)):
        parsed = _fastq_key_and_read(p)
        if parsed is None:
            continue
        key, read = parsed
        pairs.setdefault(key, {})[read] = p

    resolved: list[str] = []
    for key in sorted(pairs):
        mates = pairs[key]
        if 1 in mates and 2 in mates:
            resolved.extend([str(mates[1]), str(mates[2])])
    return resolved


def _targets(df: pd.DataFrame) -> list[dict]:
    records: list[dict] = []
    for study in STUDIES:
        sub = df[df.study_accession == study]
        if sub.empty:
            continue

        if (sub.platform.astype(str).str.lower() != "singleton").any():
            platforms = sorted(set(sub.platform.astype(str).str.lower()))
            raise WorkflowError(
                f"Study {study} contains non-singleton platform(s): {platforms}. "
                "This module currently implements CeleScope for singleton only."
            )

        if TARGET_LEVEL == "run":
            for _, r in sub.drop_duplicates(["donor_id", "sample_id", "run_id"]).iterrows():
                fq_dir = _fastq_dir(study, r.donor_id, r.sample_id, r.run_id)
                if ONLY_EXISTING_INPUTS and not fq_dir.exists():
                    continue
                records.append(
                    {
                        "study": study,
                        "donor_id": r.donor_id,
                        "sample_id": r.sample_id,
                        "target_id": r.run_id,
                    }
                )
        else:
            for _, r in sub.drop_duplicates(["donor_id", "sample_id"]).iterrows():
                if ONLY_EXISTING_INPUTS:
                    sample_runs = sub[(sub.donor_id == r.donor_id) & (sub.sample_id == r.sample_id)][
                        "run_id"
                    ].tolist()
                    if not any(
                        _fastq_dir(study, r.donor_id, r.sample_id, run_id).exists()
                        for run_id in sample_runs
                    ):
                        continue
                records.append(
                    {
                        "study": study,
                        "donor_id": r.donor_id,
                        "sample_id": r.sample_id,
                        "target_id": r.sample_id,
                    }
                )
    return records


def _input_fastqs(wc):
    df = _load_metadata()
    sub = df[(df.study_accession == wc.study) & (df.sample_id == wc.sample_id)]
    if TARGET_LEVEL == "run":
        sub = sub[sub.run_id == wc.target_id]

    if sub.empty:
        raise WorkflowError(f"No metadata rows for {wc.study} {wc.sample_id} {wc.target_id}")

    fastqs: list[str] = []
    for _, r in sub.drop_duplicates(["donor_id", "sample_id", "run_id"]).iterrows():
        fq_dir = _fastq_dir(wc.study, r.donor_id, r.sample_id, r.run_id)
        if not fq_dir.exists():
            raise WorkflowError(f"FASTQ directory not found: {fq_dir}")

        run_fastqs = _discover_paired_fastqs(fq_dir)
        if not run_fastqs:
            raise WorkflowError(f"Could not find paired FASTQs in: {fq_dir}")
        fastqs.extend(run_fastqs)

    return fastqs


TARGETS = _targets(_load_metadata())


def _all_outputs() -> list[str]:
    return [
        f"{OUT_ROOT}/{t['study']}/counts/singleton/{TARGET_LEVEL}/{t['donor_id']}/{t['sample_id']}/{t['target_id']}/.count_done"
        for t in TARGETS
    ]


rule fastq_prepros_all:
    input:
        _all_outputs()


rule celescope_count:
    input:
        fastqs=_input_fastqs,
    output:
        done=touch(
            f"{OUT_ROOT}/{{study}}/counts/singleton/{TARGET_LEVEL}/{{donor_id}}/{{sample_id}}/{{target_id}}/.count_done"
        ),
    params:
        outdir=lambda wc: str(
            Path(OUT_ROOT)
            / wc.study
            / "counts"
            / "singleton"
            / TARGET_LEVEL
            / wc.donor_id
            / wc.sample_id
            / wc.target_id
        ),
        sample=lambda wc: wc.target_id,
        genome_dir=lambda wc: str(Path(CELESCOPE_STAR_INDEX).resolve().parent),
        threads=CELESCOPE_THREADS,
        fastqs=lambda wc, input: ",".join(input.fastqs),
    log:
        f"{OUT_ROOT}/logs/fastq_prepros/celescope_count/{{study}}_{{donor_id}}_{{sample_id}}_{{target_id}}.log",
    conda:
        "../env/fastq_prepros.yaml"
    threads:
        CELESCOPE_THREADS
    resources:
        mem_mb=CELESCOPE_MEM_GB * 1024,
    shell:
        "python workflow/modules/fastq_prepros/scripts/run_celescope_multi_rna.py "
        "--outdir {params.outdir} "
        "--sample {params.sample} "
        "--fastqs {params.fastqs} "
        "--genome-dir {params.genome_dir} "
        "--threads {params.threads} "
        ">> {log} 2>&1"
