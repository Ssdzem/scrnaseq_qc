import argparse
import json
import os
from collections.abc import Sequence
from pathlib import Path

import yaml
from cleaner_registry import CLEANER_REGISTRY, normalize_cleaner_name
from rtlog import log_kv, log_line


def _load_spec(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _pipeline_cleaners(study_spec: dict) -> list[str]:
    pipeline = study_spec.get("pipeline")
    if pipeline is None:
        pipeline = study_spec.get("Pipeline")
    if pipeline is None:
        raise SystemExit("Missing required 'pipeline' in spec")

    if isinstance(pipeline, Sequence) and not isinstance(pipeline, str):
        cleaners = [normalize_cleaner_name(v) for v in pipeline if str(v).strip()]
    elif isinstance(pipeline, dict):
        # YAML mapping order is preserved.
        cleaners = [normalize_cleaner_name(v) for _k, v in pipeline.items() if str(v).strip()]
    else:
        raise SystemExit("Invalid 'pipeline' format. Use a list or mapping of cleaner names.")

    if not cleaners:
        raise SystemExit("Pipeline is empty")
    unknown = [c for c in cleaners if c not in CLEANER_REGISTRY]
    if unknown:
        raise SystemExit(f"Unknown cleaners in pipeline: {unknown}")
    return cleaners


def _run_cleaner(
    *,
    cleaner: str,
    study: str,
    study_spec: dict,
    output_h5ad: Path,
    workers: int,
    log,
) -> None:
    if cleaner == "count_matrix_10x":
        from cleaners.count_matrix_10x import run_count_matrix_10x

        run_count_matrix_10x(study=study, study_spec=study_spec, output_h5ad=output_h5ad, log=log)
    elif cleaner == "tenx_h5_multi":
        from cleaners.tenx_h5_multi import run_tenx_h5_multi

        run_tenx_h5_multi(
            study=study,
            study_spec=study_spec,
            output_h5ad=output_h5ad,
            log=log,
            workers=workers,
        )
    elif cleaner == "metadata_join":
        from cleaners.metadata_join import run_metadata_join

        run_metadata_join(study=study, study_spec=study_spec, output_h5ad=output_h5ad, log=log)
    elif cleaner == "seurat_rds":
        from cleaners.seurat_rds import run_seurat_rds

        run_seurat_rds(study=study, study_spec=study_spec, output_h5ad=output_h5ad, log=log)
    elif cleaner == "csv_handler":
        from cleaners.csv_handler import run_csv_handler

        run_csv_handler(study=study, study_spec=study_spec, output_h5ad=output_h5ad, log=log)
    elif cleaner == "count_matrix_10x_multi_runs":
        from cleaners.count_matrix_10x_multi_runs import run_count_matrix_10x_multi_runs

        run_count_matrix_10x_multi_runs(
            study=study,
            study_spec=study_spec,
            output_h5ad=output_h5ad,
            workers=workers,
            log=log,
        )
    else:
        raise SystemExit(f"Cleaner '{cleaner}' is not wired in _run_cleaner")


def _run_pipeline(*, study: str, study_spec: dict, output_h5ad: Path, workers: int, log) -> None:
    cleaners = _pipeline_cleaners(study_spec)
    inp = study_spec.setdefault("input", {})
    sample_id = str(inp.get("sample_id", study)).strip() or study

    prev_h5ad: Path | None = None
    metadata_already_joined = False
    prev_cleaner: str | None = None
    for idx, cleaner in enumerate(cleaners):
        is_last = idx == len(cleaners) - 1
        out_path = output_h5ad if is_last else output_h5ad.parent / ".tmp" / f"{study}.pipeline.{idx + 1}.{cleaner}.h5ad"
        out_path.parent.mkdir(parents=True, exist_ok=True)

        log_kv(log, "pipeline.cleaner.start", study=study, cleaner=cleaner, step=idx + 1, output=out_path)

        if cleaner == "count_matrix_10x":
            local_spec = {"input": inp}
        elif cleaner == "metadata_join":
            source_h5ad = prev_h5ad if prev_h5ad is not None else inp.get("h5ad")
            if not source_h5ad:
                raise SystemExit(
                    "Pipeline step 'metadata_join' requires a previous h5ad output or input.h5ad in spec."
                )
            local_spec = {
                "input": {
                    "h5ad": str(source_h5ad),
                    "metadata": inp.get("metadata"),
                    "join_column": inp.get("join_column"),
                    "obs_key": inp.get("obs_key"),
                }
            }
            metadata_already_joined = True
        elif cleaner == "tenx_h5_multi":
            use_explicit_samples = bool(inp.get("samples")) and (
                prev_h5ad is None or prev_cleaner in {"csv_handler", "count_matrix_10x_multi_runs"}
            )
            if use_explicit_samples:
                # Use explicit input.samples for:
                #   1) single-step tenx_h5_multi specs, or
                #   2) csv_handler -> tenx_h5_multi pipelines, where csv_handler writes
                #      a placeholder h5ad and prepares real sample inputs.
                local_spec = {
                    "input": {
                        "samples": inp.get("samples"),
                        "metadata": None if metadata_already_joined else inp.get("metadata"),
                        "qc": inp.get("qc"),
                    }
                }
            else:
                source_h5ad = prev_h5ad if prev_h5ad is not None else inp.get("h5ad")
                if not source_h5ad:
                    raise SystemExit(
                        "Pipeline step 'tenx_h5_multi' requires input.samples or a previous h5ad output (or input.h5ad)."
                    )
                local_spec = {
                    "input": {
                        "samples": [{"sample_id": sample_id, "h5ad": str(source_h5ad)}],
                        "metadata": None if metadata_already_joined else inp.get("metadata"),
                        "qc": inp.get("qc"),
                    }
                }
        else:
            local_spec = study_spec

        _run_cleaner(
            cleaner=cleaner,
            study=study,
            study_spec=local_spec,
            output_h5ad=out_path,
            workers=workers,
            log=log,
        )

        prev_h5ad = out_path
        prev_cleaner = cleaner
        log_kv(log, "pipeline.cleaner.complete", study=study, cleaner=cleaner, step=idx + 1, output=out_path)


def _resolve_workers(*, study_spec: dict, cli_workers: int | None) -> int:
    if cli_workers is not None:
        return max(1, int(cli_workers))

    env_workers = os.getenv("A1_P1_STANDARDIZE_WORKERS")
    if env_workers is not None and str(env_workers).strip():
        try:
            return max(1, int(env_workers))
        except ValueError as exc:
            raise SystemExit("A1_P1_STANDARDIZE_WORKERS must be an integer >= 1") from exc

    inp = study_spec.get("input") or {}
    spec_workers = inp.get("workers")
    if spec_workers is not None:
        try:
            return max(1, int(spec_workers))
        except (TypeError, ValueError) as exc:
            raise SystemExit("input.workers must be an integer >= 1") from exc

    return 1


def main() -> None:
    parser = argparse.ArgumentParser(description="Standardize study outputs to a canonical .h5ad")
    parser.add_argument("--spec", required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    study_spec = _load_spec(args.spec)
    out_path = Path(args.output)
    _ensure_parent(out_path)

    log_path = Path(args.log)
    _ensure_parent(log_path)

    workers = _resolve_workers(study_spec=study_spec, cli_workers=args.workers)

    with open(log_path, "w", encoding="utf-8") as log:
        log_line(log, json.dumps({"study": args.study, "spec": args.spec}), with_ts=False)
        log_kv(log, "standardize.start", study=args.study, spec=args.spec)
        log_kv(log, "standardize.workers", study=args.study, workers=workers)
        _run_pipeline(
            study=args.study,
            study_spec=study_spec,
            output_h5ad=out_path,
            workers=workers,
            log=log,
        )
        log_kv(log, "standardize.complete", study=args.study, output=args.output)


if __name__ == "__main__":
    main()
