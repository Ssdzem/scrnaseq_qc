import argparse
import json
from collections.abc import Sequence
from pathlib import Path

import yaml
from rtlog import log_kv, log_line


def _load_spec(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _normalize_cleaner_name(value: str) -> str:
    cleaner = str(value).strip()
    if cleaner.endswith(".py"):
        cleaner = cleaner[:-3]
    return cleaner


def _pipeline_cleaners(study_spec: dict) -> list[str]:
    pipeline = study_spec.get("pipeline")
    if pipeline is None:
        pipeline = study_spec.get("Pipeline")
    if pipeline is None:
        return []

    if isinstance(pipeline, Sequence) and not isinstance(pipeline, str):
        return [_normalize_cleaner_name(v) for v in pipeline if str(v).strip()]

    if isinstance(pipeline, dict):
        # Preserve explicit cleaner ordering from keys like cleaner1/cleaner2.
        return [_normalize_cleaner_name(v) for _k, v in pipeline.items() if str(v).strip()]

    raise SystemExit("Invalid 'pipeline' format. Use a list or mapping of cleaner names.")


def _run_process(
    *,
    process: str,
    study: str,
    study_spec: dict,
    output_h5ad: Path,
    log,
) -> None:
    if process == "count_matrix_10x":
        from cleaners.count_matrix_10x import run_count_matrix_10x

        run_count_matrix_10x(study=study, study_spec=study_spec, output_h5ad=output_h5ad, log=log)
    elif process == "tenx_h5_multi":
        from cleaners.tenx_h5_multi import run_tenx_h5_multi

        run_tenx_h5_multi(study=study, study_spec=study_spec, output_h5ad=output_h5ad, log=log)
    elif process == "seurat_rds":
        from cleaners.seurat_rds import run_seurat_rds

        run_seurat_rds(study=study, study_spec=study_spec, output_h5ad=output_h5ad, log=log)
    elif process == "gse202379_raw_counts":
        from cleaners.gse202379_raw_counts import run_gse202379_raw_counts

        run_gse202379_raw_counts(study=study, study_spec=study_spec, output_h5ad=output_h5ad, log=log)
    elif process == "hra007511_count_matrix_runs":
        from cleaners.hra007511_count_matrix_runs import run_hra007511_count_matrix_runs

        run_hra007511_count_matrix_runs(
            study=study,
            study_spec=study_spec,
            output_h5ad=output_h5ad,
            log=log,
        )
    else:
        raise SystemExit(f"Unknown process '{process}' for study '{study}'")


def _run_pipeline(*, study: str, study_spec: dict, output_h5ad: Path, log) -> None:
    cleaners = _pipeline_cleaners(study_spec)
    if not cleaners:
        raise SystemExit("Missing both 'process' and 'pipeline' in spec")

    inp = dict(study_spec.get("input") or {})
    sample_id = str(inp.get("sample_id", study)).strip() or study
    qc_cfg = inp.get("qc")

    prev_h5ad: Path | None = None
    for idx, cleaner in enumerate(cleaners):
        is_last = idx == len(cleaners) - 1
        out_path = output_h5ad if is_last else output_h5ad.parent / ".tmp" / f"{study}.pipeline.{idx + 1}.{cleaner}.h5ad"
        out_path.parent.mkdir(parents=True, exist_ok=True)

        log_kv(log, "pipeline.cleaner.start", study=study, cleaner=cleaner, step=idx + 1, output=out_path)

        if cleaner == "count_matrix_10x":
            local_spec = {"input": inp}
            _run_process(
                process="count_matrix_10x",
                study=study,
                study_spec=local_spec,
                output_h5ad=out_path,
                log=log,
            )
        elif cleaner == "tenx_h5_multi":
            source_h5ad = prev_h5ad if prev_h5ad is not None else inp.get("h5ad")
            if not source_h5ad:
                raise SystemExit(
                    "Pipeline step 'tenx_h5_multi' requires a previous h5ad output "
                    "or input.h5ad in spec."
                )
            local_spec = {
                "input": {
                    "samples": [{"sample_id": sample_id, "h5ad": str(source_h5ad)}],
                    "metadata": inp.get("metadata"),
                    "qc": qc_cfg,
                }
            }
            _run_process(
                process="tenx_h5_multi",
                study=study,
                study_spec=local_spec,
                output_h5ad=out_path,
                log=log,
            )
        else:
            raise SystemExit(
                f"Unsupported cleaner '{cleaner}' in pipeline for study '{study}'. "
                "Currently supported in pipeline: count_matrix_10x, tenx_h5_multi"
            )

        prev_h5ad = out_path
        log_kv(log, "pipeline.cleaner.complete", study=study, cleaner=cleaner, step=idx + 1, output=out_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Standardize study outputs to a canonical .h5ad")
    parser.add_argument("--spec", required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    study_spec = _load_spec(args.spec)
    process = study_spec.get("process")

    out_path = Path(args.output)
    _ensure_parent(out_path)

    log_path = Path(args.log)
    _ensure_parent(log_path)

    with open(log_path, "w", encoding="utf-8") as log:
        log_line(log, json.dumps({"study": args.study, "process": process, "spec": args.spec}), with_ts=False)
        log_kv(log, "standardize.start", study=args.study, process=process, spec=args.spec)

        if process:
            _run_process(
                process=str(process),
                study=args.study,
                study_spec=study_spec,
                output_h5ad=out_path,
                log=log,
            )
        else:
            _run_pipeline(study=args.study, study_spec=study_spec, output_h5ad=out_path, log=log)

        log_kv(log, "standardize.complete", study=args.study, output=args.output)


if __name__ == "__main__":
    main()
