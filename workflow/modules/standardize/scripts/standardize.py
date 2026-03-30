import argparse
import json
from pathlib import Path

import yaml
from rtlog import log_kv, log_line


def _load_spec(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Standardize study outputs to a canonical .h5ad")
    parser.add_argument("--spec", required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    study_spec = _load_spec(args.spec)
    process = study_spec.get("process")
    if not process:
        raise SystemExit(f"Missing 'process' in spec: {args.spec}")

    out_path = Path(args.output)
    _ensure_parent(out_path)

    log_path = Path(args.log)
    _ensure_parent(log_path)

    with open(log_path, "w", encoding="utf-8") as log:
        log_line(log, json.dumps({"study": args.study, "process": process, "spec": args.spec}), with_ts=False)
        log_kv(log, "standardize.start", study=args.study, process=process, spec=args.spec)

        if process == "count_matrix_10x":
            from cleaners.count_matrix_10x import run_count_matrix_10x

            run_count_matrix_10x(study=args.study, study_spec=study_spec, output_h5ad=out_path, log=log)
        elif process == "tenx_h5_multi":
            from cleaners.tenx_h5_multi import run_tenx_h5_multi

            run_tenx_h5_multi(study=args.study, study_spec=study_spec, output_h5ad=out_path, log=log)
        elif process == "seurat_rds":
            from cleaners.seurat_rds import run_seurat_rds

            run_seurat_rds(study=args.study, study_spec=study_spec, output_h5ad=out_path, log=log)
        elif process == "gse202379_raw_counts":
            from cleaners.gse202379_raw_counts import run_gse202379_raw_counts

            run_gse202379_raw_counts(
                study=args.study, study_spec=study_spec, output_h5ad=out_path, log=log
            )
        else:
            raise SystemExit(f"Unknown process '{process}' for study '{args.study}'")

        log_kv(log, "standardize.complete", study=args.study, output=args.output)


if __name__ == "__main__":
    main()
