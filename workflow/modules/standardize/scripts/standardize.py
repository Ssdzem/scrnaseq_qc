import argparse
import json
from pathlib import Path

import yaml

from cleaners.count_matrix_10x import run_count_matrix_10x


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
        log.write(json.dumps({"study": args.study, "process": process, "spec": args.spec}) + "\n")

        if process == "count_matrix_10x":
            run_count_matrix_10x(study=args.study, study_spec=study_spec, output_h5ad=out_path, log=log)
        else:
            raise SystemExit(f"Unknown process '{process}' for study '{args.study}'")


if __name__ == "__main__":
    main()
