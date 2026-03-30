from __future__ import annotations

from pathlib import Path
from typing import TextIO

import pandas as pd
from rtlog import log_kv, log_line

from cleaners.count_matrix_10x import run_count_matrix_10x
from cleaners.tenx_h5_multi import run_tenx_h5_multi


def _require_path(value: str | None, name: str) -> Path:
    if not value:
        raise SystemExit(f"Missing required input '{name}'")
    path = Path(value)
    if not path.exists():
        raise SystemExit(f"Input '{name}' not found: {path}")
    return path


def _load_run_table(*, study: str, metadata_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(metadata_csv)
    required = {"study_accession", "donor_id", "sample_id", "run_id", "include"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit(f"Metadata CSV is missing required columns: {missing}")

    study_df = df.loc[df["study_accession"].astype(str) == study].copy()
    if study_df.empty:
        raise SystemExit(f"No metadata rows found for study_accession={study}")

    include_mask = study_df["include"].astype(str).str.upper().eq("TRUE")
    study_df = study_df.loc[include_mask].copy()
    if study_df.empty:
        raise SystemExit(f"No included rows found for study_accession={study}")

    study_df["donor_id"] = study_df["donor_id"].astype(str)
    study_df["sample_id"] = study_df["sample_id"].astype(str)
    study_df["run_id"] = study_df["run_id"].astype(str)
    study_df = study_df.drop_duplicates(subset=["run_id"], keep="first").reset_index(drop=True)
    return study_df


def _run_counts_triplet_dir(*, counts_root: Path, donor_id: str, sample_id: str, run_id: str) -> Path:
    return counts_root / donor_id / sample_id / run_id / "shell" / run_id / "outs" / "filtered"


def _select_runs_for_execution(*, run_df: pd.DataFrame, inp: dict, log: TextIO) -> pd.DataFrame:
    selected = run_df.copy()

    selected_run_ids = inp.get("run_ids")
    if selected_run_ids is not None:
        if not isinstance(selected_run_ids, list) or not selected_run_ids:
            raise SystemExit("input.run_ids must be a non-empty list of run_id values")
        selected_set = {str(x) for x in selected_run_ids}
        selected = selected.loc[selected["run_id"].isin(selected_set)].copy()
        if selected.empty:
            raise SystemExit("input.run_ids did not match any included run_id for this study")
        log_line(log, f"Selected runs by input.run_ids: {sorted(selected_set)}")

    max_runs = inp.get("max_runs")
    if max_runs is not None:
        try:
            max_runs_int = int(max_runs)
        except (TypeError, ValueError) as exc:
            raise SystemExit("input.max_runs must be an integer >= 1") from exc
        if max_runs_int < 1:
            raise SystemExit("input.max_runs must be an integer >= 1")
        selected = selected.head(max_runs_int).copy()
        log_line(log, f"Selected first {len(selected)} runs by input.max_runs={max_runs_int}")

    return selected


def run_hra007511_count_matrix_runs(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    inp = study_spec.get("input", {})
    metadata_csv = _require_path(inp.get("metadata_csv"), "metadata_csv")
    counts_root = _require_path(inp.get("counts_root"), "counts_root")

    run_df = _load_run_table(study=study, metadata_csv=metadata_csv)
    run_df = _select_runs_for_execution(run_df=run_df, inp=inp, log=log)
    log_kv(log, "hra007511_count_matrix_runs.start", study=study, runs=len(run_df))

    tmp_dir = output_h5ad.parent / "_per_run_h5ad"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    per_run_samples: list[dict[str, str]] = []

    for row in run_df.itertuples(index=False):
        triplet_dir = _run_counts_triplet_dir(
            counts_root=counts_root,
            donor_id=row.donor_id,
            sample_id=row.sample_id,
            run_id=row.run_id,
        )
        matrix = triplet_dir / "matrix.mtx.gz"
        genes = triplet_dir / "features.tsv.gz"
        barcodes = triplet_dir / "barcodes.tsv.gz"
        run_h5ad = tmp_dir / f"{row.run_id}.h5ad"

        run_spec = {
            "input": {
                "matrix": str(matrix),
                "genes": str(genes),
                "barcodes": str(barcodes),
            }
        }

        log_line(
            log,
            f"[{study}] Standardizing run_id={row.run_id} donor_id={row.donor_id} "
            f"sample_id={row.sample_id} from {triplet_dir}",
        )
        run_count_matrix_10x(
            study=f"{study}:{row.run_id}",
            study_spec=run_spec,
            output_h5ad=run_h5ad,
            log=log,
        )
        per_run_samples.append({"sample_id": row.run_id, "h5ad": str(run_h5ad)})

    merge_spec = {
        "input": {
            "samples": per_run_samples,
        }
    }
    if "qc" in inp:
        merge_spec["input"]["qc"] = inp["qc"]

    log_line(log, f"[{study}] Merging {len(per_run_samples)} run-level h5ad files")
    run_tenx_h5_multi(study=study, study_spec=merge_spec, output_h5ad=output_h5ad, log=log)
    log_kv(log, "hra007511_count_matrix_runs.complete", study=study, output=output_h5ad)
