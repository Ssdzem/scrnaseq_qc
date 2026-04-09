from __future__ import annotations

from pathlib import Path
from typing import TextIO

import anndata as ad
import pandas as pd
from scipy import sparse
from rtlog import log_kv, log_line

from cleaners.count_matrix_10x import run_count_matrix_10x


def _require_path(value: str | None, name: str) -> Path:
    if not value:
        raise SystemExit(f"Missing required input '{name}'")
    path = Path(value)
    if not path.exists():
        raise SystemExit(f"Input '{name}' not found: {path}")
    return path


def _require_column(df: pd.DataFrame, column: str) -> None:
    if column not in df.columns:
        raise SystemExit(f"runs_csv is missing required column: {column}")


def _load_runs(*, study: str, inp: dict) -> pd.DataFrame:
    runs_csv = _require_path(inp.get("runs_csv"), "runs_csv")
    runs_df = pd.read_csv(runs_csv)

    filters = inp.get("run_filters") or {}
    study_col = str(filters.get("study_col", "study_accession"))
    study_value = str(filters.get("study_value", study))
    include_col = str(filters.get("include_col", "include"))
    include_true_values = filters.get("include_true_values", ["TRUE", "True", "1", "yes", "YES", "Y", "y"])
    donor_id_col = str(filters.get("donor_id_col", "donor_id"))
    sample_id_col = str(filters.get("sample_id_col", "sample_id"))
    run_id_col = str(filters.get("run_id_col", "run_id"))

    for col in [study_col, donor_id_col, sample_id_col, run_id_col, include_col]:
        _require_column(runs_df, col)

    include_set = {str(v).strip() for v in include_true_values}

    study_df = runs_df.loc[runs_df[study_col].astype(str) == study_value].copy()
    if study_df.empty:
        raise SystemExit(f"No rows found in runs_csv for {study_col}={study_value}")

    include_mask = study_df[include_col].astype(str).str.strip().isin(include_set)
    study_df = study_df.loc[include_mask].copy()
    if study_df.empty:
        raise SystemExit(f"No included rows found in runs_csv for {study_col}={study_value}")

    study_df[donor_id_col] = study_df[donor_id_col].astype(str)
    study_df[sample_id_col] = study_df[sample_id_col].astype(str)
    study_df[run_id_col] = study_df[run_id_col].astype(str)

    study_df = study_df.rename(
        columns={
            donor_id_col: "donor_id",
            sample_id_col: "sample_id",
            run_id_col: "run_id",
        }
    )
    return study_df.drop_duplicates(subset=["run_id"], keep="first").reset_index(drop=True)


def _select_runs_for_execution(*, run_df: pd.DataFrame, inp: dict, log: TextIO) -> pd.DataFrame:
    selected = run_df.copy()

    selected_run_ids = inp.get("run_ids")
    if selected_run_ids is not None:
        if not isinstance(selected_run_ids, list) or not selected_run_ids:
            raise SystemExit("input.run_ids must be a non-empty list of run_id values")
        selected_set = {str(x) for x in selected_run_ids}
        selected = selected.loc[selected["run_id"].isin(selected_set)].copy()
        if selected.empty:
            raise SystemExit("input.run_ids did not match any included run_id")
        log_line(log, f"Selected runs by input.run_ids: {sorted(selected_set)}")

    max_runs = inp.get("max_runs")
    run_filters = inp.get("run_filters") or {}
    if max_runs is None:
        max_runs = run_filters.get("max_runs")
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


def _triplet_dir(*, inp: dict, donor_id: str, sample_id: str, run_id: str) -> Path:
    triplet_root = _require_path(inp.get("triplet_root"), "triplet_root")
    triplet_layout = str(
        inp.get(
            "triplet_layout",
            "{triplet_root}/{donor_id}/{sample_id}/{run_id}/shell/{run_id}/outs/filtered",
        )
    )

    resolved = triplet_layout.format(
        triplet_root=str(triplet_root),
        donor_id=donor_id,
        sample_id=sample_id,
        run_id=run_id,
    )
    return Path(resolved)


def run_count_matrix_10x_multi_runs(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    inp = study_spec.setdefault("input", {})

    run_df = _load_runs(study=study, inp=inp)
    run_df = _select_runs_for_execution(run_df=run_df, inp=inp, log=log)
    log_kv(log, "count_matrix_10x_multi_runs.start", study=study, runs=len(run_df))

    tmp_dir = output_h5ad.parent / "_per_run_h5ad"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    per_run_samples: list[dict[str, str]] = []

    for row in run_df.itertuples(index=False):
        triplet_dir = _triplet_dir(
            inp=inp,
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

    inp["samples"] = per_run_samples

    placeholder = ad.AnnData(
        X=sparse.csr_matrix((0, 0), dtype=float),
        obs=pd.DataFrame(index=pd.Index([], name="barcode")),
        var=pd.DataFrame(index=pd.Index([], name="gene_symbol")),
    )
    placeholder.uns["standardize"] = {"study": study, "process": "count_matrix_10x_multi_runs"}
    placeholder.write_h5ad(output_h5ad)
    log_kv(log, "count_matrix_10x_multi_runs.complete", study=study, output=output_h5ad, prepared_samples=len(per_run_samples))
