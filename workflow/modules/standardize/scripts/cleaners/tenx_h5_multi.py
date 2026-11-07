from __future__ import annotations

import json
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, ThreadPoolExecutor, wait
from pathlib import Path
from typing import TextIO

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.mixture import GaussianMixture
from rtlog import log_kv, log_line


TRIMMING_TABLE_COLUMNS = [
    "study_accession",
    "sample_id",
    "qc_applied",
    "cells_before",
    "cells_after",
    "cells_kept",
    "cells_removed",
    "invalid_empty_removed",
    "qc_method",
    "qc_k",
    "qc_mito_upper_k",
    "qc_mito_gene_prefix",
    "miqc_posterior_cutoff",
    "miqc_removed",
    "miqc_skipped",
    "miqc_skip_reason",
    "scrublet_score_threshold",
    "predicted_doublets",
    "scrublet_skipped",
    "scrublet_skip_reason",
    "log1p_total_counts_median",
    "log1p_total_counts_mad",
    "log1p_total_counts_lower_bound",
    "log1p_total_counts_upper_bound",
    "log1p_total_counts_removed",
    "log1p_n_genes_by_counts_median",
    "log1p_n_genes_by_counts_mad",
    "log1p_n_genes_by_counts_lower_bound",
    "log1p_n_genes_by_counts_upper_bound",
    "log1p_n_genes_by_counts_removed",
    "pct_counts_mt_median",
    "pct_counts_mt_mad",
    "pct_counts_mt_upper_bound",
    "pct_counts_mt_upper_k",
    "pct_counts_mt_removed",
    "min_genes",
    "min_genes_removed",
    "max_pct_counts_mt",
    "max_pct_counts_mt_removed",
]


def _require_list(value, name: str) -> list:
    if not value or not isinstance(value, list):
        raise SystemExit(f"Missing required list '{name}'")
    return value


def _require_path(value: str | None, name: str) -> Path:
    if not value:
        raise SystemExit(f"Missing required input '{name}'")
    path = Path(value)
    if not path.exists():
        raise SystemExit(f"Input '{name}' not found: {path}")
    return path


def _normalize_sample_table(samples: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(samples)
    if df.empty:
        raise SystemExit("'samples' list is empty")
    if "sample_id" not in df.columns:
        raise SystemExit("Each sample entry must contain 'sample_id'")
    has_h5 = "h5" in df.columns
    has_h5ad = "h5ad" in df.columns
    if not has_h5 and not has_h5ad:
        raise SystemExit("Each sample entry must contain either 'h5' or 'h5ad'")
    if df["sample_id"].isna().any():
        raise SystemExit("Some samples are missing 'sample_id'")
    if has_h5 and has_h5ad:
        both_missing = df["h5"].isna() & df["h5ad"].isna()
        if both_missing.any():
            raise SystemExit("Some samples are missing both 'h5' and 'h5ad'")
    elif has_h5:
        if df["h5"].isna().any():
            raise SystemExit("Some samples are missing 'h5'")
    elif has_h5ad and df["h5ad"].isna().any():
        raise SystemExit("Some samples are missing 'h5ad'")
    if df["sample_id"].duplicated().any():
        dups = sorted(df.loc[df["sample_id"].duplicated(), "sample_id"].astype(str).unique())
        raise SystemExit(f"Duplicate sample_id values in spec: {dups}")
    return df


def _read_qc_config(study_spec: dict) -> dict | None:
    qc = study_spec.get("input", {}).get("qc")
    if not qc:
        return None
    method = str(qc.get("method", "")).lower()
    if method != "mad":
        raise SystemExit(f"Unsupported qc method '{method}'. Only 'mad' is supported.")

    return {
        "method": method,
        "k": float(qc.get("k", 3.0)),
        "mito_upper_k": float(qc.get("mito_upper_k", 5.0)),
        "mito_gene_prefix": str(qc.get("mito_gene_prefix", "MT-")),
        "miqc_posterior_cutoff": float(qc.get("miqc_posterior_cutoff", 0.75)),
        "doublet_score_threshold": float(qc.get("doublet_score_threshold", 0.25)),
        # Optional absolute floors on top of the per-sample MAD bounds. MAD adapts to
        # each sample, so in a poor sample it accepts what it should remove.
        "min_genes": _optional_number(qc.get("min_genes")),
        "max_pct_counts_mt": _optional_number(qc.get("max_pct_counts_mt")),
    }


def _optional_number(value) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    number = float(value)
    if number < 0:
        raise SystemExit(f"QC limits must be non-negative, got {value!r}")
    return number


def _mad_bounds(values: pd.Series, k: float) -> tuple[float, float, float, float]:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return float("-inf"), float("inf"), float("nan"), float("nan")

    median = float(np.median(finite))
    mad = float(np.median(np.abs(finite - median)))
    lower = median - (k * mad)
    upper = median + (k * mad)
    return lower, upper, median, mad


def _run_miqc_like_filter(
    a: ad.AnnData, *, qc_cfg: dict, sample_id: str, log: TextIO
) -> tuple[np.ndarray, dict[str, float | int | bool | str | None]]:
    """Approximate miQC behavior in Python with a 2-component Gaussian mixture.

    Uses the same core signals as miQC: mitochondrial percentage and library complexity.
    """

    posterior_cutoff = qc_cfg["miqc_posterior_cutoff"]
    detected = a.obs["n_genes_by_counts"].to_numpy(dtype=float)
    mito_pct = a.obs["pct_counts_mt"].to_numpy(dtype=float)

    finite_mask = np.isfinite(detected) & np.isfinite(mito_pct) & (detected >= 0)
    log_detected = np.full(a.n_obs, np.nan, dtype=float)
    log_detected[finite_mask] = np.log1p(detected[finite_mask])
    model_mask = finite_mask & np.isfinite(log_detected)

    if model_mask.sum() < 10:
        a.obs["miqc_prob_compromised"] = 1.0
        log_line(log, f"[{sample_id}] miQC-like filter skipped: <10 usable cells for GMM fit")
        return model_mask, {
            "miqc_posterior_cutoff": posterior_cutoff,
            "miqc_removed": 0,
            "miqc_skipped": True,
            "miqc_skip_reason": "<10 usable cells for GMM fit",
        }

    features = np.column_stack([log_detected[model_mask], mito_pct[model_mask]])
    if not np.isfinite(features).all():
        a.obs["miqc_prob_compromised"] = 1.0
        log_line(log, f"[{sample_id}] miQC-like filter skipped: non-finite model features detected")
        return model_mask, {
            "miqc_posterior_cutoff": posterior_cutoff,
            "miqc_removed": 0,
            "miqc_skipped": True,
            "miqc_skip_reason": "non-finite model features detected",
        }

    model = GaussianMixture(n_components=2, covariance_type="full", random_state=0)
    model.fit(features)
    posteriors = model.predict_proba(features)
    means = model.means_

    compromised_idx = int(np.argmax(means[:, 1] - means[:, 0]))
    intact_idx = 1 - compromised_idx

    prob_compromised = np.ones(a.n_obs, dtype=float)
    prob_compromised[model_mask] = posteriors[:, compromised_idx]
    a.obs["miqc_prob_compromised"] = prob_compromised

    keep = prob_compromised <= posterior_cutoff

    intact_mean_detected = float(np.expm1(means[intact_idx, 0]))
    intact_mean_mito = float(means[intact_idx, 1])
    keep |= model_mask & (detected >= intact_mean_detected) & (mito_pct <= intact_mean_mito)

    discarded_finite = model_mask & (~keep)
    if discarded_finite.any():
        min_discard_mito = float(mito_pct[discarded_finite].min())
        left_cutoff_detected = float(detected[discarded_finite].min())
        keep |= ~((detected <= left_cutoff_detected) & (mito_pct >= min_discard_mito))

    removed = int((~keep).sum())
    log_line(
        log,
        f"[{sample_id}] miQC-like filter: posterior_cutoff={posterior_cutoff:.2f}, "
        f"removed={removed}",
    )
    return keep, {
        "miqc_posterior_cutoff": posterior_cutoff,
        "miqc_removed": removed,
        "miqc_skipped": False,
        "miqc_skip_reason": None,
    }


def _annotate_doublet_scores(
    a: ad.AnnData, *, qc_cfg: dict, sample_id: str, log: TextIO
) -> dict[str, float | int | bool | str | None]:
    score_threshold = qc_cfg["doublet_score_threshold"]
    try:
        sc.pp.scrublet(a, random_state=0, threshold=score_threshold)
    except ValueError as exc:
        a.obs["doublet_score"] = np.nan
        a.obs["predicted_doublet"] = False
        log_line(log, f"[{sample_id}] scrublet annotation skipped: {exc}")
        return {
            "scrublet_score_threshold": score_threshold,
            "predicted_doublets": 0,
            "scrublet_skipped": True,
            "scrublet_skip_reason": str(exc),
        }

    if "doublet_score" not in a.obs.columns:
        a.obs["doublet_score"] = np.nan
    if "predicted_doublet" not in a.obs.columns:
        a.obs["predicted_doublet"] = False
    a.obs["doublet_score"] = pd.to_numeric(a.obs["doublet_score"], errors="coerce").astype(float)
    a.obs["predicted_doublet"] = (
        pd.Series(a.obs["predicted_doublet"]).fillna(False).astype(bool).to_numpy()
    )

    predicted = int(a.obs["predicted_doublet"].to_numpy(dtype=bool).sum())
    log_line(
        log,
        f"[{sample_id}] scrublet annotation: score_threshold={score_threshold:.3f}, "
        f"predicted_doublets={predicted}",
    )
    return {
        "scrublet_score_threshold": score_threshold,
        "predicted_doublets": predicted,
        "scrublet_skipped": False,
        "scrublet_skip_reason": None,
    }


def _apply_sample_mad_qc(
    a: ad.AnnData, *, sample_id: str, qc_cfg: dict, log: TextIO
) -> tuple[ad.AnnData, dict[str, float | int | bool | str | None]]:
    mito_prefix = qc_cfg["mito_gene_prefix"]
    k = qc_cfg["k"]
    mito_upper_k = qc_cfg["mito_upper_k"]
    stats: dict[str, float | int | bool | str | None] = {
        "sample_id": str(sample_id),
        "qc_applied": True,
        "qc_method": str(qc_cfg["method"]),
        "qc_k": float(k),
        "qc_mito_upper_k": float(mito_upper_k),
        "qc_mito_gene_prefix": str(mito_prefix),
        "cells_before": int(a.n_obs),
    }

    a.var["mt"] = a.var_names.str.upper().str.startswith(mito_prefix.upper())
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=True)

    valid_mask = (
        (a.obs["total_counts"] > 0)
        & (a.obs["n_genes_by_counts"] > 0)
        & np.isfinite(a.obs["total_counts"])
        & np.isfinite(a.obs["n_genes_by_counts"])
        & np.isfinite(a.obs["pct_counts_mt"])
    ).to_numpy()

    combined_mask = valid_mask.copy()
    threshold_report: list[str] = []

    for metric in ["log1p_total_counts", "log1p_n_genes_by_counts"]:
        lower, upper, median, mad = _mad_bounds(a.obs[metric], k)
        metric_values = a.obs[metric].to_numpy(dtype=float)
        metric_mask = np.isfinite(metric_values) & (metric_values >= lower) & (metric_values <= upper)
        combined_mask &= metric_mask
        stats[f"{metric}_median"] = median
        stats[f"{metric}_mad"] = mad
        stats[f"{metric}_lower_bound"] = lower
        stats[f"{metric}_upper_bound"] = upper
        stats[f"{metric}_removed"] = int((~metric_mask).sum())
        threshold_report.append(
            (
                f"{metric}: median={median:.3f}, mad={mad:.3f}, "
                f"bounds=[{lower:.3f}, {upper:.3f}], removed={(~metric_mask).sum()}"
            )
        )

    mt_metric = "pct_counts_mt"
    _mt_base_lower, _mt_base_upper, mt_median, mt_mad = _mad_bounds(a.obs[mt_metric], 1.0)
    mt_upper = mt_median + (mito_upper_k * mt_mad)
    mt_values = a.obs[mt_metric].to_numpy(dtype=float)
    mt_mask = np.isfinite(mt_values) & (mt_values <= mt_upper)
    combined_mask &= mt_mask
    stats["pct_counts_mt_median"] = mt_median
    stats["pct_counts_mt_mad"] = mt_mad
    stats["pct_counts_mt_upper_bound"] = mt_upper
    stats["pct_counts_mt_upper_k"] = float(mito_upper_k)
    stats["pct_counts_mt_removed"] = int((~mt_mask).sum())
    threshold_report.append(
        (
            f"{mt_metric}: median={mt_median:.3f}, mad={mt_mad:.3f}, "
            f"upper<={mt_upper:.3f} (upper_k={mito_upper_k:.2f}), removed={(~mt_mask).sum()}"
        )
    )

    for key, metric, keep_if in (
        ("min_genes", "n_genes_by_counts", lambda values, limit: values >= limit),
        ("max_pct_counts_mt", mt_metric, lambda values, limit: values <= limit),
    ):
        limit = qc_cfg.get(key)
        stats[key] = limit
        if limit is None:
            stats[f"{key}_removed"] = 0
            continue
        values = a.obs[metric].to_numpy(dtype=float)
        limit_mask = np.isfinite(values) & keep_if(values, limit)
        combined_mask &= limit_mask
        stats[f"{key}_removed"] = int((~limit_mask).sum())
        threshold_report.append(f"{key}={limit:g}, removed={(~limit_mask).sum()}")

    miqc_keep, miqc_stats = _run_miqc_like_filter(a, qc_cfg=qc_cfg, sample_id=sample_id, log=log)
    combined_mask &= miqc_keep
    stats.update(miqc_stats)

    stats.update(_annotate_doublet_scores(a, qc_cfg=qc_cfg, sample_id=sample_id, log=log))

    kept = int(combined_mask.sum())
    removed = int((~combined_mask).sum())
    invalid = int((~valid_mask).sum())
    stats["cells_after"] = kept
    stats["cells_kept"] = kept
    stats["cells_removed"] = removed
    stats["invalid_empty_removed"] = invalid

    log_line(
        log,
        f"[{sample_id}] QC summary (k={k:.2f}, mito_upper_k={mito_upper_k:.2f}):",
    )
    log_line(log, f"[{sample_id}]   invalid/empty removed={invalid}")
    for line in threshold_report:
        log_line(log, f"[{sample_id}]   {line}")
    log_line(log, f"[{sample_id}] Cells kept={kept} removed={removed}")

    return a[combined_mask].copy(), stats


def _trimming_table_path(output_h5ad: Path) -> Path:
    return output_h5ad.parent / "trimming_tables.csv"


def _write_trimming_table(
    *, study: str, output_h5ad: Path, trim_rows: list[dict], log: TextIO
) -> None:
    trim_path = _trimming_table_path(output_h5ad)
    trim_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(trim_rows)
    if df.empty:
        df = pd.DataFrame(columns=TRIMMING_TABLE_COLUMNS)

    for column in TRIMMING_TABLE_COLUMNS:
        if column not in df.columns:
            df[column] = pd.NA

    extra_columns = sorted([column for column in df.columns if column not in TRIMMING_TABLE_COLUMNS])
    df = df[TRIMMING_TABLE_COLUMNS + extra_columns]

    log_line(log, f"[{study}] Writing trimming table: {trim_path}")
    df.to_csv(trim_path, index=False)
    log_kv(log, "trimming_table.complete", study=study, output=trim_path, rows=len(df))


def _finalize_tenx_multi(
    *,
    study: str,
    inp: dict,
    adatas: list[ad.AnnData],
    output_h5ad: Path,
    log: TextIO,
) -> None:
    if not adatas:
        raise SystemExit(f"[{study}] No sample adatas were provided")

    log_line(log, f"[{study}] Concatenating {len(adatas)} samples")
    merged = ad.concat(adatas, axis=0, join="outer", label=None, merge="same", index_unique="-")
    log_line(log, f"[{study}] Concatenated: n_obs={merged.n_obs} n_vars={merged.n_vars}")

    if "doublet_score" in merged.obs.columns:
        merged.obs["doublet_score"] = pd.to_numeric(merged.obs["doublet_score"], errors="coerce").astype(float)
    if "predicted_doublet" in merged.obs.columns:
        merged.obs["predicted_doublet"] = (
            pd.Series(merged.obs["predicted_doublet"]).fillna(False).astype(bool).to_numpy()
        )

    metadata_path = inp.get("metadata")
    if metadata_path:
        md = _require_path(metadata_path, "metadata")
        log_line(log, f"[{study}] Merging metadata: {md}")
        md_df = pd.read_csv(md)
        if "sample_id" in md_df.columns:
            md_df["sample_id"] = md_df["sample_id"].astype(str)
        key = None
        for candidate in ["cell", "barcode", "Cell", "Barcode"]:
            if candidate in md_df.columns:
                key = candidate
                break
        if key is not None:
            md_df = md_df.set_index(key)
            md_df.index = md_df.index.astype(str)
            merged.obs = merged.obs.join(md_df, how="left")
        else:
            log_line(log, f"[{study}] WARNING: metadata has no obvious barcode column; skipping join")

    merged.uns["standardize"] = {"study": study, "process": "tenx_h5_multi"}
    log_line(log, f"[{study}] Writing: {output_h5ad}")
    merged.write_h5ad(output_h5ad)
    log_kv(log, "tenx_h5_multi.complete", study=study, output=output_h5ad)


def _is_missing_like(value) -> bool:
    if pd.isna(value):
        return True
    if isinstance(value, str):
        return not value.strip()
    return False


def _fill_obs_column_if_missing(
    *,
    adata: ad.AnnData,
    column: str,
    fill_value: str,
    study: str,
    log: TextIO,
) -> None:
    fill_value = str(fill_value)
    if column not in adata.obs.columns:
        adata.obs[column] = fill_value
        log_line(log, f"[{study}] Filled missing obs column '{column}' with fallback value '{fill_value}'")
        return

    col = adata.obs[column]
    mask = col.map(_is_missing_like)
    missing_count = int(mask.sum())
    if missing_count == 0:
        log_line(log, f"[{study}] Preserved existing obs column '{column}' (no missing values)")
        return

    col = col.astype(object)
    col.loc[mask] = fill_value
    adata.obs[column] = col
    log_line(log, f"[{study}] Preserved '{column}' and filled {missing_count} missing values with '{fill_value}'")


def _resolve_workers(*, workers: int, n_samples: int) -> int:
    if n_samples <= 1:
        return 1
    return max(1, min(int(workers), int(n_samples)))


def _prepare_loaded_sample(
    *,
    study: str,
    sample_id: str,
    a: ad.AnnData,
    qc_cfg: dict | None,
) -> tuple[str, ad.AnnData, list[str], dict[str, float | int | bool | str | None]]:
    messages: list[str] = []
    sample_id = str(sample_id)
    a.var_names_make_unique()
    _fill_obs_column_if_missing_in_memory(adata=a, column="sample_id", fill_value=sample_id, messages=messages)
    _fill_obs_column_if_missing_in_memory(
        adata=a, column="study_accession", fill_value=study, messages=messages
    )
    messages.append(f"Loaded {sample_id}: n_obs={a.n_obs} n_vars={a.n_vars}")
    trim_stats: dict[str, float | int | bool | str | None] = {
        "study_accession": str(study),
        "sample_id": sample_id,
        "qc_applied": False,
        "cells_before": int(a.n_obs),
        "cells_after": int(a.n_obs),
        "cells_kept": int(a.n_obs),
        "cells_removed": 0,
        "invalid_empty_removed": 0,
    }

    if qc_cfg is not None:
        before = a.n_obs
        # QC implementation is shared and deterministic per sample.
        # It mutates adata in-place and then returns filtered copy.
        buffer = _ListLog()
        a, qc_stats = _apply_sample_mad_qc(a, sample_id=sample_id, qc_cfg=qc_cfg, log=buffer)
        messages.extend(buffer.lines)
        messages.append(f"QC-filtered {sample_id}: n_obs={before} -> {a.n_obs}")
        trim_stats.update(qc_stats)

    return sample_id, a, messages, trim_stats


def _fill_obs_column_if_missing_in_memory(
    *,
    adata: ad.AnnData,
    column: str,
    fill_value: str,
    messages: list[str],
) -> None:
    fill_value = str(fill_value)
    if column not in adata.obs.columns:
        adata.obs[column] = fill_value
        messages.append(f"Filled missing obs column '{column}' with fallback value '{fill_value}'")
        return

    col = adata.obs[column]
    mask = col.map(_is_missing_like)
    missing_count = int(mask.sum())
    if missing_count == 0:
        messages.append(f"Preserved existing obs column '{column}' (no missing values)")
        return

    col = col.astype(object)
    col.loc[mask] = fill_value
    adata.obs[column] = col
    messages.append(f"Preserved '{column}' and filled {missing_count} missing values with '{fill_value}'")


class _ListLog:
    def __init__(self) -> None:
        self.lines: list[str] = []

    def write(self, text: str) -> int:
        if text:
            for line in str(text).splitlines():
                if line.strip():
                    self.lines.append(line)
        return len(text)

    def flush(self) -> None:
        return


def _prepare_sample_from_source_task(
    task: dict,
) -> tuple[int, str, str, dict[str, float | int | bool | str | None], list[str]]:
    idx = int(task["idx"])
    study = str(task["study"])
    sample_id = str(task["sample_id"])
    source_kind = str(task["source_kind"])
    source_path = str(task["source_path"])
    qc_cfg = task.get("qc_cfg")
    out_h5ad = Path(str(task["out_h5ad"]))
    out_trim_json = Path(str(task["out_trim_json"]))

    if source_kind == "h5":
        a = sc.read_10x_h5(source_path)
    elif source_kind == "h5ad":
        a = sc.read_h5ad(source_path)
    else:
        raise SystemExit(f"Unsupported source_kind for sample '{sample_id}': {source_kind}")

    sample_id_out, prepared, messages, trim_stats = _prepare_loaded_sample(
        study=study,
        sample_id=sample_id,
        a=a,
        qc_cfg=qc_cfg,
    )
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    out_trim_json.parent.mkdir(parents=True, exist_ok=True)
    prepared.write_h5ad(out_h5ad)
    out_trim_json.write_text(json.dumps(trim_stats, indent=2, sort_keys=True), encoding="utf-8")
    messages.append(f"Wrote sample checkpoint h5ad: {out_h5ad}")
    messages.append(f"Wrote sample trimming stats: {out_trim_json}")
    return idx, sample_id_out, str(out_h5ad), trim_stats, messages


def run_tenx_h5_multi(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO, workers: int = 1) -> None:
    inp = study_spec.get("input", {})
    samples = _require_list(inp.get("samples"), "input.samples")
    samples_df = _normalize_sample_table(samples)
    qc_cfg = _read_qc_config(study_spec)
    log_kv(log, "tenx_h5_multi.start", study=study, samples=len(samples_df))
    if qc_cfg is None:
        log_line(log, f"[{study}] QC disabled (input.qc missing)")
    else:
        log_line(log, f"[{study}] QC config={qc_cfg}")

    checkpoint_dir = output_h5ad.parent / ".tmp" / f"{study}.tenx_qc_samples"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    sample_checkpoints: dict[int, tuple[str, Path, list[str], dict[str, float | int | bool | str | None]]] = {}
    tasks: list[dict] = []
    effective_workers = _resolve_workers(workers=workers, n_samples=len(samples_df))
    log_line(log, f"[{study}] tenx_h5_multi worker processes={effective_workers}")
    for idx, row in enumerate(samples_df.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        h5_val = getattr(row, "h5", None)
        h5ad_val = getattr(row, "h5ad", None)
        if isinstance(h5_val, float) and pd.isna(h5_val):
            h5_val = None
        if isinstance(h5ad_val, float) and pd.isna(h5ad_val):
            h5ad_val = None

        checkpoint_h5ad = checkpoint_dir / f"{idx:03d}.{sample_id}.qc.h5ad"
        checkpoint_trim_json = checkpoint_dir / f"{idx:03d}.{sample_id}.trim.json"
        if (
            checkpoint_h5ad.exists()
            and checkpoint_h5ad.stat().st_size > 0
            and checkpoint_trim_json.exists()
            and checkpoint_trim_json.stat().st_size > 0
        ):
            log_line(log, f"[{study}] Reusing sample checkpoint: {checkpoint_h5ad} (sample_id={sample_id})")
            trim_stats = json.loads(checkpoint_trim_json.read_text(encoding="utf-8"))
            sample_checkpoints[idx] = (sample_id, checkpoint_h5ad, [], trim_stats)
            continue

        if h5_val:
            h5_path = _require_path(str(h5_val), f"samples[{sample_id}].h5")
            tasks.append(
                {
                    "idx": idx,
                    "study": study,
                    "sample_id": sample_id,
                    "source_kind": "h5",
                    "source_path": str(h5_path),
                    "qc_cfg": qc_cfg,
                    "out_h5ad": str(checkpoint_h5ad),
                    "out_trim_json": str(checkpoint_trim_json),
                }
            )
        elif h5ad_val:
            h5ad_path = _require_path(str(h5ad_val), f"samples[{sample_id}].h5ad")
            tasks.append(
                {
                    "idx": idx,
                    "study": study,
                    "sample_id": sample_id,
                    "source_kind": "h5ad",
                    "source_path": str(h5ad_path),
                    "qc_cfg": qc_cfg,
                    "out_h5ad": str(checkpoint_h5ad),
                    "out_trim_json": str(checkpoint_trim_json),
                }
            )
        else:
            raise SystemExit(f"Sample '{sample_id}' missing both h5 and h5ad")

    if tasks:
        if effective_workers == 1:
            total = len(tasks)
            for i, task in enumerate(tasks, start=1):
                log_line(log, f"[{study}] tenx_h5_multi progress queued={i}/{total} sample_id={task['sample_id']}")
                idx, sample_id_out, out_h5ad_path, trim_stats, messages = _prepare_sample_from_source_task(task)
                log_line(log, f"[{study}] tenx_h5_multi progress completed={i}/{total} sample_id={sample_id_out}")
                sample_checkpoints[idx] = (sample_id_out, Path(out_h5ad_path), messages, trim_stats)
        else:
            total = len(tasks)
            with ProcessPoolExecutor(max_workers=effective_workers) as pool:
                futures = {}
                for i, task in enumerate(tasks, start=1):
                    log_line(log, f"[{study}] tenx_h5_multi progress queued={i}/{total} sample_id={task['sample_id']}")
                    fut = pool.submit(_prepare_sample_from_source_task, task)
                    futures[fut] = task

                completed = 0
                pending = set(futures)
                while pending:
                    done_now, pending = wait(pending, timeout=30, return_when=FIRST_COMPLETED)
                    if not done_now:
                        log_line(
                            log,
                            f"[{study}] tenx_h5_multi heartbeat completed={completed}/{total} pending={len(pending)}",
                        )
                        continue
                    for fut in done_now:
                        idx, sample_id_out, out_h5ad_path, trim_stats, messages = fut.result()
                        completed += 1
                        log_line(
                            log,
                            f"[{study}] tenx_h5_multi progress completed={completed}/{total} sample_id={sample_id_out}",
                        )
                        sample_checkpoints[idx] = (sample_id_out, Path(out_h5ad_path), messages, trim_stats)

    adatas: list[ad.AnnData] = []
    trim_rows: list[dict[str, float | int | bool | str | None]] = []
    for idx in sorted(sample_checkpoints):
        sample_id, path, messages, trim_stats = sample_checkpoints[idx]
        for msg in messages:
            log_line(log, f"[{study}] {msg}")
        log_line(log, f"[{study}] Loading sample checkpoint: {path} (sample_id={sample_id})")
        adatas.append(sc.read_h5ad(str(path)))
        trim_rows.append(trim_stats)

    _finalize_tenx_multi(study=study, inp=inp, adatas=adatas, output_h5ad=output_h5ad, log=log)
    _write_trimming_table(study=study, output_h5ad=output_h5ad, trim_rows=trim_rows, log=log)


def run_tenx_adata_multi(
    *,
    study: str,
    study_spec: dict,
    sample_adatas: list[tuple[str, ad.AnnData]],
    output_h5ad: Path,
    workers: int = 1,
    log: TextIO,
) -> None:
    inp = study_spec.get("input", {})
    qc_cfg = _read_qc_config(study_spec)
    log_kv(log, "tenx_h5_multi.start", study=study, samples=len(sample_adatas))
    if qc_cfg is None:
        log_line(log, f"[{study}] QC disabled (input.qc missing)")
    else:
        log_line(log, f"[{study}] QC config={qc_cfg}")

    adatas: list[ad.AnnData] = []
    trim_rows: list[dict[str, float | int | bool | str | None]] = []
    effective_workers = _resolve_workers(workers=workers, n_samples=len(sample_adatas))
    log_line(log, f"[{study}] tenx_adata_multi worker threads={effective_workers}")

    if effective_workers == 1:
        total = len(sample_adatas)
        for idx, (sample_id, a) in enumerate(sample_adatas, start=1):
            log_line(log, f"[{study}] tenx_adata_multi progress queued={idx}/{total} sample_id={sample_id}")
            sample_id_out, processed, messages, trim_stats = _prepare_loaded_sample(
                study=study,
                sample_id=str(sample_id),
                a=a,
                qc_cfg=qc_cfg,
            )
            log_line(log, f"[{study}] tenx_adata_multi progress completed={idx}/{total} sample_id={sample_id_out}")
            for msg in messages:
                log_line(log, f"[{study}] {msg}")
            adatas.append(processed)
            trim_rows.append(trim_stats)
    else:
        total = len(sample_adatas)
        results_by_index: dict[
            int, tuple[ad.AnnData, list[str], dict[str, float | int | bool | str | None]]
        ] = {}
        with ThreadPoolExecutor(max_workers=effective_workers) as pool:
            futures = {}
            for idx, item in enumerate(sample_adatas, start=1):
                sample_id = str(item[0])
                log_line(log, f"[{study}] tenx_adata_multi progress queued={idx}/{total} sample_id={sample_id}")
                fut = pool.submit(
                    _prepare_loaded_sample,
                    study=study,
                    sample_id=sample_id,
                    a=item[1],
                    qc_cfg=qc_cfg,
                )
                futures[fut] = idx

            completed = 0
            pending = set(futures)
            while pending:
                done_now, pending = wait(pending, timeout=30, return_when=FIRST_COMPLETED)
                if not done_now:
                    log_line(
                        log,
                        f"[{study}] tenx_adata_multi heartbeat completed={completed}/{total} pending={len(pending)}",
                    )
                    continue
                for fut in done_now:
                    idx = futures[fut]
                    sample_id_out, processed, messages, trim_stats = fut.result()
                    completed += 1
                    log_line(
                        log,
                        f"[{study}] tenx_adata_multi progress completed={completed}/{total} sample_id={sample_id_out}",
                    )
                    results_by_index[idx] = (processed, messages, trim_stats)

        for idx in sorted(results_by_index):
            processed, messages, trim_stats = results_by_index[idx]
            for msg in messages:
                log_line(log, f"[{study}] {msg}")
            adatas.append(processed)
            trim_rows.append(trim_stats)

    _finalize_tenx_multi(study=study, inp=inp, adatas=adatas, output_h5ad=output_h5ad, log=log)
    _write_trimming_table(study=study, output_h5ad=output_h5ad, trim_rows=trim_rows, log=log)
