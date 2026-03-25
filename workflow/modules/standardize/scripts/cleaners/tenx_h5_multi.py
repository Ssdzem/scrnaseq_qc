from __future__ import annotations

from pathlib import Path
from typing import TextIO

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.mixture import GaussianMixture


def _log_line(log: TextIO, msg: str) -> None:
    log.write(msg.rstrip() + "\n")
    log.flush()


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
    if "h5" not in df.columns:
        raise SystemExit("Each sample entry must contain 'h5'")
    if df["sample_id"].isna().any():
        raise SystemExit("Some samples are missing 'sample_id'")
    if df["h5"].isna().any():
        raise SystemExit("Some samples are missing 'h5'")
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
    }


def _mad_bounds(values: pd.Series, k: float) -> tuple[float, float, float, float]:
    median = float(values.median())
    mad = float(np.median(np.abs(values.to_numpy() - median)))
    lower = median - (k * mad)
    upper = median + (k * mad)
    return lower, upper, median, mad


def _run_miqc_like_filter(a: ad.AnnData, *, qc_cfg: dict, sample_id: str, log: TextIO) -> np.ndarray:
    """Approximate miQC behavior in Python with a 2-component Gaussian mixture.

    Uses the same core signals as miQC: mitochondrial percentage and library complexity.
    """

    posterior_cutoff = qc_cfg["miqc_posterior_cutoff"]
    detected = a.obs["n_genes_by_counts"].to_numpy(dtype=float)
    mito_pct = a.obs["pct_counts_mt"].to_numpy(dtype=float)
    features = np.column_stack([np.log1p(detected), mito_pct])

    model = GaussianMixture(n_components=2, covariance_type="full", random_state=0)
    model.fit(features)
    posteriors = model.predict_proba(features)
    means = model.means_

    compromised_idx = int(np.argmax(means[:, 1] - means[:, 0]))
    intact_idx = 1 - compromised_idx

    prob_compromised = posteriors[:, compromised_idx]
    a.obs["miqc_prob_compromised"] = prob_compromised

    keep = prob_compromised <= posterior_cutoff

    intact_mean_detected = float(np.expm1(means[intact_idx, 0]))
    intact_mean_mito = float(means[intact_idx, 1])
    keep |= (detected >= intact_mean_detected) & (mito_pct <= intact_mean_mito)

    if (~keep).any():
        min_discard_mito = float(mito_pct[~keep].min())
        left_cutoff_detected = float(detected[np.argmin(np.where(~keep, mito_pct, np.inf))])
        keep |= ~((detected <= left_cutoff_detected) & (mito_pct >= min_discard_mito))

    removed = int((~keep).sum())
    _log_line(
        log,
        f"[{sample_id}] miQC-like filter: posterior_cutoff={posterior_cutoff:.2f}, "
        f"removed={removed}",
    )
    return keep


def _run_doublet_filter(a: ad.AnnData, *, qc_cfg: dict, sample_id: str, log: TextIO) -> np.ndarray:
    score_threshold = qc_cfg["doublet_score_threshold"]
    sc.pp.scrublet(a, random_state=0)

    if "doublet_score" not in a.obs.columns:
        raise SystemExit(f"[{sample_id}] scanpy scrublet did not produce 'doublet_score'")

    keep = a.obs["doublet_score"].to_numpy(dtype=float) <= score_threshold
    removed = int((~keep).sum())
    _log_line(
        log,
        f"[{sample_id}] scrublet filter: score_threshold={score_threshold:.3f}, "
        f"removed={removed}",
    )
    return keep


def _apply_sample_mad_qc(a: ad.AnnData, *, sample_id: str, qc_cfg: dict, log: TextIO) -> ad.AnnData:
    mito_prefix = qc_cfg["mito_gene_prefix"]
    k = qc_cfg["k"]
    mito_upper_k = qc_cfg["mito_upper_k"]

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
        metric_mask = (a.obs[metric] >= lower) & (a.obs[metric] <= upper)
        combined_mask &= metric_mask.to_numpy()
        threshold_report.append(
            (
                f"{metric}: median={median:.3f}, mad={mad:.3f}, "
                f"bounds=[{lower:.3f}, {upper:.3f}], removed={(~metric_mask).sum()}"
            )
        )

    mt_metric = "pct_counts_mt"
    _mt_lower, mt_upper, mt_median, mt_mad = _mad_bounds(a.obs[mt_metric], mito_upper_k)
    mt_mask = a.obs[mt_metric] <= mt_upper
    combined_mask &= mt_mask.to_numpy()
    threshold_report.append(
        (
            f"{mt_metric}: median={mt_median:.3f}, mad={mt_mad:.3f}, "
            f"upper<={mt_upper:.3f} (k={mito_upper_k:.2f}), removed={(~mt_mask).sum()}"
        )
    )

    miqc_keep = _run_miqc_like_filter(a, qc_cfg=qc_cfg, sample_id=sample_id, log=log)
    combined_mask &= miqc_keep

    doublet_keep = _run_doublet_filter(a, qc_cfg=qc_cfg, sample_id=sample_id, log=log)
    combined_mask &= doublet_keep

    kept = int(combined_mask.sum())
    removed = int((~combined_mask).sum())
    invalid = int((~valid_mask).sum())

    _log_line(log, f"[{sample_id}] QC summary (k={k:.2f}, mito_upper_k={mito_upper_k:.2f}):")
    _log_line(log, f"[{sample_id}]   invalid/empty removed={invalid}")
    for line in threshold_report:
        _log_line(log, f"[{sample_id}]   {line}")
    _log_line(log, f"[{sample_id}] Cells kept={kept} removed={removed}")

    return a[combined_mask].copy()


def run_tenx_h5_multi(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    inp = study_spec.get("input", {})
    samples = _require_list(inp.get("samples"), "input.samples")
    samples_df = _normalize_sample_table(samples)
    qc_cfg = _read_qc_config(study_spec)
    _log_line(log, f"[{study}] Method=tenx_h5_multi sample_count={len(samples_df)}")
    if qc_cfg is not None:
        _log_line(log, f"[{study}] QC config={qc_cfg}")

    adatas: list[ad.AnnData] = []
    for row in samples_df.itertuples(index=False):
        sample_id = str(row.sample_id)
        h5_path = _require_path(str(row.h5), f"samples[{sample_id}].h5")
        _log_line(log, f"[{study}] Reading 10x h5: {h5_path} (sample_id={sample_id})")
        a = sc.read_10x_h5(str(h5_path))
        a.var_names_make_unique()
        a.obs["sample_id"] = sample_id
        a.obs["study_accession"] = study
        _log_line(log, f"[{study}] Loaded {sample_id}: n_obs={a.n_obs} n_vars={a.n_vars}")

        if qc_cfg is not None:
            before = a.n_obs
            a = _apply_sample_mad_qc(a, sample_id=sample_id, qc_cfg=qc_cfg, log=log)
            _log_line(log, f"[{study}] QC-filtered {sample_id}: n_obs={before} -> {a.n_obs}")

        adatas.append(a)

    _log_line(log, f"[{study}] Concatenating {len(adatas)} samples")
    merged = ad.concat(adatas, axis=0, join="outer", label=None, merge="same", index_unique="-")
    _log_line(log, f"[{study}] Concatenated: n_obs={merged.n_obs} n_vars={merged.n_vars}")

    metadata_path = inp.get("metadata")
    if metadata_path:
        md = _require_path(metadata_path, "metadata")
        _log_line(log, f"[{study}] Merging metadata: {md}")
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
            _log_line(log, f"[{study}] WARNING: metadata has no obvious barcode column; skipping join")

    merged.uns["standardize"] = {"study": study, "process": "tenx_h5_multi"}
    _log_line(log, f"[{study}] Writing: {output_h5ad}")
    merged.write_h5ad(output_h5ad)
