from __future__ import annotations

from pathlib import Path
from typing import TextIO

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.mixture import GaussianMixture
from rtlog import log_kv, log_line


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
    }


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


def _run_miqc_like_filter(a: ad.AnnData, *, qc_cfg: dict, sample_id: str, log: TextIO) -> np.ndarray:
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
        return model_mask

    features = np.column_stack([log_detected[model_mask], mito_pct[model_mask]])
    if not np.isfinite(features).all():
        a.obs["miqc_prob_compromised"] = 1.0
        log_line(log, f"[{sample_id}] miQC-like filter skipped: non-finite model features detected")
        return model_mask

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
    return keep


def _annotate_doublet_scores(a: ad.AnnData, *, qc_cfg: dict, sample_id: str, log: TextIO) -> None:
    score_threshold = qc_cfg["doublet_score_threshold"]
    try:
        sc.pp.scrublet(a, random_state=0, threshold=score_threshold)
    except ValueError as exc:
        a.obs["doublet_score"] = np.nan
        a.obs["predicted_doublet"] = False
        log_line(log, f"[{sample_id}] scrublet annotation skipped: {exc}")
        return

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
        metric_values = a.obs[metric].to_numpy(dtype=float)
        metric_mask = np.isfinite(metric_values) & (metric_values >= lower) & (metric_values <= upper)
        combined_mask &= metric_mask
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
    threshold_report.append(
        (
            f"{mt_metric}: median={mt_median:.3f}, mad={mt_mad:.3f}, "
            f"upper<={mt_upper:.3f} (upper_k={mito_upper_k:.2f}), removed={(~mt_mask).sum()}"
        )
    )

    miqc_keep = _run_miqc_like_filter(a, qc_cfg=qc_cfg, sample_id=sample_id, log=log)
    combined_mask &= miqc_keep

    _annotate_doublet_scores(a, qc_cfg=qc_cfg, sample_id=sample_id, log=log)

    kept = int(combined_mask.sum())
    removed = int((~combined_mask).sum())
    invalid = int((~valid_mask).sum())

    log_line(
        log,
        f"[{sample_id}] QC summary (k={k:.2f}, mito_upper_k={mito_upper_k:.2f}):",
    )
    log_line(log, f"[{sample_id}]   invalid/empty removed={invalid}")
    for line in threshold_report:
        log_line(log, f"[{sample_id}]   {line}")
    log_line(log, f"[{sample_id}] Cells kept={kept} removed={removed}")

    return a[combined_mask].copy()


def run_tenx_h5_multi(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    inp = study_spec.get("input", {})
    samples = _require_list(inp.get("samples"), "input.samples")
    samples_df = _normalize_sample_table(samples)
    qc_cfg = _read_qc_config(study_spec)
    log_kv(log, "tenx_h5_multi.start", study=study, samples=len(samples_df))
    if qc_cfg is None:
        log_line(log, f"[{study}] QC disabled (input.qc missing)")
    else:
        log_line(log, f"[{study}] QC config={qc_cfg}")

    adatas: list[ad.AnnData] = []
    for row in samples_df.itertuples(index=False):
        sample_id = str(row.sample_id)
        h5_val = getattr(row, 'h5', None)
        h5ad_val = getattr(row, 'h5ad', None)
        if isinstance(h5_val, float) and pd.isna(h5_val):
            h5_val = None
        if isinstance(h5ad_val, float) and pd.isna(h5ad_val):
            h5ad_val = None

        if h5_val:
            h5_path = _require_path(str(h5_val), f"samples[{sample_id}].h5")
            log_line(log, f"[{study}] Reading 10x h5: {h5_path} (sample_id={sample_id})")
            a = sc.read_10x_h5(str(h5_path))
        elif h5ad_val:
            h5ad_path = _require_path(str(h5ad_val), f"samples[{sample_id}].h5ad")
            log_line(log, f"[{study}] Reading h5ad: {h5ad_path} (sample_id={sample_id})")
            a = sc.read_h5ad(str(h5ad_path))
        else:
            raise SystemExit(f"Sample '{sample_id}' missing both h5 and h5ad")
        a.var_names_make_unique()
        a.obs["sample_id"] = sample_id
        a.obs["study_accession"] = study
        log_line(log, f"[{study}] Loaded {sample_id}: n_obs={a.n_obs} n_vars={a.n_vars}")

        if qc_cfg is not None:
            before = a.n_obs
            a = _apply_sample_mad_qc(a, sample_id=sample_id, qc_cfg=qc_cfg, log=log)
            log_line(log, f"[{study}] QC-filtered {sample_id}: n_obs={before} -> {a.n_obs}")

        adatas.append(a)

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


def run_tenx_adata_multi(
    *,
    study: str,
    study_spec: dict,
    sample_adatas: list[tuple[str, ad.AnnData]],
    output_h5ad: Path,
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
    for sample_id, a in sample_adatas:
        sample_id = str(sample_id)
        a.var_names_make_unique()
        a.obs["sample_id"] = sample_id
        a.obs["study_accession"] = study
        log_line(log, f"[{study}] Loaded {sample_id}: n_obs={a.n_obs} n_vars={a.n_vars}")

        if qc_cfg is not None:
            before = a.n_obs
            a = _apply_sample_mad_qc(a, sample_id=sample_id, qc_cfg=qc_cfg, log=log)
            log_line(log, f"[{study}] QC-filtered {sample_id}: n_obs={before} -> {a.n_obs}")

        adatas.append(a)

    _finalize_tenx_multi(study=study, inp=inp, adatas=adatas, output_h5ad=output_h5ad, log=log)


def run_tenx_h5_multi(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    inp = study_spec.get("input", {})
    samples = _require_list(inp.get("samples"), "input.samples")
    samples_df = _normalize_sample_table(samples)
    sample_adatas: list[tuple[str, ad.AnnData]] = []
    for row in samples_df.itertuples(index=False):
        sample_id = str(row.sample_id)
        h5_path = _require_path(str(row.h5), f"samples[{sample_id}].h5")
        log_line(log, f"[{study}] Reading 10x h5: {h5_path} (sample_id={sample_id})")
        sample_adatas.append((sample_id, sc.read_10x_h5(str(h5_path))))
    run_tenx_adata_multi(
        study=study,
        study_spec=study_spec,
        sample_adatas=sample_adatas,
        output_h5ad=output_h5ad,
        log=log,
    )
