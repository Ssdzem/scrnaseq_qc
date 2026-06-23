#!/usr/bin/env python3
"""Recreate the ordered a1_p1_qc figure set from pipeline objects.

The historical files under ``out/figures`` were Marimo/interactive ``annotate-*``
PNGs.  This script gives those figures stable, semantic filenames and records
which old file each generated image replaces in ``manifest.json``.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence


@dataclass(frozen=True)
class FigureSpec:
    order: int
    filename: str
    legacy_reference: str
    title: str
    source_object: str
    structure: str
    represents: str
    pipeline_position: str


FIGURES: tuple[FigureSpec, ...] = (
    FigureSpec(
        order=1,
        filename="01_umap_study_accession.png",
        legacy_reference="annotate-2026-06-22_20-12-51.png",
        title="study_accession UMAP",
        source_object="integration_annotation H5AD: .obs['study_accession'] + .obsm['X_umap']",
        structure="cell-level AnnData object with UMAP coordinates and study metadata",
        represents="Integrated myeloid/liver cell embedding colored by source study accession.",
        pipeline_position="After integration/annotation has produced the all-study annotated object.",
    ),
    FigureSpec(
        order=2,
        filename="02_unique_samples_by_fibrosis_stage.png",
        legacy_reference="annotate-2026-06-22_20-13-36.png",
        title="Unique samples by fibrosis stage",
        source_object="configured metadata CSV: config/liver_scRNAseq_metadata_everything.csv filtered to include truthy rows",
        structure="included sample-level metadata table grouped by fibrosis_stage",
        represents="Number of included distinct samples per fibrosis stage.",
        pipeline_position="Metadata summary / cohort composition diagnostic.",
    ),
    FigureSpec(
        order=3,
        filename="03_unique_samples_by_study_accession.png",
        legacy_reference="annotate-2026-06-22_20-13-55.png",
        title="Unique samples by study accession",
        source_object="configured metadata CSV: config/liver_scRNAseq_metadata_everything.csv filtered to include truthy rows",
        structure="included sample-level metadata table grouped by study_accession",
        represents="Number of included distinct samples contributed by each study accession.",
        pipeline_position="Metadata summary / cohort composition diagnostic.",
    ),
    FigureSpec(
        order=4,
        filename="04_cell_type_study_composition.png",
        legacy_reference="annotate-2026-06-22_20-14-33.png",
        title="Cell-type composition by study accession",
        source_object="integration_annotation H5AD: cell-type label column + .obs['study_accession']",
        structure="cell-level AnnData obs cross-tabulated by label and study",
        represents="Within-cell-type study contribution proportions, exposing study-specific composition biases.",
        pipeline_position="After annotation labels are assigned to the integrated object.",
    ),
)


def _import_plotting():
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    return np, pd, plt


def _read_h5ad_backed(path: Path):
    import anndata as ad

    if not path.exists():
        raise FileNotFoundError(f"H5AD source not found: {path}")
    return ad.read_h5ad(path, backed="r")


def _as_bool_series(series):
    lowered = series.astype("string").str.strip().str.lower()
    return lowered.isin({"true", "1", "yes", "y", "t"})


def _choose_column(columns: Iterable[str], candidates: Sequence[str], purpose: str) -> str:
    available = set(columns)
    for candidate in candidates:
        if candidate in available:
            return candidate
    raise KeyError(f"Could not find a {purpose} column. Tried: {', '.join(candidates)}")


def _close_h5ad(adata) -> None:
    file_obj = getattr(adata, "file", None)
    if file_obj is not None:
        close = getattr(file_obj, "close", None)
        if callable(close):
            close()


def _load_obs_umap(h5ad_path: Path):
    np, pd, _plt = _import_plotting()
    adata = _read_h5ad_backed(h5ad_path)
    try:
        obs = adata.obs.copy()
        if "X_umap" not in adata.obsm:
            raise KeyError(f"{h5ad_path} does not contain .obsm['X_umap']")
        umap = np.asarray(adata.obsm["X_umap"])
        if umap.shape[1] < 2:
            raise ValueError(f"{h5ad_path} .obsm['X_umap'] must have at least two columns")
        obs = obs.reset_index().rename(columns={"index": "_index"})
        return obs, umap[:, :2]
    finally:
        _close_h5ad(adata)


def _load_metadata_counts(metadata_csv: Path, group_col: str, pd):
    if not metadata_csv.exists():
        raise FileNotFoundError(f"Metadata source not found: {metadata_csv}")
    meta = pd.read_csv(metadata_csv)
    if "sample_id" not in meta.columns or group_col not in meta.columns or "include" not in meta.columns:
        raise KeyError(f"{metadata_csv} must contain sample_id, include, and {group_col}")
    meta = meta.loc[_as_bool_series(meta["include"])].copy()
    if meta.empty:
        raise ValueError(f"{metadata_csv} has no include-truthy rows for {group_col}")
    counts = (
        meta.dropna(subset=[group_col, "sample_id"])
        .drop_duplicates([group_col, "sample_id"])
        .groupby(group_col, observed=True)["sample_id"]
        .nunique()
        .rename("unique_sample_count")
        .reset_index()
    )
    return counts


def _save_bar_with_labels(df, *, x_col: str, y_col: str, title: str, xlabel: str, ylabel: str, out: Path, color: str):
    _np, _pd, plt = _import_plotting()
    plot_df = df.copy()
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
    bars = ax.bar(plot_df[x_col].astype(str), plot_df[y_col], color=color, edgecolor="black", linewidth=0.4)
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    ax.set_axisbelow(True)
    max_val = float(plot_df[y_col].max()) if len(plot_df) else 1.0
    ax.set_ylim(0, max_val * 1.18 if max_val else 1.0)
    for bar, value in zip(bars, plot_df[y_col]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max_val * 0.02,
            f"{int(value):,}",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )
    if len(plot_df) > 3:
        ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_umap_study_accession(args: argparse.Namespace, spec: FigureSpec) -> Path:
    np, pd, plt = _import_plotting()
    obs, umap = _load_obs_umap(args.h5ad)
    study_col = _choose_column(obs.columns, ["study_accession", "study", "dataset"], "study accession")
    categories = sorted(obs[study_col].astype("string").fillna("NA").unique())
    palette = plt.get_cmap("tab10")

    fig, ax = plt.subplots(figsize=(7.1, 4.5), dpi=150)
    for i, category in enumerate(categories):
        mask = obs[study_col].astype("string").fillna("NA").eq(category).to_numpy()
        ax.scatter(
            umap[mask, 0],
            umap[mask, 1],
            s=1.0,
            alpha=0.45,
            linewidths=0,
            color=palette(i % 10),
            label=str(category),
            rasterized=True,
        )
    ax.set_title("study_accession")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=5, fontsize=8)
    fig.tight_layout()
    out = args.out_dir / spec.filename
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_unique_samples_by_fibrosis(args: argparse.Namespace, spec: FigureSpec) -> Path:
    _np, pd, _plt = _import_plotting()
    df = _load_metadata_counts(args.metadata_csv, "fibrosis_stage", pd)
    order = ["F0", "F1", "F2", "F3", "F4"]
    df["_order"] = df["fibrosis_stage"].astype(str).map({v: i for i, v in enumerate(order)}).fillna(999)
    df = df.sort_values(["_order", "fibrosis_stage"]).drop(columns="_order")
    out = args.out_dir / spec.filename
    _save_bar_with_labels(
        df,
        x_col="fibrosis_stage",
        y_col="unique_sample_count",
        title="Unique samples by fibrosis stage",
        xlabel="Fibrosis stage",
        ylabel="Unique sample_id count",
        out=out,
        color="#4C78A8",
    )
    return out


def plot_unique_samples_by_study(args: argparse.Namespace, spec: FigureSpec) -> Path:
    _np, pd, _plt = _import_plotting()
    df = _load_metadata_counts(args.metadata_csv, "study_accession", pd)
    df = df.sort_values("unique_sample_count", ascending=False)
    out = args.out_dir / spec.filename
    _save_bar_with_labels(
        df,
        x_col="study_accession",
        y_col="unique_sample_count",
        title="Unique samples by study accession",
        xlabel="Study accession",
        ylabel="Unique sample_id count",
        out=out,
        color="#F58518",
    )
    return out


def plot_cell_type_study_composition(args: argparse.Namespace, spec: FigureSpec) -> Path:
    _np, pd, plt = _import_plotting()
    obs, _umap = _load_obs_umap(args.h5ad)
    label_col = _choose_column(
        obs.columns,
        ["scanvi_pred", "tier_1_cell_type", "scanvi_predicted_cell_type", "cell_type", "annot"],
        "cell-type label",
    )
    study_col = _choose_column(obs.columns, ["study_accession", "study", "dataset"], "study accession")
    work = obs[[label_col, study_col]].dropna().copy()
    work[label_col] = work[label_col].astype(str)
    work[study_col] = work[study_col].astype(str)
    tab = pd.crosstab(work[label_col], work[study_col], normalize="index").sort_index()

    fig, ax = plt.subplots(figsize=(8.5, max(4, 0.32 * len(tab))), dpi=150)
    tab.plot(kind="barh", stacked=True, ax=ax, width=0.85)
    ax.set_xlabel("Proportion within cell type")
    ax.set_ylabel(label_col)
    ax.set_xlim(0, 1)
    ax.legend(title="study_accession", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
    fig.tight_layout()
    out = args.out_dir / spec.filename
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


PLOTTERS: dict[int, Callable[[argparse.Namespace, FigureSpec], Path]] = {
    1: plot_umap_study_accession,
    2: plot_unique_samples_by_fibrosis,
    3: plot_unique_samples_by_study,
    4: plot_cell_type_study_composition,
}


def _write_manifest(args: argparse.Namespace, outputs: list[Path]) -> Path:
    manifest = {
        "module": "a1_p1_qc.figures",
        "out_dir": str(args.out_dir),
        "h5ad": str(args.h5ad),
        "metadata_csv": str(args.metadata_csv),
        "metadata_filter": "include column truthy: true/1/yes/y/t",
        "figures": [
            {
                **asdict(spec),
                "output": str(args.out_dir / spec.filename),
            }
            for spec in FIGURES
        ],
        "generated": [str(p) for p in outputs],
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=False) + "\n")
    return args.manifest


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--h5ad",
        type=Path,
        default=Path("out/integration_annotation/all_studies/all_annotated.h5ad"),
    )
    parser.add_argument("--metadata-csv", type=Path, default=Path("config/liver_scRNAseq_metadata_everything.csv"))
    parser.add_argument("--out-dir", type=Path, default=Path("out/figures"))
    parser.add_argument("--manifest", type=Path, default=Path("out/figures/manifest.json"))
    parser.add_argument("--only", type=int, choices=[spec.order for spec in FIGURES], action="append")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    selected = set(args.only or [spec.order for spec in FIGURES])
    args.out_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for spec in FIGURES:
        if spec.order not in selected:
            continue
        outputs.append(PLOTTERS[spec.order](args, spec))
    _write_manifest(args, outputs)
    print(f"Generated {len(outputs)} a1_p1_qc figure(s) into {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
