"""Plot raw (pre-trim) QC metric distributions for a standardized study.

Computes standard per-cell QC metrics from out/{study}/raw.h5ad (which, for
studies with no `input.qc` block in their standardize spec, has not been
filtered at all) and writes a per-cell metrics table plus histograms split by
sample_id, so cutoffs can be picked by eye instead of preset in the spec.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


METRICS = [
    ("log1p_total_counts", "log1p(total_counts)"),
    ("log1p_n_genes_by_counts", "log1p(n_genes_by_counts)"),
    ("pct_counts_mt", "% mitochondrial counts"),
]


def compute_qc_metrics(adata, mito_gene_prefix: str):
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(mito_gene_prefix.upper())
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=True)
    cols = ["sample_id", "total_counts", "n_genes_by_counts", "pct_counts_mt",
            "log1p_total_counts", "log1p_n_genes_by_counts"]
    return adata.obs[cols].copy()


def plot_metric_by_sample(df, metric: str, xlabel: str, out_path: Path) -> None:
    samples = sorted(df["sample_id"].astype(str).unique())
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for sample_id in samples:
        values = df.loc[df["sample_id"] == sample_id, metric].to_numpy()
        ax.hist(values, bins=60, histtype="step", linewidth=1.5, label=f"{sample_id} (n={len(values)})")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("cells")
    ax.set_title(f"{xlabel} — raw, pre-trim")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_counts_vs_genes(df, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    scatter = ax.scatter(
        df["total_counts"], df["n_genes_by_counts"], c=df["pct_counts_mt"],
        s=4, cmap="viridis", alpha=0.6, linewidths=0,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("total_counts (log scale)")
    ax.set_ylabel("n_genes_by_counts (log scale)")
    ax.set_title("counts vs genes, colored by %mito — raw, pre-trim")
    fig.colorbar(scatter, ax=ax, label="pct_counts_mt")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True)
    parser.add_argument("--input", required=True, help="Path to out/{study}/raw.h5ad")
    parser.add_argument("--mito-gene-prefix", default="MT-")
    parser.add_argument("--output-metrics", required=True)
    parser.add_argument("--output-counts-fig", required=True)
    parser.add_argument("--output-genes-fig", required=True)
    parser.add_argument("--output-mito-fig", required=True)
    parser.add_argument("--output-scatter-fig", required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    if "sample_id" not in adata.obs.columns:
        raise SystemExit(f"[{args.study}] raw.h5ad is missing obs['sample_id']; cannot split QC by sample")

    df = compute_qc_metrics(adata, args.mito_gene_prefix)
    df.insert(0, "study_accession", args.study)

    metrics_path = Path(args.output_metrics)
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(metrics_path, index=False)
    print(f"[{args.study}] wrote per-cell QC metrics ({len(df)} cells): {metrics_path}")

    fig_paths = {
        "log1p_total_counts": Path(args.output_counts_fig),
        "log1p_n_genes_by_counts": Path(args.output_genes_fig),
        "pct_counts_mt": Path(args.output_mito_fig),
    }
    for metric, xlabel in METRICS:
        plot_metric_by_sample(df, metric, xlabel, fig_paths[metric])
        print(f"[{args.study}] wrote {metric} figure: {fig_paths[metric]}")

    plot_counts_vs_genes(df, Path(args.output_scatter_fig))
    print(f"[{args.study}] wrote counts-vs-genes scatter: {args.output_scatter_fig}")


if __name__ == "__main__":
    main()
