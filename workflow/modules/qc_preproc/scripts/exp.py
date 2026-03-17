from __future__ import annotations
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("WebAgg")

sc.settings.set_figure_params(dpi=50, facecolor="white")

exp_path = "out/GSE256398/raw/GSM8648658_R5_CB_raw_feature_bc_matrix_filtered.h5"

sample_adata = sc.read_10x_h5(exp_path)
sample_adata.var_names_make_unique()
adata = sample_adata

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

adata2 = sc.read_10x_mtx("")
