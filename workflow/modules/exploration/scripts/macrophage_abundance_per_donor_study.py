import marimo

__generated_with = "0.23.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import anndata as ad
    import scanpy as sc
    import numpy as np, scipy.sparse as sp
    import pandas as pd
    import plotly.express as px
    import marimo as mo

    return ad, pd


@app.cell
def _(ad):
    adata = ad.read_h5ad("lSMIF/a1_p1_qc/neighbors_umap_only.h5ad")
    return (adata,)


@app.cell
def _(adata, pd):
    tab = pd.crosstab(
        adata.obs["scanvi_pred"],
        adata.obs["study_accession"],
        normalize="index"
    )

    tab.plot(kind="barh", stacked=True, figsize=(10,6))
    return


@app.cell
def _(adata, pd):
    obs = adata.obs.copy()

    # global study proportions
    global_prop = obs["study_accession"].value_counts(normalize=True)

    # celltype x study proportions
    ct = pd.crosstab(
        obs["scanvi_pred"],
        obs["study_accession"],
        normalize="index"
    )

    # observed / expected
    oe = ct.div(global_prop, axis=1)

    oe
    return obs, oe


@app.cell
def _(oe):
    import seaborn as sns
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10,7))
    sns.heatmap(
        oe,
        cmap="coolwarm",
        center=1,
        annot=True,
        fmt=".2f"
    )
    plt.title("Observed / Expected study contribution")
    return plt, sns


@app.cell
def _(obs):
    obs["scanvi_pred"].value_counts().loc["Macrophages"]
    return


@app.cell
def _(adata):
    adata.obs["scanvi_pred"].value_counts().loc["Mono+mono derived cells"]
    return


@app.cell
def _(adata):
    adata.obs["scanvi_pred"].value_counts().loc["Mig.cDCs"]
    return


@app.cell
def _(adata):
    adata.obs["scanvi_pred"].value_counts().loc["cDC1s"]
    return


@app.cell
def _(adata):
    adata.obs["scanvi_pred"].value_counts().loc["cDC2s"]
    return


@app.cell
def _(adata):
    adata.obs["scanvi_pred"].value_counts().loc["pDCs"]
    return


@app.cell
def _(obs, pd):
    total = obs.groupby("donor_id").size()
    mac = obs[obs["scanvi_pred"]=="Macrophages"].groupby("donor_id").size()

    res = pd.DataFrame({
        "total_cells": total,
        "macrophages": mac
    }).fillna(0)

    res["mac_fraction"] = res["macrophages"] / res["total_cells"]

    study = obs.groupby("donor_id")["study_accession"].first()
    res["study_accession"] = study.values
    res["donor_id"] = res.index
    return (res,)


@app.cell
def _(plt, res, sns):
    plt.figure(figsize=(16,6))
    sns.barplot(
        data=res.sort_values("mac_fraction", ascending=False),
        x="donor_id",
        y="mac_fraction",
        hue="study_accession",
        dodge=False
    )
    plt.xticks(rotation=90)
    plt.title("Macrophage fraction per donor")
    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(res, sns):
    sns.histplot(res["macrophages"], bins=40)
    return


@app.cell
def _(plt, res, sns):
    plt.figure(figsize=(10,6))
    sns.scatterplot(
        data=res,
        x="total_cells",
        y="mac_fraction",
        hue="study_accession",
        s=80
    )
    plt.xscale("log")
    plt.show()
    return


@app.cell
def _():
    from scipy.stats import kruskal
    from statsmodels.stats.multicomp import pairwise_tukeyhsd

    return kruskal, pairwise_tukeyhsd


@app.cell
def _(kruskal, obs, pairwise_tukeyhsd, pd):
    # --- donor-level macrophage fraction table ---
    total2 = obs.groupby("donor_id").size()
    mac2 = obs[obs["scanvi_pred"] == "Macrophages"].groupby("donor_id").size()

    res2 = pd.DataFrame({
        "total_cells": total2,
        "macrophages": mac2
    }).fillna(0)

    res2["mac_fraction"] = res2["macrophages"] / res2["total_cells"]
    res2["study_accession"] = obs.groupby("donor_id")["study_accession"].first()

    # --- Kruskal-Wallis across studies ---
    groups = [
        g["mac_fraction"].values
        for _, g in res2.groupby("study_accession")
    ]

    stat, pval = kruskal(*groups)

    print(f"Kruskal-Wallis H={stat:.4f}, p={pval:.3e}")

    # --- descriptive stats ---
    print(
        res2.groupby("study_accession")["mac_fraction"]
        .agg(["count","median","mean","std","min","max"])
        .sort_values("median", ascending=False)
    )

    # --- optional pairwise posthoc (parametric-ish, exploratory) ---
    posthoc = pairwise_tukeyhsd(
        endog=res2["mac_fraction"],
        groups=res2["study_accession"],
        alpha=0.05
    )

    print(posthoc)
    return


@app.cell
def _(adata, pd, res):
    tmp = res.copy()
    tmp["fibrosis"] = adata.obs.groupby("donor_id")["fibrosis_stage"].first()

    pd.crosstab(
        tmp["study_accession"],
        tmp["fibrosis"],
        values=tmp["macrophages"],
        aggfunc="count"
    )
    return (tmp,)


@app.cell
def _(tmp):
    tmp.groupby(["study_accession", "fibrosis"])["macrophages"].agg(
        ["count", "median", "min", "max"]
    )
    return


@app.cell
def _(res):
    res["usable_macrophage_donor"] = res["macrophages"] >= 50
    return


@app.cell
def _(adata, res):
    tmp3 = res.copy()
    tmp3["fibrosis"] = adata.obs.groupby("donor_id")["fibrosis_stage"].first()
    tmp3["usable_macrophage_donor"] = res["usable_macrophage_donor"]
    return (tmp3,)


@app.cell
def _(pd, tmp3):
    pd.crosstab(
        tmp3["study_accession"],
        tmp3["fibrosis"],
        values=tmp3["usable_macrophage_donor"],
        aggfunc="sum"
    )
    return


@app.cell
def _(pd, plt, sns, tmp3):
    tab3 = pd.crosstab(
        tmp3["study_accession"],
        tmp3["fibrosis"],
        values=tmp3["usable_macrophage_donor"],
        aggfunc="sum"
    ).fillna(0)

    plt.figure(figsize=(10,5))
    sns.heatmap(tab3, annot=True, fmt=".0f", cmap="viridis")
    plt.title("Usable macrophage donors by Study × Fibrosis")
    plt.ylabel("Study")
    plt.xlabel("Fibrosis")
    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(adata):
    keep_labels = [
        "Macrophages",
        "Mono+mono derived cells",
        "cDC1s",
        "cDC2s",
        "Mig.cDCs",
        "pDCs",
        "Neutrophils",   # optional, useful as contaminant/boundary
    ]

    adata_myeloid = adata[adata.obs["scanvi_pred"].isin(keep_labels)].copy()
    return


if __name__ == "__main__":
    app.run()
