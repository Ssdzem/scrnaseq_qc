import marimo

__generated_with = "0.23.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import os
    import sys
    from pathlib import Path

    print("cwd:", os.getcwd())                 # where script is running from
    print("__file__:", __file__)              # script file path
    print("script dir:", Path(__file__).parent.resolve())
    print("argv[0]:", sys.argv[0])            # how script was invoked
    return


@app.cell
def _():
    import anndata as ad
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import plotly.express as px
    import marimo as mo

    return ad, mo, pd, px


@app.cell
def _(ad):
    adata = ad.read_h5ad("lSMIF/a1_p1_qc/neighbors_umap_only.h5ad")
    return (adata,)


@app.cell
def _(adata):
    adata
    return


@app.cell
def _(adata):
    adata.obs
    return


@app.cell
def _(adata):
    cols = [
        "scanvi_pred",
        "scanvi_pred_confidence",
        "scanvi_pred_entropy",
        "scanvi_pred_margin",
        "sample_id",
        "study_accession",
        "donor_id",
    ]

    obs = adata.obs[cols].copy()
    obs["cell_id"] = obs.index.astype(str)

    obs.to_parquet("scanvi_qc_obs.parquet", index=False)
    return


@app.cell
def _(pd, px):
    df = pd.read_parquet("scanvi_qc_obs.parquet")

    fig = px.violin(
        df,
        x="scanvi_pred",
        y="scanvi_pred_confidence",
        box=True,
        points=False,
    )

    fig
    return (df,)


@app.cell
def _(mo):
    metric = mo.ui.dropdown(
        options=[
            "scanvi_pred_confidence",
            "scanvi_pred_entropy",
            "scanvi_pred_margin",
        ],
        value="scanvi_pred_confidence",
        label="Metric",
    )

    metric
    return (metric,)


@app.cell
def _(df, metric, px):
    order = (
        df.groupby("scanvi_pred")[metric.value]
        .median()
        .sort_values()
        .index
    )

    fig2 = px.violin(
        df,
        x="scanvi_pred",
        y=metric.value,
        category_orders={"scanvi_pred": list(order)},
        box=True,
        points=False,
        height=800,
    )

    fig2.update_layout(xaxis_tickangle=45)
    fig2
    return


if __name__ == "__main__":
    app.run()
