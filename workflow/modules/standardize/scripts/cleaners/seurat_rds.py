from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path
from typing import TextIO

import anndata as ad


def _require_path(value: str | None, name: str) -> Path:
    if not value:
        raise SystemExit(f"Missing required input '{name}'")
    path = Path(value)
    if not path.exists():
        raise SystemExit(f"Input '{name}' not found: {path}")
    return path


def run_seurat_rds(*, study: str, study_spec: dict, output_h5ad: Path, log: TextIO) -> None:
    inp = study_spec.get("input", {})
    rds_path = _require_path(inp.get("rds"), "rds")

    log.write(f"[{study}] Converting Seurat RDS to h5ad: {rds_path}\n")

    with tempfile.TemporaryDirectory(prefix=f"standardize_{study}_") as tmp:
        tmpdir = Path(tmp)
        out_h5ad = tmpdir / "converted.h5ad"
        stdout_json = tmpdir / "summary.json"

        r_script = r"""
args <- commandArgs(trailingOnly=TRUE)
rds_path <- args[[1]]
out_h5ad <- args[[2]]
out_json <- args[[3]]

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(zellkonverter))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(jsonlite))

obj <- readRDS(rds_path)

# Best-effort: accept Seurat object or SingleCellExperiment.
sce <- NULL
if (inherits(obj, 'Seurat')) {
  sce <- as.SingleCellExperiment(obj)
} else if (inherits(obj, 'SingleCellExperiment')) {
  sce <- obj
} else {
  stop(paste('Unsupported R object class:', class(obj)[1]))
}

zellkonverter::writeH5AD(sce, out_h5ad)

summary <- list(
  n_cells = ncol(sce),
  n_genes = nrow(sce),
  assays = names(SummarizedExperiment::assays(sce)),
  coldata_cols = colnames(SummarizedExperiment::colData(sce))
)
write(toJSON(summary, auto_unbox=TRUE, pretty=TRUE), file=out_json)
"""

        cmd = [
            "Rscript",
            "-e",
            r_script,
            str(rds_path),
            str(out_h5ad),
            str(stdout_json),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            log.write(proc.stdout)
            log.write(proc.stderr)
            raise SystemExit(f"R conversion failed for {study}: {rds_path}")

        if stdout_json.exists():
            log.write(f"[{study}] R conversion summary: {json.loads(stdout_json.read_text())}\n")

        log.write(f"[{study}] Reading converted h5ad\n")
        adata = ad.read_h5ad(out_h5ad)
        adata.obs["study_accession"] = study
        adata.uns["standardize"] = {"study": study, "process": "seurat_rds"}

        log.write(f"[{study}] Writing: {output_h5ad}\n")
        adata.write_h5ad(output_h5ad)

