suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratDisk)
})

option_list <- list(
  make_option(c("--input-rds"), type = "character", help = "Input Seurat .rds"),
  make_option(c("--output-h5ad"), type = "character", help = "Output .h5ad"),
  make_option(c("--log"), type = "character", default = NA, help = "Log file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$log)) {
  dir.create(dirname(opt$log), showWarnings = FALSE, recursive = TRUE)
  zz <- file(opt$log, open = "wt")
  sink(zz, type = "output")
  sink(zz, type = "message")
}

cat("Reading:", opt$`input-rds`, "\n")
obj <- readRDS(opt$`input-rds`)

# Some older Seurat objects (especially those with SCT assays) can fail
# conversion in SeuratDisk due to legacy SCTModel slots. Dropping SCT
# keeps the main RNA assay for downstream Scanpy QC.
if ("SCT" %in% names(obj@assays)) {
  cat("Dropping SCT assay to improve conversion compatibility\n")
  obj[["SCT"]] <- NULL
}

tmp_h5seurat <- sub("\\.h5ad$", ".h5seurat", opt$`output-h5ad`)
cat("Writing temporary:", tmp_h5seurat, "\n")
SaveH5Seurat(obj, filename = tmp_h5seurat, overwrite = TRUE)

cat("Converting to h5ad:", opt$`output-h5ad`, "\n")
Convert(tmp_h5seurat, dest = "h5ad", overwrite = TRUE)

cat("Done.\n")
