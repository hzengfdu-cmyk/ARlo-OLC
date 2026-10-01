# Seurat input helpers

source("code/00_settings.R")
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
rna_layer <- function(object, layer) {
  if (packageVersion("SeuratObject") >= "5.0.0") {
    SeuratObject::LayerData(object,assay="RNA",layer=layer)
  } else {
    Seurat::GetAssayData(object,assay="RNA",slot=layer)
  }
}
read_rna <- function(path) {
  object <- if (grepl("\\.qs$", path, ignore.case=TRUE)) qs::qread(path) else readRDS(path)
  DefaultAssay(object) <- "RNA"
  if (inherits(object[["RNA"]], "Assay5")) {
    object <- JoinLayers(object, assay="RNA")
  }
  object
}
