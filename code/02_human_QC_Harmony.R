# Human single-cell QC and Harmony integration

source("code/00_setup.R")
library(harmony)
out <- file.path(RESULTS, "02_human_QC"); dir.create(out, showWarnings=FALSE)
sce <- read_rna(HUMAN_RAW)
counts <- rna_layer(sce,"counts")
mt <- grep("^MT-", rownames(counts), value=TRUE)
meta <- sce[[]][colnames(counts), , drop=FALSE]
meta$nCount_raw <- Matrix::colSums(counts)
meta$nFeature_raw <- Matrix::colSums(counts > 0)
meta$percent_mt_raw <- 100 * Matrix::colSums(counts[mt,,drop=FALSE]) / meta$nCount_raw
calls_file <- if (file.exists(HUMAN_SCRUBLET_FILE)) {
  HUMAN_SCRUBLET_FILE
} else {
  file.path(RESULTS, "01_scrublet/doublet_calls.csv")
}
calls <- read.csv(calls_file, check.names=FALSE)
if (!"cell_id" %in% colnames(calls) && colnames(calls)[1] %in% c("", "X")) {
  colnames(calls)[1] <- "cell_id"
}
flag <- tolower(as.character(calls$predicted_doublet[match(rownames(meta),calls$cell_id)]))
meta$predicted_doublet <- flag == "true"
meta$QC_keep <- meta$nFeature_raw >= 250 & meta$nCount_raw >= 500 &
  is.finite(meta$percent_mt_raw) & meta$percent_mt_raw <= 20 & !meta$predicted_doublet
counts <- counts[,meta$QC_keep,drop=FALSE]
remove_gene <- trimws(readLines(HUMAN_GENE_REMOVE_FILE, warn=FALSE))
pseudogene <- trimws(readLines(HUMAN_PSEUDOGENE_FILE, warn=FALSE))
drop_gene <- unique(c(remove_gene, pseudogene))
counts <- counts[!rownames(counts) %in% drop_gene,,drop=FALSE]
symbol <- limma::alias2SymbolTable(rownames(counts), species="Hs")
original <- rownames(counts)
symbol[is.na(symbol) | symbol==""] <- original[is.na(symbol) | symbol==""]
target <- unique(symbol)
if (anyDuplicated(symbol)) {
  aggregate_map <- Matrix::sparseMatrix(i=match(symbol,target),
    j=seq_along(symbol),x=1,dims=c(length(target),length(symbol)))
  counts <- aggregate_map %*% counts
}
rownames(counts) <- target
counts <- counts[!rownames(counts) %in% drop_gene,,drop=FALSE]
keep_gene <- Matrix::rowSums(counts > 0) >= 50
counts <- counts[keep_gene,,drop=FALSE]
sce <- CreateSeuratObject(counts=counts, meta.data=meta[colnames(counts),,drop=FALSE])
sce <- NormalizeData(sce, normalization.method="LogNormalize", scale.factor=10000)
sce <- FindVariableFeatures(sce, selection.method="vst", nfeatures=2000)
sce <- ScaleData(sce, features=VariableFeatures(sce))
sce <- RunPCA(sce, features=VariableFeatures(sce), npcs=30, seed.use=SEED)
sce <- harmony::RunHarmony(sce, group.by.vars=SAMPLE_COLUMN, dims.use=1:30)
sce <- FindNeighbors(sce, reduction="harmony", dims=HUMAN_GRAPH_DIMS,
                     graph.name=c("RNA_nn","RNA_snn"))
sce <- FindClusters(sce, graph.name="RNA_snn", resolution=HUMAN_CLUSTER_RES, random.seed=SEED)
sce <- RunTSNE(sce, reduction="harmony", dims=HUMAN_GRAPH_DIMS, seed.use=SEED)
saveRDS(sce, file.path(out,"human_QC_Harmony.rds"))
write.csv(sce[[]], file.path(out,"cell_metadata_for_annotation.csv"))
