# Mouse single-cell QC and Harmony integration

source("code/00_setup.R")
library(harmony)
out <- file.path(RESULTS,"04_mouse"); dir.create(out,showWarnings=FALSE)
sce <- read_rna(MOUSE_RAW)
counts <- rna_layer(sce,"counts")
mt <- grep("^mt-",rownames(counts),value=TRUE)
meta <- sce[[]][colnames(counts),,drop=FALSE]
meta$nFeature_raw <- Matrix::colSums(counts>0)
meta$percent_mt_raw <- 100*Matrix::colSums(counts[mt,,drop=FALSE])/Matrix::colSums(counts)
meta$QC_keep <- meta$nFeature_raw>200 & is.finite(meta$percent_mt_raw) & meta$percent_mt_raw<=25
sce <- subset(sce,cells=rownames(meta)[meta$QC_keep])
sce <- NormalizeData(sce,scale.factor=10000)
sce <- FindVariableFeatures(sce,selection.method="vst",nfeatures=2000)
sce <- ScaleData(sce,features=VariableFeatures(sce))
sce <- RunPCA(sce,npcs=30,seed.use=SEED)
sce <- harmony::RunHarmony(sce,group.by.vars=SAMPLE_COLUMN,dims.use=1:30)
sce <- FindNeighbors(sce,reduction="harmony",dims=MOUSE_GRAPH_DIMS,graph.name=c("RNA_nn","RNA_snn"))
sce <- FindClusters(sce,graph.name="RNA_snn",resolution=1,random.seed=SEED)
sce <- RunTSNE(sce,reduction="harmony",dims=MOUSE_GRAPH_DIMS,seed.use=SEED)
saveRDS(sce,file.path(out,"mouse_QC_Harmony.rds"))
write.csv(sce[[]],file.path(out,"mouse_metadata_for_annotation.csv"))
