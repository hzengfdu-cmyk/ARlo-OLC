# T-cell and myeloid subclustering

source("code/00_setup.R")
library(harmony)
out <- file.path(RESULTS,"03b_subclusters"); dir.create(out,showWarnings=FALSE)
sce <- read_rna(file.path(RESULTS,"03_annotations/human_annotated.rds"))
populations <- list(Tcell=c("Tcell","T cell","CD4T","CD8T","Treg"),
                    Myeloid=c("Myeloid","Mph/Mono","Osteoclast","DC","pDC","Neutrophil","Mast"))
for(pop in names(populations)) {
  cells <- colnames(sce)[sce$celltype %in% populations[[pop]]]
  x <- subset(sce,cells=cells)
  x <- NormalizeData(x,scale.factor=10000)
  x <- FindVariableFeatures(x,selection.method="vst",nfeatures=2000)
  x <- ScaleData(x,features=VariableFeatures(x))
  x <- RunPCA(x,npcs=30,seed.use=SEED)
  x <- harmony::RunHarmony(x,group.by.vars="sample_id",dims.use=1:30)
  x <- FindNeighbors(x,reduction="harmony",dims=1:20,graph.name=c("RNA_nn","RNA_snn"))
  x <- FindClusters(x,graph.name="RNA_snn",resolution=1,random.seed=SEED)
  x <- RunTSNE(x,reduction="harmony",dims=1:20,seed.use=SEED)
  markers <- FindAllMarkers(x,assay="RNA",slot="data",test.use="wilcox",only.pos=TRUE,
                            min.pct=.5,logfc.threshold=1)
  write.csv(markers,file.path(out,paste0(pop,"_annotation_markers.csv")),row.names=FALSE)
  write.csv(x[[]],file.path(out,paste0(pop,"_metadata.csv")))
  saveRDS(x,file.path(out,paste0(pop,"_subclusters.rds")))
}
