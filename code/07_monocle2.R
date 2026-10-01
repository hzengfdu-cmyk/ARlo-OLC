# Monocle2 OLC trajectory

source("code/00_setup.R")
library(monocle)
out <- file.path(RESULTS,"07_Monocle2"); dir.create(out,showWarnings=FALSE)
sce <- read_rna(file.path(RESULTS,"05_AUCell/human_scored.rds"))
oc <- subset(sce,cells=colnames(sce)[!is.na(sce$AR_state)])
counts <- rna_layer(oc,"counts")
pd <- new("AnnotatedDataFrame",data=oc[[]][colnames(counts),,drop=FALSE])
fd <- new("AnnotatedDataFrame",data=data.frame(gene_short_name=rownames(counts),row.names=rownames(counts)))
cds <- newCellDataSet(counts,phenoData=pd,featureData=fd,expressionFamily=negbinomial.size())
cds <- detectGenes(cds,min_expr=1)
cds <- cds[fData(cds)$num_cells_expressed>10,]
cds <- estimateSizeFactors(cds); cds <- estimateDispersions(cds)
oc <- FindVariableFeatures(oc,selection.method="vst",nfeatures=2000)
ordering_genes <- intersect(VariableFeatures(oc),rownames(cds))
write.csv(data.frame(gene=ordering_genes),file.path(out,"ordering_genes.csv"),row.names=FALSE)
cds <- setOrderingFilter(cds,ordering_genes)
cds <- reduceDimension(cds,max_components=2,reduction_method="DDRTree")
cds <- orderCells(cds)
if(!is.na(MONOCLE_ROOT_STATE)) {
  cds <- orderCells(cds,root_state=MONOCLE_ROOT_STATE)
}
saveRDS(cds,file.path(out,"OLC_monocle2.rds"))
d <- pData(cds); d$cell_id <- rownames(d); d$root_confirmed <- !is.na(MONOCLE_ROOT_STATE)
write.csv(d,file.path(out,"cell_pseudotime.csv"),row.names=FALSE)
