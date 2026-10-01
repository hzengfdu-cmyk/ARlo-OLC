# Export counts for scTour

source("code/00_setup.R")
out <- file.path(RESULTS,"08_scTour"); dir.create(out,showWarnings=FALSE)
sce <- read_rna(file.path(RESULTS,"03_annotations/human_annotated.rds"))
cells <- colnames(sce)[sce$celltype %in% c("Mph/Mono","Osteoclast") & sce$specimen_type=="solid_lesion"]
counts <- rna_layer(sce,"counts")[,cells,drop=FALSE]
Matrix::writeMM(counts,file.path(out,"counts.mtx"))
write.csv(data.frame(gene=rownames(counts)),file.path(out,"genes.csv"),row.names=FALSE)
m <- sce[[]][cells,,drop=FALSE]; m$cell_id <- rownames(m)
write.csv(m,file.path(out,"metadata.csv"),row.names=FALSE)
if("tsne" %in% Reductions(sce)) write.csv(Embeddings(sce,"tsne")[cells,,drop=FALSE],file.path(out,"tsne.csv"))
