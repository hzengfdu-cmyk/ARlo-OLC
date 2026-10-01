# Cell composition and Ro/e analysis

source("code/00_setup.R")
out <- file.path(RESULTS,"15_composition"); dir.create(out,showWarnings=FALSE)
sce <- read_rna(file.path(RESULTS,"03_annotations/human_annotated.rds"))
m <- sce[[]]
m <- m[m$specimen_type=="solid_lesion" | m$specimen_type %in% MFUZZ_MARROW_TYPES,,drop=FALSE]
m$analysis_group <- ifelse(m$specimen_type=="solid_lesion",m$stage,"Marrow_comparator")
counts <- table(m$celltype,m$analysis_group)
expected <- outer(rowSums(counts),colSums(counts))/sum(counts)
roe <- counts/expected
label <- matrix(as.character(cut(as.vector(roe),breaks=c(-Inf,.1,.5,1.5,2,Inf),
  labels=c("-","+/-","+","++","+++"),right=FALSE)),nrow=nrow(roe),dimnames=dimnames(roe))
write.csv(counts,file.path(out,"pooled_counts.csv"))
write.csv(roe,file.path(out,"pooled_RoE.csv"))
write.csv(label,file.path(out,"pooled_RoE_labels.csv"))
sample_counts <- table(m$sample_id,m$celltype)
write.csv(sample_counts,file.path(out,"per_sample_counts.csv"))
write.csv(prop.table(sample_counts,1),file.path(out,"per_sample_proportions.csv"))
