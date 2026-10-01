# AUCell scoring and ARhi/ARlo OLC states

source("code/00_setup.R")
library(AUCell)
out <- file.path(RESULTS,"05_AUCell"); dir.create(out,showWarnings=FALSE)
sce <- read_rna(file.path(RESULTS,"03_annotations/human_annotated.rds"))
meta <- sce[[]]
genes_long <- read.csv(file.path(CONFIG,"genesets_human_long.csv"))
geneSets <- split(genes_long$gene,genes_long$pathway)
geneSets <- lapply(geneSets,unique)
counts <- rna_layer(sce,"counts")
aucMaxRank <- ceiling(AUCMAX_FRACTION*nrow(counts))
cells <- rownames(meta)[meta$celltype %in% c("Osteoclast","Mph/Mono","CD4T","CD8T","Treg")]
scores <- matrix(NA_real_,length(geneSets),length(cells),dimnames=list(names(geneSets),cells))
blocks <- split(seq_along(cells),ceiling(seq_along(cells)/2000))
set.seed(SEED)
for(idx in blocks) {
  ranks <- AUCell_buildRankings(counts[,cells[idx],drop=FALSE],plotStats=FALSE,
                               BPPARAM=BiocParallel::SerialParam(),verbose=FALSE)
  a <- getAUC(AUCell_calcAUC(geneSets,ranks,aucMaxRank=aucMaxRank,verbose=FALSE))
  scores[rownames(a),colnames(a)] <- a
}
for(g in rownames(scores)) {
  v <- setNames(rep(NA_real_,ncol(sce)),colnames(sce)); v[colnames(scores)] <- scores[g,]
  sce[[paste0("auc_",g)]] <- v
}
eligible <- meta$celltype=="Osteoclast" & meta$specimen_type=="solid_lesion"
ar_cells <- rownames(meta)[eligible]
auc <- scores["Androgen_response",ar_cells]
ar_cutoff <- median(auc)
AR_state <- ifelse(auc>=ar_cutoff,"ARhi","ARlo")
v <- setNames(rep(NA_character_,ncol(sce)),colnames(sce)); v[ar_cells] <- AR_state
sce$AR_state <- v[colnames(sce)]
write.csv(data.frame(cell_id=ar_cells,auc=auc,AR_state=AR_state,cutoff=ar_cutoff),
          file.path(out,"OLC_AR_states.csv"),row.names=FALSE)
saveRDS(scores,file.path(out,"AUCell_matrix.rds"))
saveRDS(sce,file.path(out,"human_scored.rds"))
