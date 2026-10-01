# Cytokine-array and single-cell integration

source("code/00_setup.R")
out <- file.path(RESULTS,"11_cytokines"); dir.create(out,showWarnings=FALSE)
array <- read.csv(file.path(INPUT,"cytokine_array_gene_level.csv"),check.names=FALSE)
needed <- c("human_gene","category","pair","D_OCM","OCM")
array$human_gene <- trimws(as.character(array$human_gene))
array$D_OCM <- suppressWarnings(as.numeric(as.character(array$D_OCM)))
array$OCM <- suppressWarnings(as.numeric(as.character(array$OCM)))
array$valid_ratio <- is.finite(array$D_OCM)&is.finite(array$OCM)&array$D_OCM>0&array$OCM>0
array$log2_D_OCM_over_OCM <- NA_real_
good <- array$valid_ratio
array$log2_D_OCM_over_OCM[good] <- log2(array$D_OCM[good]/array$OCM[good])
summary_rows <- list()
for(g in unique(array$human_gene)) {
  d <- array[array$human_gene==g,,drop=FALSE]
  v <- d$valid_ratio
  summary_rows[[g]] <- data.frame(human_gene=g,category=d$category[1],n_pairs=nrow(d),
    n_valid_pairs=sum(v),mean_D_OCM=mean(d$D_OCM,na.rm=TRUE),
    mean_OCM=mean(d$OCM,na.rm=TRUE),
    mean_log2_D_OCM_over_OCM=if(any(v)) mean(d$log2_D_OCM_over_OCM[v]) else NA_real_)
}
arr_summary <- do.call(rbind,summary_rows)
arr_summary$log2_ratio_of_means_D_OCM_over_OCM <- NA_real_
ok <- is.finite(arr_summary$mean_D_OCM)&is.finite(arr_summary$mean_OCM)&
  arr_summary$mean_D_OCM>0&arr_summary$mean_OCM>0
arr_summary$log2_ratio_of_means_D_OCM_over_OCM[ok] <-
  log2(arr_summary$mean_D_OCM[ok]/arr_summary$mean_OCM[ok])
sce <- read_rna(file.path(RESULTS,"05_AUCell/human_scored.rds"))
oc <- subset(sce,cells=colnames(sce)[!is.na(sce$AR_state)])
genes <- intersect(unique(arr_summary$human_gene),rownames(oc))
de <- FindMarkers(oc,ident.1="ARhi",ident.2="ARlo",group.by="AR_state",assay="RNA",slot="data",
                  features=genes,test.use="wilcox",min.pct=0,logfc.threshold=0,only.pos=FALSE)
de$human_gene <- rownames(de)
x <- rna_layer(oc,"data")[genes,,drop=FALSE]
for(state in c("ARhi","ARlo")) {
  e <- x[,oc$AR_state==state,drop=FALSE]
  de[[paste0("mean_log_normalized_",state)]] <- Matrix::rowMeans(e)[rownames(de)]
  de[[paste0("pct_expressing_",state)]] <- Matrix::rowMeans(e>0)[rownames(de)]
}
joined <- merge(arr_summary,de,by="human_gene",all.x=TRUE)
write.csv(joined,file.path(out,"array_vs_OLC_ARhi_ARlo_source.csv"),row.names=FALSE)
