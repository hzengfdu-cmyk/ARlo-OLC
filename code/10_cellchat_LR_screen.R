# State-associated ligand-receptor filtering

source("code/00_setup.R")
out <- file.path(RESULTS,"10_LR_screen"); dir.create(out,showWarnings=FALSE)
sce <- read_rna(file.path(RESULTS,"05_AUCell/human_scored.rds"))
oc <- subset(sce,cells=colnames(sce)[!is.na(sce$AR_state)])
Idents(oc) <- "AR_state"
ligand_de <- FindAllMarkers(oc,assay="RNA",slot="data",test.use="wilcox",only.pos=FALSE,
                            min.pct=0,logfc.threshold=0,return.thresh=Inf,random.seed=SEED)
ligand_hi <- ligand_de[as.character(ligand_de$cluster)=="ARhi",,drop=FALSE]
keep <- sce$specimen_type=="solid_lesion"
tme <- subset(sce,cells=colnames(sce)[keep])
Idents(tme) <- "treatment"
receptor_de <- FindAllMarkers(tme,assay="RNA",slot="data",test.use="wilcox",only.pos=FALSE,
                              min.pct=0,logfc.threshold=0,return.thresh=Inf,random.seed=SEED)
receptor_pre <- receptor_de[as.character(receptor_de$cluster)=="Pre",,drop=FALSE]
write.csv(ligand_de,file.path(out,"OLC_ARhi_ARlo_DE.csv"),row.names=FALSE)
write.csv(receptor_de,file.path(out,"TME_Pre_ADT_DE.csv"),row.names=FALSE)
db <- readRDS(file.path(RESULTS,"09_CellChat/CellChat_database.rds"))
net <- read.csv(file.path(RESULTS,"09_CellChat/Pooled_LR_network.csv"),check.names=FALSE)
net <- net[net$source %in% c("OLC_ARhi","OLC_ARlo"),,drop=FALSE]
net$pass_Methods <- FALSE
net$complex_rule <- "all_subunits_same_direction_and_adjusted_P_lt_0.05"
for(i in seq_len(nrow(net))) {
  components <- list(ligand=net$ligand[i],receptor=net$receptor[i])
  for(side in names(components)) {
    symbol <- components[[side]]
    if(symbol %in% rownames(db$complex)) {
      symbols <- as.character(unlist(db$complex[symbol,,drop=FALSE]))
      components[[side]] <- symbols[!is.na(symbols)&nzchar(symbols)]
    }
  }
  ld <- ligand_hi[match(components$ligand,ligand_hi$gene),,drop=FALSE]
  rd <- receptor_pre[match(components$receptor,receptor_pre$gene),,drop=FALSE]
  L <- ld$avg_log2FC; Lp <- ld$p_val_adj
  R <- rd$avg_log2FC; Rp <- rd$p_val_adj
  hi <- length(L)>0 && length(R)>0 && all(is.finite(c(L,Lp,R,Rp))) &&
    all(L>0 & Lp<.05) && all(R>0 & Rp<.05)
  lo <- length(L)>0 && length(R)>0 && all(is.finite(c(L,Lp,R,Rp))) &&
    all(L<0 & Lp<.05) && all(R<0 & Rp<.05)
  net$pass_Methods[i] <- if(net$source[i]=="OLC_ARhi") hi else lo
}
write.csv(net,file.path(out,"LR_screen_all.csv"),row.names=FALSE)
write.csv(net[net$pass_Methods,,drop=FALSE],file.path(out,"LR_screen_passed.csv"),row.names=FALSE)
