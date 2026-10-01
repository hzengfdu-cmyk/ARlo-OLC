# Differential expression and enrichment analysis

source("code/00_setup.R")
library(clusterProfiler); library(org.Hs.eg.db)
out <- file.path(RESULTS,"14_DE"); dir.create(out,showWarnings=FALSE)
sce <- read_rna(file.path(RESULTS,"05_AUCell/human_scored.rds"))
oc <- subset(sce,cells=colnames(sce)[!is.na(sce$AR_state)])
state_de <- FindMarkers(oc,ident.1="ARhi",ident.2="ARlo",group.by="AR_state",test.use="wilcox",
                        assay="RNA",slot="data",min.pct=0,logfc.threshold=0,only.pos=FALSE)
state_de$gene <- rownames(state_de)
write.csv(state_de,file.path(out,"ARhi_vs_ARlo_OLC.csv"),row.names=FALSE)
cells <- colnames(sce)[sce$specimen_type=="solid_lesion" & sce$celltype %in% c("Osteoclast","Mph/Mono")]
myeloid <- subset(sce,cells=cells)
myeloid$signature_group <- ifelse(!is.na(myeloid$AR_state)&myeloid$AR_state=="ARlo","ARlo","Other_OLC_MphMono")
signature_de <- FindMarkers(myeloid,ident.1="ARlo",ident.2="Other_OLC_MphMono",group.by="signature_group",
                            assay="RNA",slot="data",test.use="wilcox",min.pct=0,logfc.threshold=0,only.pos=FALSE)
signature_de$gene <- rownames(signature_de)
write.csv(signature_de,file.path(out,"ARlo_vs_other_OLC_and_MphMono.csv"),row.names=FALSE)
map <- bitr(state_de$gene,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db)
universe <- unique(map$ENTREZID)
for(direction in c("ARhi","ARlo")) {
  yes <- is.finite(state_de$p_val_adj)&state_de$p_val_adj<DE_ENRICH_PADJ &
    if(direction=="ARhi") state_de$avg_log2FC>DE_ENRICH_LOGFC else state_de$avg_log2FC< -DE_ENRICH_LOGFC
  ids <- unique(map$ENTREZID[map$SYMBOL %in% state_de$gene[yes]])
  if(!length(ids)) {
    write.csv(data.frame(status="no_selected_genes"),file.path(out,paste0(direction,"_enrichment_status.csv")),row.names=FALSE)
    next
  }
  go <- enrichGO(ids,universe=universe,OrgDb=org.Hs.eg.db,ont="BP",keyType="ENTREZID",
                 pAdjustMethod="BH",pvalueCutoff=.05,qvalueCutoff=.2,readable=TRUE)
  write.csv(as.data.frame(go),file.path(out,paste0(direction,"_GO_BP.csv")),row.names=FALSE)
  kegg <- tryCatch(enrichKEGG(ids,universe=universe,organism="hsa",pAdjustMethod="BH",pvalueCutoff=.05),error=identity)
  if(inherits(kegg,"error")) {
    writeLines(conditionMessage(kegg),file.path(out,paste0(direction,"_KEGG_error.txt")))
  } else write.csv(as.data.frame(kegg),file.path(out,paste0(direction,"_KEGG.csv")),row.names=FALSE)
}
