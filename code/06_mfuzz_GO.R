# Mfuzz expression modules and GO analysis

source("code/00_settings.R")
library(Matrix); library(ClusterGVis); library(clusterProfiler); library(org.Hs.eg.db)
out <- file.path(RESULTS,"06_Mfuzz"); dir.create(out,showWarnings=FALSE)
bundle <- readRDS(file.path(RESULTS,"03_annotations/human_RNA_portable.rds"))
m <- bundle$meta
keep <- m$celltype=="Osteoclast" & (m$specimen_type=="solid_lesion" | m$specimen_type %in% MFUZZ_MARROW_TYPES)
m <- m[keep,,drop=FALSE]
m$mfuzz_group <- ifelse(m$specimen_type %in% MFUZZ_MARROW_TYPES,"Marrow_comparator",m$stage)
groups <- c("Marrow_comparator","Treatment_naive","ADT_HSPC","CRPC")
counts <- bundle$counts[,rownames(m),drop=FALSE]
detected <- Matrix::rowMeans(counts>0)>0.1
x <- bundle$data[,rownames(m),drop=FALSE]
averages <- sapply(groups,function(g) Matrix::rowMeans(expm1(x[detected,m$mfuzz_group==g,drop=FALSE])))
averages <- averages[apply(averages,1,sd)>0,,drop=FALSE]
write.csv(averages,file.path(out,"group_average_linear_normalized_expression.csv"))
set.seed(MFUZZ_SEED)
if("obj" %in% names(formals(ClusterGVis::clusterData))) {
  cm <- ClusterGVis::clusterData(obj=as.data.frame(averages),cluster.method="mfuzz",cluster.num=5,seed=MFUZZ_SEED)
} else {
  cm <- ClusterGVis::clusterData(exp=as.data.frame(averages),cluster.method="mfuzz",cluster.num=5,seed=MFUZZ_SEED)
}
saveRDS(cm,file.path(out,"mfuzz_5_modules.rds"))
write.csv(cm$wide.res,file.path(out,"module_assignment.csv"),row.names=FALSE)
map <- bitr(rownames(averages),fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db)
universe <- unique(map$ENTREZID)
terms <- list()
for(k in sort(unique(cm$wide.res$cluster))) {
  genes <- cm$wide.res$gene[cm$wide.res$cluster==k]
  ids <- unique(map$ENTREZID[map$SYMBOL %in% genes])
  if(!length(ids)) next
  go <- enrichGO(ids,OrgDb=org.Hs.eg.db,keyType="ENTREZID",universe=universe,
                 ont="BP",pAdjustMethod="BH",pvalueCutoff=.05,qvalueCutoff=.2,readable=TRUE)
  d <- as.data.frame(go)
  if(!nrow(d)) next
  d$module <- k; terms[[as.character(k)]] <- d
}
if(length(terms)) {
  all_terms <- do.call(rbind,terms)
  top6 <- do.call(rbind,lapply(terms,function(d) head(d[order(d$pvalue,d$ID),],6)))
  write.csv(all_terms,file.path(out,"GO_BP_all_significant.csv"),row.names=FALSE)
  write.csv(top6,file.path(out,"GO_BP_top6_per_module.csv"),row.names=FALSE)
} else write.csv(data.frame(status="no_significant_term"),file.path(out,"GO_status.csv"),row.names=FALSE)
