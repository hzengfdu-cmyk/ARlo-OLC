# CellChat communication analysis

source("code/00_settings.R")
library(Matrix); library(CellChat)
out <- file.path(RESULTS,"09_CellChat"); dir.create(out,showWarnings=FALSE)
bundle <- readRDS(file.path(RESULTS,"03_annotations/human_RNA_portable.rds"))
meta <- bundle$meta
ar <- read.csv(file.path(RESULTS,"05_AUCell/OLC_AR_states.csv"))
meta$AR_state <- ar$AR_state[match(rownames(meta),ar$cell_id)]
meta$cc_group <- meta$celltype
meta$cc_group[!is.na(meta$AR_state)] <- paste0("OLC_",meta$AR_state[!is.na(meta$AR_state)])
results <- list()
for(condition in c("Pre","ADT","Pooled")) {
  take <- meta$specimen_type=="solid_lesion"
  if(condition!="Pooled") take <- take & meta$treatment==condition
  m <- droplevels(meta[take,,drop=FALSE])
  n <- table(m$cc_group)
  m <- m[m$cc_group %in% names(n)[n>=CELLCHAT_MIN_CELLS],,drop=FALSE]
  x <- bundle$data[,rownames(m),drop=FALSE]
  cc <- createCellChat(object=x,meta=m,group.by="cc_group")
  cc@DB <- subsetDB(CellChatDB.human,search="Secreted Signaling")
  saveRDS(cc@DB,file.path(out,"CellChat_database.rds"))
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc,raw.use=TRUE,population.size=CELLCHAT_POPULATION_SIZE,seed.use=SEED)
  cc <- filterCommunication(cc,min.cells=CELLCHAT_MIN_CELLS)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc,slot.name="netP")
  saveRDS(cc,file.path(out,paste0(condition,"_cellchat.rds")))
  write.csv(subsetCommunication(cc),file.path(out,paste0(condition,"_LR_network.csv")),row.names=FALSE)
  write.csv(cc@net$count,file.path(out,paste0(condition,"_interaction_count.csv")))
  write.csv(cc@net$weight,file.path(out,paste0(condition,"_interaction_weight.csv")))
  for(path in names(cc@netP$centr)) {
    z <- cc@netP$centr[[path]]
    results[[paste(condition,path)]] <- data.frame(condition=condition,pathway=path,
      celltype=names(z$outdeg),outgoing=z$outdeg,incoming=z$indeg)
  }
}
cent <- do.call(rbind,results)
write.csv(cent,file.path(out,"centrality_long.csv"),row.names=FALSE)
pre <- cent[cent$condition=="Pre",c("celltype","pathway","outgoing","incoming")]
post <- cent[cent$condition=="ADT",c("celltype","pathway","outgoing","incoming")]
delta <- merge(pre,post,by=c("celltype","pathway"),all=TRUE,suffixes=c("_Pre","_ADT"))
delta$outgoing_Pre_minus_ADT <- delta$outgoing_Pre-delta$outgoing_ADT
delta$incoming_Pre_minus_ADT <- delta$incoming_Pre-delta$incoming_ADT
write.csv(delta,file.path(out,"centrality_Pre_minus_ADT.csv"),row.names=FALSE)
