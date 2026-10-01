# Human cell annotations and sample metadata

source("code/00_setup.R")
out <- file.path(RESULTS,"03_annotations"); dir.create(out,showWarnings=FALSE)
sce <- read_rna(HUMAN_ANNOTATED)
meta <- sce[[]]
manifest <- read.csv(file.path(INPUT,"human_sample_manifest.csv"),check.names=FALSE)
meta$sample_id <- trimws(as.character(meta[[SAMPLE_COLUMN]]))
aliases_file <- file.path(INPUT,"sample_aliases.csv")
if(file.exists(aliases_file)) {
  aliases <- read.csv(aliases_file)
  idx <- match(meta$sample_id,aliases$object_sample_id)
  meta$sample_id[!is.na(idx)] <- aliases$sample_id[idx[!is.na(idx)]]
}
missing <- setdiff(unique(meta$sample_id),manifest$sample_id)
idx <- match(meta$sample_id,manifest$sample_id)
for(nm in c("donor_id","dataset","specimen_type","stage","treatment")) meta[[nm]] <- manifest[[nm]][idx]
meta$celltype <- trimws(as.character(meta[[CELLTYPE_COLUMN]]))
meta$celltype[meta$celltype=="Mph_Mono"] <- "Mph/Mono"
sce <- AddMetaData(sce,meta)
sce <- NormalizeData(sce,assay="RNA",normalization.method="LogNormalize",scale.factor=10000)
write.csv(sce[[]],file.path(out,"all_cells_metadata.csv"))
write.csv(as.data.frame(table(sce$sample_id,sce$celltype)),file.path(out,"sample_celltype_counts.csv"),row.names=FALSE)
saveRDS(sce,file.path(out,"human_annotated.rds"))
saveRDS(list(counts=rna_layer(sce,"counts"),
             data=rna_layer(sce,"data"),meta=sce[[]]),
        file.path(out,"human_RNA_portable.rds"))
