# Import bulk expression and clinical data

source("code/00_settings.R")
out <- file.path(RESULTS,"12_bulk"); dir.create(out,showWarnings=FALSE)
specs <- read.csv(file.path(INPUT,"bulk_expression_manifest.csv"),check.names=FALSE)
clinical <- read.csv(file.path(INPUT,"bulk_clinical.csv"),check.names=FALSE,stringsAsFactors=FALSE)
clinical$selected <- tolower(as.character(clinical$selected))=="true"
clinical$OS_time <- suppressWarnings(as.numeric(as.character(clinical$OS_time)))
clinical$OS_event <- suppressWarnings(as.numeric(as.character(clinical$OS_event)))
clinical$eligible <- clinical$selected & tolower(trimws(clinical$site))=="bone" &
  is.finite(clinical$OS_time)&clinical$OS_time>0&clinical$OS_event %in% c(0,1)
clinical <- clinical[clinical$eligible,,drop=FALSE]
clinical$OS_months <- ifelse(clinical$time_unit=="days",clinical$OS_time/(365.25/12),clinical$OS_time)
for(dataset in c("Merged","SU2C","WCDT")) {
  p <- specs$expression_file[match(dataset,specs$dataset)]
  x <- as.matrix(read.csv(p,row.names=1,check.names=FALSE))
  storage.mode(x) <- "double"
  ids <- if(dataset=="Merged") clinical$sample_id else clinical$sample_id[clinical$cohort==dataset]
  x <- x[,ids,drop=FALSE]
  saveRDS(x,file.path(out,paste0(dataset,"_supplied_log2TPM.rds")))
}
write.csv(clinical,file.path(out,"clinical_eligible.csv"),row.names=FALSE)
