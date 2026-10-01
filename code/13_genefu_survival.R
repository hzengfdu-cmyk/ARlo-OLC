# ARlo signature and survival analysis

source("code/00_settings.R")
if(dir.exists(GENEFU_LIB)) .libPaths(c(GENEFU_LIB,.libPaths()))
library(genefu); library(survival)
out <- file.path(RESULTS,"13_survival"); dir.create(out,showWarnings=FALSE)
markers <- read.csv(file.path(RESULTS,"14_DE/ARlo_vs_other_OLC_and_MphMono.csv"),check.names=FALSE)
signature <- markers[is.finite(markers$avg_log2FC)&markers$avg_log2FC>=1,,drop=FALSE]
clinical_all <- read.csv(file.path(RESULTS,"12_bulk/clinical_eligible.csv"),check.names=FALSE)
expression_list <- list(
  Merged=readRDS(file.path(RESULTS,"12_bulk/Merged_supplied_log2TPM.rds")),
  SU2C=readRDS(file.path(RESULTS,"12_bulk/SU2C_supplied_log2TPM.rds")),
  WCDT=readRDS(file.path(RESULTS,"12_bulk/WCDT_supplied_log2TPM.rds")))
genes_used <- Reduce(intersect,c(list(signature$gene),lapply(expression_list,rownames)))
sig <- signature[match(genes_used,signature$gene),,drop=FALSE]
signature$used_in_all <- signature$gene %in% genes_used
models <- list()
for(dataset in names(expression_list)) {
  x <- expression_list[[dataset]]
  clinical <- if(dataset=="Merged") clinical_all else clinical_all[clinical_all$cohort==dataset,,drop=FALSE]
  x <- x[genes_used,clinical$sample_id,drop=FALSE]
  signature_matrix <- data.frame(probe=sig$gene,EntrezGene.ID=NA_character_,coefficient=sig$avg_log2FC)
  annot <- data.frame(EntrezGene.ID=rep(NA_character_,nrow(x)),row.names=rownames(x))
  score <- genefu::sig.score(x=signature_matrix,data=t(x),annot=annot,
                            do.mapping=FALSE,signed=FALSE,verbose=FALSE)$score
  clinical$score <- as.numeric(score[clinical$sample_id])
  valid <- is.finite(clinical$score)&is.finite(clinical$OS_months)&clinical$OS_months>0&clinical$OS_event %in% c(0,1)
  clinical <- clinical[valid,,drop=FALSE]
  if(dataset=="Merged") {
    cutoff <- median(clinical$score)
    clinical$group <- factor(ifelse(clinical$score>=cutoff,"High","Low"),levels=c("Low","High"))
    km <- survfit(Surv(OS_months,OS_event)~group,data=clinical)
    lr <- survdiff(Surv(OS_months,OS_event)~group,data=clinical)
    p <- pchisq(lr$chisq,df=1,lower.tail=FALSE)
    saveRDS(km,file.path(out,"pooled_KM.rds"))
    write.csv(data.frame(cutoff=cutoff,n=nrow(clinical),events=sum(clinical$OS_event),logrank_p=p),
              file.path(out,"pooled_logrank.csv"),row.names=FALSE)
    write.csv(data.frame(time=km$time,n_risk=km$n.risk,n_event=km$n.event,
      survival=km$surv,lower=km$lower,upper=km$upper,stratum=rep(names(km$strata),km$strata)),
      file.path(out,"KM_curve_source.csv"),row.names=FALSE)
  } else {
    clinical$score_z <- as.numeric(scale(clinical$score))
    fit <- coxph(Surv(OS_months,OS_event)~score_z,data=clinical,ties="efron",x=TRUE)
    s <- summary(fit); ph <- cox.zph(fit)
    models[[dataset]] <- data.frame(cohort=dataset,n=s$n,events=s$nevent,
      HR=exp(coef(fit))[1],lower=s$conf.int[1,"lower .95"],upper=s$conf.int[1,"upper .95"],
      p=s$coefficients[1,"Pr(>|z|)"],PH_p=ph$table[1,"p"],genes_used=length(genes_used))
    saveRDS(list(cox=fit,PH=ph),file.path(out,paste0(dataset,"_Cox.rds")))
  }
  write.csv(clinical,file.path(out,paste0(dataset,"_source_data.csv")),row.names=FALSE)
}
write.csv(do.call(rbind,models),file.path(out,"cohort_specific_Cox.csv"),row.names=FALSE)
