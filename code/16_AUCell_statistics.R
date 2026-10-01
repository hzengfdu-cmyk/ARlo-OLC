# AUCell group comparisons

source("code/00_setup.R")
out <- file.path(RESULTS,"16_AUCell_statistics"); dir.create(out,showWarnings=FALSE)
sce <- read_rna(file.path(RESULTS,"05_AUCell/human_scored.rds"))
m <- sce[[]]; columns <- grep("^auc_",names(m),value=TRUE)
tests <- list(); stats <- list(); k <- 0L
for(pop in c("OLC_AR","OLC_categories","MphMono_TNFR1","Tcell")) {
  if(pop=="OLC_AR") {
    d <- m[!is.na(m$AR_state),,drop=FALSE]; d$comparison_group <- d$AR_state
  } else if(pop=="OLC_categories") {
    d <- m[m$celltype=="Osteoclast" & (m$specimen_type=="solid_lesion" | m$specimen_type %in% MFUZZ_MARROW_TYPES),,drop=FALSE]
    d$comparison_group <- ifelse(d$specimen_type=="solid_lesion",d$stage,"Marrow_comparator")
  } else if(pop=="MphMono_TNFR1") {
    d <- m[m$celltype=="Mph/Mono" & m$specimen_type=="solid_lesion",,drop=FALSE]
    count <- rna_layer(sce,"counts")
    d$comparison_group <- ifelse(as.numeric(count["TNFRSF1A",rownames(d)])>0,"Positive","Negative")
  } else {
    d <- m[m$celltype %in% c("CD4T","CD8T","Treg") & m$specimen_type=="solid_lesion",,drop=FALSE]
    d$comparison_group <- d$treatment
  }
  for(col in columns) {
    z <- d[is.finite(d[[col]]) & !is.na(d$comparison_group),,drop=FALSE]
    if(!nrow(z)) next
    groups <- unique(z$comparison_group); k <- k+1L
    p <- NA_real_; method <- "not_testable"
    if(length(groups)==2 && all(table(z$comparison_group)>=3)) {
      p <- wilcox.test(z[[col]][z$comparison_group==groups[1]],z[[col]][z$comparison_group==groups[2]],exact=FALSE)$p.value
      method <- "two_sided_cell_level_Wilcoxon"
    } else if(length(groups)>2) {
      p <- kruskal.test(z[[col]],as.factor(z$comparison_group))$p.value; method <- "cell_level_Kruskal_Wallis"
    }
    tests[[k]] <- data.frame(population=pop,pathway=col,n_cells=nrow(z),n_donors=length(unique(z$donor_id)),method=method,p_value=p)
    for(g in groups) {
      v <- z[[col]][z$comparison_group==g]
      stats[[paste(pop,col,g)]] <- data.frame(population=pop,pathway=col,group=g,n=length(v),mean=mean(v),SD=sd(v))
    }
    donor <- aggregate(z[[col]],by=z[c("dataset","donor_id","comparison_group")],FUN=mean)
    names(donor)[ncol(donor)] <- "mean_AUC"
    write.csv(donor,file.path(out,paste0(pop,"_",col,"_donor_means.csv")),row.names=FALSE)
  }
}
d <- do.call(rbind,tests)
d$p_adj_BH_this_table <- p.adjust(d$p_value,"BH")
write.csv(d,file.path(out,"cell_level_tests.csv"),row.names=FALSE)
write.csv(do.call(rbind,stats),file.path(out,"mean_SD.csv"),row.names=FALSE)
