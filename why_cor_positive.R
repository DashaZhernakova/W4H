library(ggplotify)
library(RColorBrewer)
d <- read.delim("../../../../UMCG/data/Olink/Data/rawProteinData/CVD3_olinkNormal_1447_LLDsamples_t_ProtNames.txt", as.is = T, check.names = F, sep = "\t")
d <- read.delim("../../../../UMCG/data/Olink/Data/correctedForAgeGenderSmokingContracCellCounts/CVD3_olinkNormal_1447_LLDsamples_t_ENSBL.ProbesCentered.SamplesZTransformed.CovariatesRemoved.protnames.txt.gz", as.is = T, check.names = F, sep = "\t", row.names = 1)

d <- as.data.frame(t(d))
corm <- cor(d)

d2 <- d_wide[d_wide$TP == '1',]
d2$TP <- NULL
row.names(d2) <- d2$ID
d2$ID <- NULL

prots <- intersect(colnames(d2), colnames(d))

corm_2 <- cor(d2[,prots])
corm <- cor(d[,prots])

pheatmap(corm, fontsize = 6, cluster_rows = F, cluster_cols = F)
pheatmap(corm_2, fontsize = 6, cluster_rows = F, cluster_cols = F)

corm_3 <- cor(d2, use = 'complete.obs')
pheatmap(corm_3, fontsize = 6)


tmp <-my_pivot_wider(prot_lmm, "prot1", "prot2", "r")
tmp[is.na(tmp)] <- 0
pheatmap(cor(tmp[,prots]), fontsize = 6, cluster_rows = F, cluster_cols = F)


run_cor_tp <- function(tp, prots){
  d2 <- d_wide_adj_covar[d_wide_adj_covar$TP == tp,]
  d2$TP <- NULL
  row.names(d2) <- d2$ID
  d2$ID <- NULL
  
  corm_2 <- cor(d2[,prots], use = 'complete.obs')
  return(corm_2)
}

breaksList = seq(-1, 1, by = 0.2)
colorList <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))

p <- pheatmap(run_cor_tp('1', prots), fontsize = 6, treeheight_row = 0, treeheight_col =0, color = colorList, breaks = breaksList, main = 'TP1')
prot_ord <- prots[p$tree_row$order]

p1 <- as.ggplot(p)
p2 <- as.ggplot(pheatmap(run_cor_tp('2', prot_ord), fontsize = 6, cluster_rows = F, cluster_cols = F, color = colorList, breaks = breaksList, main = 'TP2'))
p3 <- as.ggplot(pheatmap(run_cor_tp('3', prot_ord), fontsize = 6, cluster_rows = F, cluster_cols = F, color = colorList, breaks = breaksList, main = 'TP3'))
p4 <- as.ggplot(pheatmap(run_cor_tp('4', prot_ord), fontsize = 6, cluster_rows = F, cluster_cols = F, color = colorList, breaks = breaksList, main = 'TP4'))

pdf("../plots/test_heatmap_corr_per_tp_all_prot.pdf", width = 10, height = 10)
(p1 + p2) / (p3 + p4)
dev.off()

pheatmap(corm[prot_ord,prot_ord], fontsize = 6, cluster_rows = F, cluster_cols = F, color = colorList, breaks = breaksList, main = 'LLD')



pheatmap(cor(d), fontsize = 6, cluster_rows = F, cluster_cols = F, color = colorList, breaks = breaksList, main = 'LLD all')



prots = c('IL6', 'IL18',  'IL11', 'IL4', 'IL10', 'TGFB1')
pheatmap(run_cor_tp('1', prots), fontsize = 6, cluster_rows = F, cluster_cols = F, color = colorList, breaks = breaksList, main = 'LLD IL')




d <- read.delim("../../../../UMCG/data/Olink/Data/rawProteinData/CVD3_olinkNormal_1447_LLDsamples_t_ProtNames.txt", as.is = T, check.names = F, sep = "\t", row.names = 1)
c <- read.delim("../../../../UMCG/data/Olink/Data/correctedForAgeGenderSmokingContracCellCounts/age_gender_smk_contrac_cell_counts.txt", as.is = T, check.names = F, sep = "\t", row.names = 1)

d <- as.data.frame(t(d))
c <- as.data.frame(t(c))

d_new <- regress_covariates_lld(d, c)
row.names(d_new) <- d_new$ID
d_new$ID <- NULL
pheatmap(cor(d_new), fontsize = 6, cluster_rows = F, cluster_cols = F, color = colorList, breaks = breaksList, main = 'LLD adj cov', treeheight_row = 0, treeheight_col =0)
pheatmap(cor(d), fontsize = 6, cluster_rows = F, cluster_cols = F, color = colorList, breaks = breaksList, main = 'LLD adj cov', treeheight_row = 0, treeheight_col =0)


regress_covariates_lld <- function(data_tp, covariates){
  
  ids <- intersect(row.names(covariates), row.names(data_tp))
  data_tp <- data_tp[ids,]
  covariates <- covariates[ids,]
  
  new_d <- data.frame(matrix(ncol = ncol(data_tp), nrow = nrow(data_tp)))
  colnames(new_d) <- colnames(data_tp)
  row.names(new_d) <- row.names(data_tp)
  for (col in 1:ncol(data_tp)){
    tmp <- as.data.frame(cbind(data_tp[,col], covariates))
    colnames(tmp)[1] <- "prot"
    form <- as.formula(paste("prot ~", paste(colnames(covariates), collapse = "+")))
    lm_fit <- lm(form, data = tmp)
    new_d[,col] <- residuals(lm_fit)
  }
  new_d <- new_d %>%
    rownames_to_column(var = 'ID')
  
  return(new_d)
}


d2 <- d_wide[d_wide$TP == '1',]
d2$TP <- NULL
row.names(d2) <- d2$ID
d2$ID <- NULL

p <- pheatmap(cor(d2, use = 'complete.obs'), fontsize = 6, treeheight_row = 0, treeheight_col =0, color = colorList, breaks = breaksList)
prot_ord <- prots[p$tree_row$order]
pheatmap(run_cor_tp('1',prot_ord), fontsize = 6, treeheight_row = 0, treeheight_col =0, color = colorList, breaks = breaksList, main = 'TP1')



####
for (i in 1: nrow(prot_lmm)){
  p1 <- prot_lmm[i, "prot1"]
  p2 <- prot_lmm[i, "prot2"]
  
  for (tp in 1:4){
    cor(d_wide[d_wide$TP == tp, p1], d_wide[d_wide$TP == tp, p2], use = 'complete.obs')
  }

}



