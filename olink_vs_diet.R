
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/")
source("scripts/utility_functions.R")


library(ggplot2)
library(rmcorr)
library(dplyr)
library(lme4)
library(grid)
library(gridExtra)
library(pheatmap)
library(corrplot)
library(patchwork)

out_basedir <- "results/pheno_batch2_prot_rm_outliers_4sd/"
d_wide <- read.delim("data/olink_clean_CVD+INF_rm_outliers_4sd.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))
d_wide$TP <- as.numeric(d_wide$TP)

covariates <- read.delim("data/covariates_age_bmi_storage_preg.txt", sep = "\t", check.names = F, as.is = T, colClasses = c(ID = "character"))
covariates[] <- lapply(covariates, function(col) {
  if (length(unique(col)) < 3) {
    return(factor(col))
  } else {
    return(col)
  }
})

covariate_names <- c("Age","BMI","storage_months", "Pregnancy_category")
d_wide_adj_covar = read.delim(paste0(out_basedir, "olink_clean_CVD+INF_adj_covariates.txt"), sep = "\t", as.is = T)

pca <- read.delim(paste0(out_basedir, "olink_clean_CVD+INF_rm_outliers_4sd.PCA.txt"), as.is = T, check.names = F, sep = "\t",  colClasses = c(ID = "character"))

ffq <- read.delim("../phenotypes/DietFreeze_DataCleaned_Feb2025/fmt_dasha/FFQ.txt", as.is = T, sep = "\t", colClasses = c(ID = "character"))
diet_weekly <- read.delim("../phenotypes/DietFreeze_DataCleaned_Feb2025/fmt_dasha/FFQ_DIARI_ms_long.txt", as.is = T, sep = "\t", colClasses = c(ID = "character"))
diet_3gg <- read.delim("../phenotypes/DietFreeze_DataCleaned_Feb2025/fmt_dasha/FFQ_DIARI_3gg_long.txt", as.is = T, sep = "\t", colClasses = c(ID = "character"))

################################################################################
# FFQ Diet for TP = 1
################################################################################
pca_1 <- pca[pca$TP == 1,1:7]
pca_1$TP <- NULL

d_wide_1 <- d_wide[d_wide$TP == 1,]
d_wide_1$TP <- NULL
d_wide_1$SampleID <- NULL

ffq_pca_lm <- data.frame()
for (diet_item in colnames(ffq)[2:ncol(ffq)]){
  for (pc in colnames(pca_1)[2:ncol(pca_1)]){
    d_subs <- full_join(ffq[,c("ID", diet_item)], pca_1[,c("ID", pc)], by = "ID")
    colnames(d_subs) <- c("ID", "ffq", "pc")
    fit <- lm(pc ~ ffq, data = d_subs)
    ffq_pca_lm <- rbind(ffq_pca_lm, c(diet_item, pc, summary(fit)$coefficients['ffq', c(1,2,4)]))
  }
}
colnames(ffq_pca_lm) <- c("FFQ_item", "PC", "beta", "SE", "pval")
ffq_pca_lm <- na.omit(ffq_pca_lm) %>%
  mutate(across(-c( FFQ_item, PC), as.numeric)) 

write.table(ffq_pca_lm, file = paste0(out_basedir, "correlations_with_covariates/PCA_vs_FFQ_tp1_lm.txt"), quote = F, sep = "\t", row.names = FALSE)


ffq_prot_lm <- data.frame()
for (diet_item in colnames(ffq)[2:ncol(ffq)]){
  for (prot in colnames(d_wide_1)[2:ncol(d_wide_1)]){
    d_subs <- full_join(ffq[,c("ID", diet_item)], d_wide_1[,c("ID", prot)], by = "ID")
    colnames(d_subs) <- c("ID", "ffq", "prot")
    fit <- lm(prot ~ ffq, data = d_subs)
    ffq_prot_lm <- rbind(ffq_prot_lm, c(diet_item, prot, summary(fit)$coefficients['ffq', c(1,2,4)]))
  }
}
colnames(ffq_prot_lm) <- c("FFQ_item", "prot", "beta", "SE", "pval")
ffq_prot_lm <- na.omit(ffq_prot_lm) %>%
  mutate(across(-c( FFQ_item, prot), as.numeric)) 

ffq_prot_lm$BH_pval <- p.adjust(ffq_prot_lm$pval, method = 'BH')
ffq_prot_lm$bonf_sign <- ifelse(ffq_prot_lm$pval < 0.05/(39*length(unique(ffq_prot_lm$FFQ_item))), T, F)
write.table(ffq_prot_lm, file = paste0(out_basedir, "correlations_with_covariates/prot_vs_FFQ_tp1_lm.txt"), quote = F, sep = "\t", row.names = FALSE)



################################################################################
# Longitudinal data
################################################################################
diet <- diet_weekly
#diet <- diet_3gg
diet <- cbind(data.frame(SampleID = paste0(diet$ID, "_", diet$TP)), diet)


  gam_res_diet <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(diet) -4), ncol = 10))
  colnames(gam_res_diet) <- c("prot", "diet_item", "gam_pval", "gam_est", "gam_se", "n", " n_samples",
                              "lmm_est", "lmm_pval", "lmm_se")
  
  cnt <- 1
  for (diet_item in colnames(diet)[4:ncol(diet)]) {
    cat(diet_item, "\n")
    for (prot in colnames(d_wide)[4:ncol(d_wide)]){
      res_gam <- gam_prot_pheno_adj_covar(d_wide, diet, prot, diet_item, covariates, scale = T, adjust_timepoint = 'spline', anova_pval = F)
      res_lmm <- lmm_pheno_prot_adj_covar(d_wide, diet, prot, diet_item, covariates, scale = T, adjust_timepoint = 'cubic')
      
      gam_res_diet[cnt,] <- c(prot, diet_item, unlist(res_gam), unlist(res_lmm)[1:3])
      cnt <- cnt + 1
    }
  }
  
  gam_res_diet <- na.omit(gam_res_diet) %>%
    mutate(across(-c(diet_item, prot), as.numeric)) 
  gam_res_diet$BH_gam_pval <- p.adjust(gam_res_diet$gam_pval, method = 'BH')
  gam_res_diet$BH_lmm_pval <- p.adjust(gam_res_diet$lmm_pval, method = 'BH')
  bonf_cutoff <- 0.5 / (length(unique(gam_res_diet$diet_item)) * 39)
  gam_res_diet$BH_gam_bonferroni <- ifelse(gam_res_diet$gam_pval < bonf_cutoff, T, F)
  gam_res_diet$BH_lmm_bonferroni <- ifelse(gam_res_diet$lmm_pval < bonf_cutoff, T, F)
  
#write.table(gam_res_diet, file = paste0(out_basedir, "correlations_with_covariates/prot_vs_diet_3gg_gam_lmm_adj_covar.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(gam_res_diet, file = paste0(out_basedir, "correlations_with_covariates/prot_vs_diet_weekly_gam_lmm_adj_covar.txt"), quote = F, sep = "\t", row.names = FALSE)



#### Corrplots
gam_res_diet_w = read.delim(paste0(out_basedir, "correlations_with_covariates/prot_vs_diet_weekly_gam_lmm_adj_covar.txt"), as.is = T, check.names = F, sep = "\t")
gam_res_diet_3gg = read.delim(paste0(out_basedir, "correlations_with_covariates/prot_vs_diet_3gg_gam_lmm_adj_covar.txt"), as.is = T, check.names = F, sep = "\t")
gam_res_diet_3gg$diet_item <- gsub("3g", "3gg", gam_res_diet_3gg$diet_item)

prots <- unique(c(gam_res_diet_3gg[gam_res_diet_3gg$lmm_pval < 0.01, "prot"],gam_res_diet_w[gam_res_diet_w$lmm_pval < 0.01, "prot"]))

res_wide <- my_pivot_wider(gam_res_diet_3gg[gam_res_diet_3gg$prot %in% prots,], "diet_item", "prot", "lmm_est")
pval_wide <- my_pivot_wider(gam_res_diet_3gg[gam_res_diet_3gg$prot %in% prots,], "diet_item", "prot", "lmm_pval")
sign_labels <- my_pivot_wider(gam_res_diet_3gg[gam_res_diet_3gg$prot %in% prots,], "diet_item", "prot", "BH_lmm_bonferroni")
sign_labels[sign_labels == TRUE] <- "*"
sign_labels[sign_labels == FALSE] <- ""

max_val <- max(abs(min(res_wide)), max(res_wide))
breaksList = seq(-max_val, max_val, by = 0.01)
colorList <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
h <- pheatmap(res_wide, display_numbers = sign_labels, fontsize_number = 14, 
         color = colorList, breaks = breaksList, fontsize_col = 6,
         cellwidth=8, cellheight=12, 
         filename = paste0(out_basedir, "correlations_with_covariates/plots/prot_vs_diet_3gg_lmm_adj_covar_pval0.01.pdf"))

row_order <- row.names(res_wide)[h$tree_row$order]
col_order <- colnames(res_wide)[h$tree_col$order]
row_order <- gsub("_3gg", "", row_order)


res_wide_w <- my_pivot_wider(gam_res_diet_w[gam_res_diet_w$prot %in% prots,], "diet_item", "prot", "lmm_est")
pval_wide_w <- my_pivot_wider(gam_res_diet_w[gam_res_diet_w$prot %in% prots,], "diet_item", "prot", "lmm_pval")
sign_labels_w <- my_pivot_wider(gam_res_diet_w[gam_res_diet_w$prot %in% prots,], "diet_item", "prot", "BH_lmm_bonferroni")
sign_labels_w[sign_labels_w == TRUE] <- "*"
sign_labels_w[sign_labels_w == FALSE] <- ""

pheatmap(res_wide_w[row_order, col_order], display_numbers = sign_labels_w[row_order, col_order], fontsize_number = 14, 
         color = colorList, breaks = breaksList, fontsize_col = 6,
         cellwidth=8, cellheight=12, cluster_rows = F, cluster_cols = F,
         filename = paste0(out_basedir, "correlations_with_covariates/plots/prot_vs_diet_weekly_lmm_adj_covar_pval0.01.pdf"))

sign_combined <- rbind(
  gam_res_diet[gam_res_diet$BH_lmm_bonferroni == T, c("prot", "diet_item", "lmm_est", "lmm_pval")],
  gam_res_diet_w[gam_res_diet_w$BH_lmm_bonferroni == T, c("prot", "diet_item", "lmm_est", "lmm_pval")])


attribs <- rbind(
  data.frame(node = sign_combined[,"prot"], type = 'protein'),
  data.frame(node = sign_combined[,"diet_item"], type = 'diet')
)
attribs <- unique(attribs)
attribs_colors <- c("protein" = "#56B4E9", "diet" = "#E69F00")

combined_edges <- sign_combined[,c("prot", "diet_item")]
net <- graph_from_data_frame(d=combined_edges, directed=F) 

E(net)$weight <- sign_combined$lmm_est
V(net)$type <- attribs$type[match(V(net)$name, attribs$node)]
V(net)$color <- attribs_colors[V(net)$type]

pdf(paste0(out_basedir, "correlations_with_covariates/plots/diet_network_plot_r.pdf"), width = 10, height = 10)
plot(net, 
     edge.width = abs(E(net)$weight) * 5, # Scale edge width by weight
     edge.color = ifelse(E(net)$weight > 0, "#D55E00", "#009E73"), 
     vertex.size = 16, edge.width = 2,
     vertex.label.cex = 0.8,
     vertex.label.color = "black")
dev.off()



