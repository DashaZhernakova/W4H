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

out_basedir <- "results/pheno_batch2_prot_rm_outliers_4sd/"
d_wide <- read.delim("data/olink_clean_CVD+INF_rm_outliers_4sd.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))
covariates <- read.delim("data/covariates_age_bmi_storage_preg.txt", sep = "\t", check.names = F, as.is = T, colClasses = c(ID = "character"))

d_wide$TP <- as.numeric(d_wide$TP)

covariates[] <- lapply(covariates, function(col) {
  if (length(unique(col)) < 3) {
    return(factor(col))
  } else {
    return(col)
  }
})


# Make a dataframe with proteins adjusted for all covariates per visit
covariate_names <- c("Age","BMI","storage_months", "Pregnancy_category")

joined_data <- full_join(covariates, d_wide, by = c("ID"), relationship = 'one-to-many')
d_wide_adj_covar <- regress_covariates_lmm(d_wide, covariates, covars_longitudinal = F)
joined_data_adj_covar <- full_join(covariates, d_wide_adj_covar, by = c("ID"), relationship = 'one-to-many')

write.table(d_wide_adj_covar, file = paste0(out_basedir, "olink_clean_CVD+INF_adj_covariates.txt"), quote = F, sep = "\t", row.names = FALSE)


################################################################################
# ICC for each protein
################################################################################

icc <- data.frame(matrix(nrow = (ncol(d_wide) -4), ncol = 2))
colnames(icc) <- c("prot", "ICC")
cnt <- 1
for (prot in colnames(d_wide)[4:ncol(d_wide)]){
  res <- get_ICC(d_wide, prot)
  icc[cnt,] <- c(prot, res)
  cnt <- cnt + 1
}

icc <- na.omit(icc) %>%
  mutate(across(-c( prot), as.numeric)) 
icc <- icc[order(icc$ICC, decreasing = F),]

cat("ICC ranges from", min(icc$ICC), "to", max(icc$ICC), "with a median of", median(icc$ICC), "\n")


write.table(icc, file = paste0(out_basedir, "ICC_per_protein.txt"), quote = F, sep = "\t", row.names = FALSE)

pdf(paste0(out_basedir, "plots/ICC_per_protein.pdf"))
p1 <- ggplot(icc, aes(x = ICC)) + 
  geom_density() + 
  theme_minimal() 

p2 <- ggplot(icc, aes(y = ICC)) + 
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.text.x=element_blank())
p1 + p2
dev.off()

################################################################################
# Protein vs TP GAM and LMM
################################################################################

gam_res_prot_tp <- data.frame(matrix(nrow = (ncol(d_wide) -4), ncol = 7))
colnames(gam_res_prot_tp) <- c("prot", "gam_pval", "gam_edf", "gam_fval", "n", "n_samples", "lmm_pval")


cnt <- 1
for (prot in colnames(d_wide)[4:ncol(d_wide)]){
  res_gam <- gam_prot_tp_adj_covar(d_wide, prot, covariates, scale = T, predict = F)
  res_lmm <- lmm_prot_tp_poly3_adj_covar(d_wide, prot, covariates)
  gam_res_prot_tp[cnt,] <- c(prot, unlist(res_gam), res_lmm)
  cnt <- cnt + 1
}

gam_res_prot_tp <- na.omit(gam_res_prot_tp) %>%
  mutate(across(-c( prot), as.numeric)) 

gam_res_prot_tp$gam_BH_pval <- p.adjust(gam_res_prot_tp$gam_pval, method = 'BH')
gam_res_prot_tp$lmm_BH_pval <- p.adjust(gam_res_prot_tp$lmm_pval, method = 'BH')

gam_res_prot_tp <- gam_res_prot_tp[order(gam_res_prot_tp$gam_pval),]
gam_res_prot_tp$gam_bonf_sign <- ifelse(gam_res_prot_tp$gam_pval < 0.05/38,T,F)
gam_res_prot_tp$gam_edf_round <- round(gam_res_prot_tp$gam_edf)
gam_res_prot_tp$gam_BH_sign <- ifelse(gam_res_prot_tp$gam_BH_pval < 0.05,T,F)

write.table(gam_res_prot_tp, file = paste0(out_basedir, "prot_vs_tp_gam_adj_age_bmi_preg_storage.txt"), quote = F, sep = "\t", row.names = FALSE)
signif <- gam_res_prot_tp[gam_res_prot_tp$gam_BH_pval < 0.05,]

cat ("Number of proteins that change significantly with time:", nrow(signif), "\n")
cat("Of them, the number of proteins with a non-linear change: ", nrow(signif[signif$gam_edf_round > 1,]), "\n")
cat("Of them, the number of proteins showing a significant association with time also in LMMs:", nrow(signif[signif$lmm_BH_pval < 0.05,]))

################################################################################
# Cluster longitudinal trajectories: GAM
################################################################################

n_points = 100

# Take only non-linear trajectories:
all_prots <- signif[signif$gam_edf_round > 1, "prot"]
# Take all significant trajectories
all_prots <- signif$prot
n_prots <- length(all_prots)
prot_trajs <- data.frame(matrix(nrow = length(all_prots) , ncol = n_points))
row.names(prot_trajs) <- all_prots
colnames(prot_trajs) <- seq(1,4, length.out = n_points)

prot_dist <- data.frame(matrix(nrow = length(all_prots) * length(all_prots), ncol = 3))
colnames(prot_dist) <- c("prot1", "prot2", "eucl")

pb = txtProgressBar(min = 0, max = length(all_prots), initial = 0) 
stepi = 0
cnt <- 1
prot1_passed <- c()
for (prot1 in all_prots){
  gam_fit1 <- gam_prot_tp_adj_covar(d_wide, prot1, covariates, scale = T, predict = T)
  prot_trajs[prot1,] <- gam_fit1$predicted
  
  setTxtProgressBar(pb,stepi)
  stepi <- stepi + 1
  
  prot1_passed <- c(prot1_passed, prot1)
  
  for (prot2 in all_prots_gam){
    if (prot2 %in% prot1_passed) next
    
    gam_fit2 <- gam_prot_tp_adj_covar(d_wide, prot2, covariates, scale = T, predict = T)
    
    # Euclidean distance between trajectories
    eucl_dist <- sum(abs(as.numeric(gam_fit1$predicted) - as.numeric(gam_fit2$predicted)))/length(gam_fit1$predicted)
    
    prot_dist[cnt,] <- c(prot1, prot2, eucl_dist)
    cnt <- cnt + 1
  }
  
}

cat("\n")
close(pb)

prot_dist <- na.omit(prot_dist) %>%
  mutate(across(-c( prot1, prot2), as.numeric)) 

# select only non-linear trajectories
all_prots <- signif[signif$gam_edf_round > 1, "prot"]
prot_dist_all_sign <- prot_dist
prot_dist <- prot_dist[prot_dist$prot1 %in% all_prots & prot_dist$prot2 %in% all_prots, ]
n_prots = length(all_prots)

prot_trajs_all_sign <- prot_trajs
prot_trajs <- prot_trajs[all_prots, ]

close_prots <- prot_dist[prot_dist$eucl < quantile(prot_dist$eucl, 0.05),]
close_prots <- close_prots[order(close_prots$eucl),]
#dir.create(file.path(out_basedir, "trajectories_gam"))
pdf(paste0(out_basedir, "trajectories_gam/gam_prot_prot_nonlinear_eucl_q0.05.pdf"), height = 20, width = 15)
plot_list <- list()
cnt <- 1
for (i in 1:nrow(close_prots)){
  plot_list[[cnt]] <- plot_together(prot = close_prots[i, "prot1"], ph = close_prots[i, "prot2"], trajectories = prot_trajs)
  cnt <- cnt + 1
  if(cnt == 21){
    grid.arrange(grobs = plot_list, ncol = 4, nrow = 5)
    plot_list <- list()
    cnt <- 1
  }
}
grid.arrange(grobs = plot_list, ncol = 4, nrow = 5)
dev.off()

far_prots <- prot_dist[prot_dist$eucl > quantile(prot_dist$eucl, 0.95),]
far_prots <- far_prots[order(far_prots$eucl, decreasing = T),]
pdf(paste0(out_basedir, "trajectories_gam/gam_prot_prot_nonlinear_eucl_q0.95.pdf"), height = 20, width = 15)
plot_list <- list()
cnt <- 1
for (i in 1:nrow(far_prots)){
  plot_list[[cnt]] <- plot_together(prot = far_prots[i, "prot1"], ph = far_prots[i, "prot2"], trajectories = prot_trajs)
  cnt <- cnt + 1
  if(cnt == 21){
    grid.arrange(grobs = plot_list, ncol = 4, nrow = 5)
    plot_list <- list()
    cnt <- 1
  }
}
grid.arrange(grobs = plot_list, ncol = 4, nrow = 5)
dev.off()

# fix the issue that PROK1 has always high distances
mean_dist <- prot_dist  %>% 
  group_by(prot1) %>%
  mutate(mean = mean(eucl)) %>% 
  distinct(prot1, mean)
mean_dif <- mean_dist[mean_dist$prot1 == 'PROK1', "mean"] - max(mean_dist[mean_dist$prot1 != 'PROK1', "mean"])

#convert distance to similarity metric
prot_dist$eucl_similarity <- 1 / (1 + prot_dist$eucl)

prot_dist_rescaled <- prot_dist
prot_dist_rescaled[prot_dist_rescaled$prot1 == 'PROK1', "eucl"] <- prot_dist_rescaled[prot_dist_rescaled$prot1 == 'PROK1', "eucl"] - as.numeric(mean_dif)

tmp <- data.frame("prot1" = prot_dist_rescaled$prot2, "prot2" = prot_dist_rescaled$prot1, "eucl" = prot_dist_rescaled$eucl, "eucl_similarity" = prot_dist_rescaled$eucl_similarity)
dist_matrix <- my_pivot_wider(rbind(prot_dist_rescaled, tmp), "prot1", "prot2", "eucl")
dist_matrix <- dist_matrix[colnames(dist_matrix),]

sim_matrix <- my_pivot_wider(rbind(prot_dist_rescaled, tmp), "prot1", "prot2", "eucl_similarity")
sim_matrix <- sim_matrix[colnames(sim_matrix),]

write.table(sim_matrix, file = paste0(out_basedir, "trajectories_gam/similarity_matrix_eucl_", n_prots, "_prots.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(dist_matrix, file = paste0(out_basedir, "trajectories_gam/distance_matrix_eucl_", n_prots, "_prots.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(prot_dist, file = paste0(out_basedir, "trajectories_gam/protein_distance_eucl_", n_prots, "_prots.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(prot_trajs, file = paste0(out_basedir, "trajectories_gam/protein_trajectories_", n_prots, "_prots.txt"), quote = F, sep = "\t", col.names = NA, row.names = T)

write.table(prot_dist_all_sign, file = paste0(out_basedir, "trajectories_gam/protein_distance_eucl_all_sign_prots.txt"), quote = F, sep = "\t", col.names = NA, row.names = T)
write.table(prot_trajs_all_sign, file = paste0(out_basedir, "trajectories_gam/protein_trajectories_all_sign_prots.txt"), quote = F, sep = "\t", col.names = NA, row.names = T)

# Heatmap and clustering on euclidean distance between curves
# method = 'hclust_my_eucl_dist'
# num_k = 5
# hm <- pheatmap(dist_matrix, cutree_rows = num_k, cutree_cols = num_k, filename = paste0(out_basedir,"trajectories_gam/", method, "_k", num_k, "_", n_prots, "_prots_heatmap.pdf"))
# dev.off()
# cl  = cutree(hm$tree_row, k = num_k)
# plot_clusters(cl, method, num_k, prot_trajs = prot_trajs, out_path = paste0(out_basedir,"trajectories_gam/", method, "_k", num_k, "_", n_prots, "_prots.pdf"))
# write.table(as.data.frame(cl), file = paste0(out_basedir, "trajectories_gam/", method, "_k", num_k, "_", n_prots, "_prots.txt"), quote = F, sep = "\t")

# Get pairwise differences between i and i+1 timepoint and cluster according to them

prot_deltas <- data.frame(matrix(nrow = nrow(prot_trajs), ncol = ncol(prot_trajs) - 1))

row.names(prot_deltas) <- row.names(prot_trajs)

for (c in 1:(ncol(prot_trajs) - 1)){
  prot_deltas[,c] <- prot_trajs[,(c + 1)] - prot_trajs[,c]
}
method = 'pam_delta'
num_k = 8
#num_k = 15
pam_res <- cluster::pam(prot_deltas, num_k, diss = F)
cl <- pam_res$clustering
plot_clusters(cl, method, num_k, prot_trajs = prot_trajs, add_cluster_name = T, out_path = paste0(out_basedir,"trajectories_gam/", method, "_k", num_k, "_", n_prots, "_prots.pdf"))
write.table(as.data.frame(cl), file = paste0(out_basedir, "trajectories_gam/", method, "_k", num_k, "_", n_prots, "_prots.txt"), quote = F, sep = "\t")
non_linear_clusters <- as.data.frame(cl)
non_linear_clusters$cl <- paste0("nonlinear_", non_linear_clusters$cl)

# Linear trajectories
linear_prots <- signif[signif$gam_edf_round == 1, "prot"]
method = 'pam_traj_linear'
num_k = 3
#num_k = 15
pam_res <- cluster::pam(prot_trajs_all_sign[linear_prots,], num_k, diss = F)
cl <- pam_res$clustering
plot_clusters(cl, method, num_k, prot_trajs = prot_trajs_all_sign, out_path = paste0(out_basedir,"trajectories_gam/", method, "_k", num_k, "_", length(linear_prots), "_prots.pdf"))

cl_2 <- cl
cl_2[cl_2 == 2] <- 1
cl_2[cl_2 == 3] <- 2
num_k = "2_2"
plot_clusters(cl_2, method, num_k, add_cluster_name = T, prot_trajs = prot_trajs_all_sign, out_path = paste0(out_basedir,"trajectories_gam/", method, "_k", num_k, "_", length(linear_prots), "_prots.pdf"))
write.table(as.data.frame(cl_2), file = paste0(out_basedir, "trajectories_gam/", method, "_k", num_k, "_", n_prots, "_prots.txt"), quote = F, sep = "\t")

linear_clusters <- as.data.frame(cl_2)
colnames(linear_clusters)[1] <- "cl"
linear_clusters$cl <- paste0("linear_", linear_clusters$cl)

all_clusters <- rbind(linear_clusters, non_linear_clusters)
colnames(all_clusters)[1] <- "trajectory_cluster"
all_clusters <- all_clusters %>% rownames_to_column(var = "prot")
gam_res_prot_tp <- full_join(gam_res_prot_tp, all_clusters, by = "prot")
write.table(gam_res_prot_tp, file = paste0(out_basedir, "prot_vs_tp_gam_adj_age_bmi_preg_storage.txt"), quote = F, sep = "\t", row.names = FALSE)


################################################################################
# Plot all significant protein trajectories 
################################################################################
prot_trajs_all_sign <- read.delim(paste0(out_basedir, "trajectories_gam/protein_trajectories_all_sign_prots.txt"), as.is = T, row.names = 1, sep = "\t")

pdf(paste0(out_basedir,"/plots/proteins_vs_tp_gam_withpoints0.pdf"), height = 20, width = 15)

plot_list <- list()
cnt <- 1 
for (prot in row.names(prot_trajs_all_sign)){
  traj <- data.frame(TP = as.numeric(colnames(prot_trajs_all_sign)), prot = as.numeric(prot_trajs_all_sign[prot,]))
  traj <- full_join(d_wide_adj_covar[,c("TP", "ID", prot)], traj, by = "TP")
  colnames(traj)[3] <- "values"
  traj$values <- scale(traj$values)
  pval <- gam_res_prot_tp[gam_res_prot_tp$prot == prot, "gam_pval"]
  plot_list[[prot]] <- ggplot(traj) + 
    geom_boxplot(aes(x = TP, y = values, group = TP), width = 0.1, color = my_colors[3], outliers = F) +
    geom_line(aes(x = TP, y = prot), color = my_colors[4]) + 
    geom_jitter(aes(x = TP, y = values), alpha = 0.2, width = 0.1, color = my_colors[3]) +
    labs(x = "Timepoint ", y = prot) +
    ggtitle(paste0(prot, "; GAM P = ", formatC(pval, digits = 3))) + 
    theme_minimal() 
  cnt <- cnt + 1
  
  if (cnt == 21){
    grid.arrange(grobs = plot_list, ncol = 4, nrow = 5)  
    plot_list <- list()
    cnt <- 1
  }
}

grid.arrange(grobs = plot_list, ncol = 4, nrow = 5)  
dev.off()

################################################################################
# Classify using median values and paired wilcoxon test
################################################################################
# source("scripts/classify_proteins_wilcox_delta.R")
# 
# prots <- signif$prot
# classified <- classify_median_wilcox(d_wide, prots, wilcox_p_threshold = 0.05)
# cl <- classified$Pattern
# names(cl) <- classified$prot
# plot_clusters(cl, prot_trajs = prot_trajs, out_path = paste0(out_basedir, "/trajectories_gam/classify_median_wilcox_", length(all_prots), "_prots.pdf"), add_cluster_name = T)  
# 
# classified$test_pattern <- classified$Pattern
# classified$test_pattern <- gsub("-", "same",classified$test_pattern)
# classified$test_pattern <- gsub("/", "up",classified$test_pattern)
# classified$test_pattern <- gsub("\\", "down",classified$test_pattern,  fixed=TRUE)
# classified$test_pattern <- gsub(" ", "_",classified$test_pattern)
# write.table(classified, file = paste0(out_basedir, "/trajectories_gam/classification_median_wilcox_", length(all_prots), "_prots.txt"), quote = F, sep = "\t", row.names = FALSE)

################################################################################
# Run trajectory analysis using LMM
################################################################################
source("scripts/run_trajectory_analysis_lmm.R")
run_trajectory_analysis_lmm(d_wide, gam_res_prot_tp[gam_res_prot_tp$lmm_BH_pval < 0.05, "prot"], out_basedir, n_points = 100)

################################################################################
# Protein - protein association 
################################################################################
 
#(? do we need it ?)

################################################################################
# Differentially expressed proteins between visits
################################################################################
library("limma")

run_limma<-function(joined_data, tp1, tp2) {
  df <-joined_data[joined_data$TP %in% c(tp1, tp2),]
  df$ID <- as.factor(df$ID)
  df$SampleID <- NULL
  
  # design a model 
  design<-model.matrix(~0 + as.factor(df$TP) + df$Age + df$BMI + df$storage_months + df$Pregnancy_category)
  colnames(design)<-c("TP1", "TP2", "Age", "BMI", "storage_months", "Pregnancy_category")
  
  # specify the pairing
  corfit <- duplicateCorrelation(t(df[,7:ncol(df)]), design, block = df$ID)
  
  # make contrast - what to compare
  contrast<- makeContrasts(Diff = TP1 - TP2, levels=design)
  
  # apply linear model to each protein
  # Robust regression provides an alternative to least squares regression that works with less restrictive assumptions. Specifically, it provides much better regression coefficient estimates when outliers are present in the data
  fit<-lmFit(t(df[,7:ncol(df)]), design=design,  method="robust", correlation =
               corfit$consensus )
  # apply contrast
  contrast_fit<-contrasts.fit(fit, contrast)
  # apply empirical Bayes smoothing to the SE
  ebays_fit<-eBayes(contrast_fit)
  # summary
  print(summary(decideTests(ebays_fit)))
  # extract DE results
  DE_results<-topTable(ebays_fit, n=(ncol(df) - 8), adjust.method="BH", confint=TRUE)
  
  return(DE_results)
}

run_wilcox <- function(joined_data_adj_covar, tp1, tp2) {
  joined_data_adj_covar$SampleID <- NULL
  wilcox_pvals <- data.frame(matrix(ncol = 3))
  colnames(wilcox_pvals) <- c("TP1_TP2", "prot", "wilcox_pval")
  cnt <- 1
  for (prot in colnames(joined_data_adj_covar)[7:ncol(joined_data_adj_covar)]){
    df <-joined_data_adj_covar[joined_data_adj_covar$TP %in% c(tp1, tp2), c("ID", "TP", prot)]
    df_wide <- na.omit(my_pivot_wider(df, row_names = "ID", names_from = "TP", values_from = prot))
    pval <- wilcox.test(df_wide[,1], df_wide[,2], paired = T)$p.value
    wilcox_pvals[cnt,] <- c(paste0(tp1, "_", tp2), prot, pval)
    cnt <- cnt + 1
  }
  wilcox_pvals$wilcox_pval <- as.numeric(wilcox_pvals$wilcox_pval)
  wilcox_pvals$BH_qval <- p.adjust(wilcox_pvals$wilcox_pval, method = 'BH')
  return(wilcox_pvals)
}

limma_res_all <- data.frame()
wilcox_res_all <- data.frame()
for (tp1 in 1:4){
  for (tp2 in (tp1 + 1):4){
    if(tp2 < 5 & tp1 != tp2){
      tp1 <- as.character(tp1)
      tp2 <- as.character(tp2)
      
      # limma
      limma_res <- run_limma(joined_data, tp1, tp2) %>%
        rownames_to_column(var = 'prot')
      if (nrow(limma_res[limma_res$adj.P.Val < 0.05,]) > 0) {
        limma_res_all <- rbind(limma_res_all, cbind(paste0(tp1, "_", tp2), limma_res))
      }
      
      # wilcoxon
      wilcox_res <- run_wilcox(joined_data_adj_covar, tp1, tp2)
      wilcox_res_all <- rbind(wilcox_res_all, wilcox_res)
    }
  }
}
colnames(limma_res_all)[1] <- "TP1_TP2"
limma_res_all$sign <- ifelse(limma_res_all$adj.P.Val < 0.05, T, F)

all(limma_res_all[limma_res_all$adj.P.Val < 0.05, "prot"] %in% signif$prot)
write.table(limma_res_all, file = paste0(out_basedir, "limma_DEPs.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(wilcox_res_all, file = paste0(out_basedir, "wilcox_DEPs.txt"), quote = F, sep = "\t", row.names = FALSE)

# Plot significant DEPs

# stacked barplot
results <- limma_res_all %>%
  mutate(
    Regulation = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Up-regulated",
      adj.P.Val < 0.05 & logFC < 0 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

summary_data <- results %>%
  filter(Regulation != "Not significant") %>%  # Exclude non-significant proteins
  group_by(TP1_TP2, Regulation) %>%
  summarize(Count = n(), .groups = "drop")

pdf(paste0(out_basedir, "limma_DEPs_barplot.pdf"))
ggplot(summary_data, aes(x = TP1_TP2, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Timepoint Comparison",
    y = "Number of DEP",
    fill = "Efect direction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = my_colors[c(3,2)])
dev.off()

################################################################################
# GAM  phenotype vs TP
################################################################################

pheno <- read.delim("../phenotypes/blood_pheno_13122024_log_adj_batch_storage_rm_outiers_4sds.txt.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))
pheno_adjusted <- regress_covariates_lmm(pheno, subset(covariates, select = -storage_months), covars_longitudinal = F)
## Hormone and lipid trajectories
all_pheno <- colnames(pheno)[4:ncol(pheno)]

# pheno_cov_table <- data.frame(matrix(nrow = length(all_pheno) * 3, ncol = 3))
# cnt <-1
# for (ph in all_pheno){
#   for(cov in c("Age", "BMI", "Pregnancy_category")){
#    d_subs <- inner_join(pheno[,c("ID", "TP", ph)], covariates[, c("ID", cov)], by = "ID")
#    colnames(d_subs) <- c("ID", "TP", "pheno", "cov")
#    model <- lmer(pheno ~ cov + (1|ID), data = d_subs)
#    coef <- summary(model)$coefficients
#    row.names(coef) <- gsub("cov2", "cov", row.names(coef))
#    pval <- coef["cov", "Pr(>|t|)"]
#    pheno_cov_table[cnt,] <- c(ph, cov, pval)
#    cnt <- cnt + 1
#   }
# }
# pheno_cov_table$X3 <- as.numeric(pheno_cov_table$X3)
# pheno_cov_table$BH_pval <- p.adjust(pheno_cov_table$X3)

gam_res_pheno_tp <- data.frame(matrix(nrow = length(all_pheno), ncol = 7))
colnames(gam_res_pheno_tp) <- c("pheno", "gam_pval", "gam_edf", "gam_fval", "n", "n_samples", "lmm_pval")

plot_list <- list()
cnt <- 1
for (ph in all_pheno){
  res_gam <- gam_prot_tp_adj_covar(pheno, ph, subset(covariates, select = -storage_months), scale = T, predict = T)
  res_lmm <- lmm_prot_tp_poly3_adj_covar(pheno, ph, subset(covariates, select = -storage_months))
  gam_res_pheno_tp[cnt,] <- c(ph, unlist(res_gam[1:5]), res_lmm)
  cnt <- cnt + 1
  
  traj <- data.frame(TP = seq(1,4, length.out = length(res_gam$predicted)), pheno = res_gam$predicted, lower = res_gam$lower, upper = res_gam$upper)
  traj <- full_join(pheno_adjusted[,c("TP", "ID", ph)], traj, by = "TP")
  colnames(traj)[3] <- "values"
  traj$values <- scale(traj$values)
  plot_list[[ph]] <- ggplot(traj) + 
    geom_boxplot(aes(x = TP, y = values, group = TP), width = 0.1, color = my_colors[3], outliers = F) +
    geom_line(aes(x = TP, y = pheno), color = my_colors[4]) + 
    geom_ribbon(aes(x = TP, ymin = lower, ymax = upper), alpha = 0.2, fill = my_colors[4]) +
    geom_jitter(aes(x = TP, y = values), alpha = 0.3, width = 0.1, color = my_colors[3]) +
    labs(x = "Timepoint ", y = ph) +
    ggtitle(paste0(ph, "; GAM P = ", formatC(res_gam$pval, digits = 3))) + 
    theme_minimal() 
}
pdf(paste0(out_basedir,"/plots/lipids_hormones_gam_withpoints.pdf"), height = 15, width = 15)
grid.arrange(grobs = plot_list, ncol = 4, nrow = 4)  
dev.off()

gam_res_pheno_tp <- na.omit(gam_res_pheno_tp) %>%
  mutate(across(-c( pheno), as.numeric)) 

gam_res_pheno_tp$gam_BH_pval <- p.adjust(gam_res_pheno_tp$gam_pval, method = 'BH')
gam_res_pheno_tp$lmm_BH_pval <- p.adjust(gam_res_pheno_tp$lmm_pval, method = 'BH')

gam_res_pheno_tp <- gam_res_pheno_tp[order(gam_res_pheno_tp$gam_pval),]
gam_res_pheno_tp$gam_edf_round <- round(gam_res_pheno_tp$gam_edf)
gam_res_pheno_tp$gam_BH_sign <- ifelse(gam_res_pheno_tp$gam_BH_pval < 0.05,T,F)

write.table(gam_res_pheno_tp, file = paste0(out_basedir, "pheno_vs_tp_gam_adj_age_bmi_preg.txt"), quote = F, sep = "\t", row.names = FALSE)

################################################################################
# GAM protein vs phenotype 
################################################################################
lmm_res <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol =14))
colnames(lmm_res) <- c("prot", "pheno", 
                       paste( c("estimate", "pval", "se", "tval", "N", "N_unique"), "withTP", sep = "_"),
                       paste( c("estimate", "pval", "se", "tval", "N", "N_unique"), "noTP", sep = "_")
                      )

gam_res <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol = 17))
colnames(gam_res) <- c("prot", "pheno", 
                       paste(c("pval", "edf_round", "fval", "n", " n_samples"), "spline", sep =  "_"),
                       paste(c("pval", "edf_round", "fval", "n", " n_samples"), "linear", sep =  "_"),
                       paste(c("pval", "edf_round", "fval", "n", " n_samples"), "noTP", sep =  "_")
                      )

cnt <- 1
for (ph in colnames(pheno)[4:ncol(pheno)]) {
  cat(ph, "\n")
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res_gam <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'spline', anova_pval = F)
    res_gam_lin <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'linear', anova_pval = F)
    res_gam_noTP <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'none', anova_pval = F)
    
    gam_res[cnt,] <- c(prot, ph, unlist(res_gam), unlist(res_gam_lin), unlist(res_gam_noTP))
    
    res_lmm <- lmm_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'cubic')
    res_lmm_noTP <- lmm_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'none')
    lmm_res[cnt,] <- c(prot, ph, unlist(res_lmm), unlist(res_lmm_noTP))
    cnt <- cnt + 1
  }
}

gam_res <- na.omit(gam_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
gam_res$BH_pval_spline <- p.adjust(gam_res$pval_spline)
gam_res$BH_pval_linear <- p.adjust(gam_res$pval_linear)
gam_res$BH_pval_noTP <- p.adjust(gam_res$pval_noTP)
gam_res <- gam_res[order(gam_res$pval_spline),]

cat("Number of BH significant associations:\n")
cat(" - with s(TP):", nrow(gam_res[gam_res$BH_pval_spline < 0.05,]), "; of them non-linear:",  nrow(gam_res[gam_res$BH_pval_spline < 0.05 & gam_res$edf_round_spline > 1,]), "\n")
cat(" - with linear TP:", nrow(gam_res[gam_res$BH_pval_linear < 0.05,]), "; of them non-linear:",  nrow(gam_res[gam_res$BH_pval_linear < 0.05 & gam_res$edf_round_linear > 1,]), "\n")
cat(" - without correcting for TP:", nrow(gam_res[gam_res$BH_pval_noTP < 0.05,]), "; of them non-linear:",  nrow(gam_res[gam_res$BH_pval_noTP < 0.05 & gam_res$edf_round_noTP > 1,]), "\n")

write.table(gam_res, file = paste0(out_basedir, "prot_vs_pheno_spline_gam_adj_covar.txt"), quote = F, sep = "\t", row.names = FALSE)


lmm_res <- na.omit(lmm_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
lmm_res$BH_pval_withTP <- p.adjust(lmm_res$pval_withTP)
lmm_res$BH_pval_noTP <- p.adjust(lmm_res$pval_noTP)
lmm_res <- lmm_res[order(lmm_res$pval_withTP),]
nrow(lmm_res[lmm_res$BH_pval_withTP < 0.05,])
nrow(lmm_res[lmm_res$BH_pval_noTP < 0.05,])

write.table(lmm_res, file = paste0(out_basedir, "prot_vs_pheno_lmm_adj_covar.txt"), quote = F, sep = "\t", row.names = FALSE)

##### Linear effect of phenotype using GAMs #####
gam_res_lin <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol = 17))
colnames(gam_res_lin) <- c("prot", "pheno", 
                       paste(c("pval", "est", "se", "n", " n_samples"), "spline", sep =  "_"),
                       paste(c("pval", "est", "se", "n", " n_samples"), "linear", sep =  "_"),
                       paste(c("pval", "est", "se",  "n", " n_samples"), "noTP", sep =  "_")
)

cnt <- 1
for (ph in colnames(pheno)[4:ncol(pheno)]) {
  cat(ph, "\n")
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res_gam <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'spline', adjust_pheno = 'linear', anova_pval = F)
    res_gam_lin <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'linear', adjust_pheno = 'linear', anova_pval = F)
    res_gam_noTP <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'none', adjust_pheno = 'linear', anova_pval = F)
    
    gam_res_lin[cnt,] <- c(prot, ph, unlist(res_gam), unlist(res_gam_lin), unlist(res_gam_noTP))
    
    cnt <- cnt + 1
  }
}

gam_res_lin <- na.omit(gam_res_lin) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
gam_res_lin$BH_pval_spline <- p.adjust(gam_res_lin$pval_spline)
gam_res_lin$BH_pval_linear <- p.adjust(gam_res_lin$pval_linear)
gam_res_lin$BH_pval_noTP <- p.adjust(gam_res_lin$pval_noTP)
gam_res_lin <- gam_res_lin[order(gam_res_lin$pval_spline),]

cat("Number of BH significant associations:\n")
cat(" - with s(TP):", nrow(gam_res_lin[gam_res_lin$BH_pval_spline < 0.05,]), "\n")
cat(" - with linear TP:", nrow(gam_res_lin[gam_res_lin$BH_pval_linear < 0.05,]), "\n")
cat(" - without correcting for TP:", nrow(gam_res_lin[gam_res_lin$BH_pval_noTP < 0.05,]), "\n")

write.table(gam_res_lin, file = paste0(out_basedir, "prot_vs_pheno_linear_gam_adj_covar.txt"), quote = F, sep = "\t", row.names = FALSE)



#### Comparisons ####
d1 <- lmm_res
gr1 = 'noTP'
cols1 <-colnames(d1)[grepl(gr1, colnames(d1))]

d2 <- gam_res_lin
gr2 <- "noTP"
cols2 <-colnames(d2)[grepl(gr2, colnames(d2))]

tmp <- full_join(d1[,c("prot", "pheno", cols1)], d2[,c("prot", "pheno", cols2)], by = c("prot", "pheno"))
#tmp[tmp[,paste0("pval_", gr1, ".x")] == 0,  paste0("pval_", gr1, ".x")] = 1e-08
c1 <- sym(paste0("pval_", gr1, ".x"))
c2 <- sym(paste0("pval_", gr2, ".y"))

ggplot(tmp, aes(x = -log10(!!c1), y = -log10(!!c2))) + 
  geom_point() +
  theme_minimal() +
  xlab("LMM visit: none") + 
  ylab("GAM phenotype: linear, visit: none") 


cat("BH significant in group 1, but not in group 2:\n")
tmp[tmp[,paste0("BH_pval_", gr1, ".x")] < 0.05 & tmp[,paste0("BH_pval_", gr2, ".y")] > 0.05, c("prot", "pheno")]

cat("BH significant in group 2, but not in group 1:\n")
tmp[tmp[,paste0("BH_pval_", gr1, ".x")] > 0.05 & tmp[,paste0("BH_pval_", gr2, ".y")] < 0.05, c("prot", "pheno")]





