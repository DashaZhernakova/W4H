
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/data/")
source("../scripts/traj_functions.R")

library(ggplot2)
library(rmcorr)
library(dplyr)
library(lme4)
library(grid)
library(gridExtra)
library(pheatmap)
library(corrplot)

pheno <- read.delim("../../phenotypes/blood_pheno_03102024_log.txt", as.is = T, check.names = F, sep = "\t")
d_wide <- read.delim("olink_clean_CVD+INF.txt", as.is = T, check.names = F, sep = "\t")
d_long <- read.delim("olink_clean_CVD+INF_long.txt", as.is = T, check.names = F, sep = "\t")
covariates <- read.delim("../data/covariates_age_bmi_storage_preg.txt", sep = "\t", check.names = F, as.is = T, colClasses = c(ID = "character"))

ID <- gsub("_.*", "", d_wide$SampleID)
TP <- gsub(".*_", "", d_wide$SampleID)

d_wide <- cbind(ID, TP, d_wide)
#d_wide$SampleID <- NULL

pheno$Record.ID <- gsub("ID_", "", pheno$Record.ID)
colnames(pheno)[1] <- "ID"
pheno <- cbind(SampleID = paste0(pheno$ID, "_", pheno$TP), pheno)
pheno$Age <- NULL

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

tmp1 <- cbind(1, regress_covariates(d_wide[d_wide$TP == 1,], covariates))
tmp2 <- cbind(2, regress_covariates(d_wide[d_wide$TP == 2,], covariates))
tmp3 <- cbind(3, regress_covariates(d_wide[d_wide$TP == 3,], covariates))
tmp4 <- cbind(4, regress_covariates(d_wide[d_wide$TP == 4,], covariates))
colnames(tmp1)[1] = colnames(tmp2)[1]= colnames(tmp3)[1] = colnames(tmp4)[1] = 'TP'
d_wide_adj_covar <- rbind(tmp1, tmp2, tmp3, tmp4)
joined_data_adj_covar <- full_join(covariates, d_wide_adj_covar, by = c("ID"), relationship = 'one-to-many')



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
ggplot(icc[1:100,], aes(x = reorder(prot, -ICC), y = ICC)) + 
  geom_point(size = 2) + 
  theme_minimal() + xlab("protein")

ggplot(d_wide, aes(x = TP, y = PROK1, group = ID)) + 
  geom_line(alpha = 0.2) + 
  geom_point(alpha = 0.4) + 
  theme_minimal() + 
  theme(legend.position = "")

################################################################################
# Protein vs TP LMM
################################################################################

lmm_res_prot_tp <- data.frame(matrix(nrow = (ncol(d_wide) -4), ncol = 2))
colnames(lmm_res_prot_tp) <- c("prot", "pval")
cnt <- 1
for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res <- lmm_prot_tp_poly3_adj_covar(d_wide, prot, covariates)
    lmm_res_prot_tp[cnt,] <- c(prot, res)
    cnt <- cnt + 1
}

lmm_res_prot_tp <- na.omit(lmm_res_prot_tp) %>%
  mutate(across(-c( prot), as.numeric)) 

lmm_res_prot_tp$BH_pval <- p.adjust(lmm_res_prot_tp$pval, method = 'BH')
lmm_res_prot_tp <- lmm_res_prot_tp[order(lmm_res_prot_tp$pval),]
lmm_res_prot_tp$bonf_sign <- ifelse(lmm_res_prot_tp$pval < 0.05/38,T,F)

write.table(lmm_res_prot_tp, file = "../results/prot_vs_tp_poly3_lmm_adj_age_bmi_preg_storage.txt", quote = F, sep = "\t", row.names = FALSE)
signif <- lmm_res_prot_tp[lmm_res_prot_tp$BH_pval < 0.05,]
nrow(signif)

lmm_res_prot_tp <- read.delim("../results/prot_vs_tp_poly3_lmm_adj_age_bmi_preg_storage.txt", sep = "\t", as.is = T, check.names = F)

# plots
prot= 'CHRDL2'
prot= 'PROK1'
prot_ensym <- ensym(prot)
ggplot(d_wide, aes(x = TP, y = !!prot_ensym)) + 
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.3) + 
  theme_bw() +
  ggtitle(paste0(prot, ", P = ", formatC(lmm_res_prot_tp[lmm_res_prot_tp$prot == prot, 2], digits = 3))) 

# radian plot

pdf("../plots/radian_BH.pdf", width = 10, height = 10)
make_radian_plot(d_wide, signif$prot)
dev.off()



################################################################################
# Cluster longitudinal trajectories
################################################################################


n_points = 100

all_prots <- colnames(d_wide)[4:ncol(d_wide)]
all_prots <- signif$prot

prot_coefs <- data.frame(matrix(nrow = length(all_prots) , ncol = 3))
row.names(prot_coefs) <- all_prots

prot_coefs_raw <- data.frame(matrix(nrow = length(all_prots) , ncol = 3))
row.names(prot_coefs_raw) <- all_prots

prot_trajs <- data.frame(matrix(nrow = length(all_prots) , ncol = n_points))
row.names(prot_trajs) <- all_prots
colnames(prot_trajs) <- seq(1,4, length.out = n_points)

prot_dist <- data.frame(matrix(nrow = length(all_prots) * length(all_prots), ncol = 4))
colnames(prot_dist) <- c("prot1", "prot2", "eucl",  "eucl_coef")

prot1_passed <- ""

pb = txtProgressBar(min = 0, max = length(all_prots), initial = 0) 
stepi = 0
cnt <- 1
for (prot1 in all_prots){
  #print(prot1)
  lmm_fit1 <- fit_lmm_poly3_adj_covar(d_wide, prot1, n = n_points, covariates, scale = T)
  prot_coefs[prot1,] <- lmm_fit1$coefficients
  prot_trajs[prot1,] <- lmm_fit1$predicted
  
  lmm_fit_raw <- fit_lmm_poly3_adj_covar(d_wide, prot1, n = n_points, covariates, scale = T, poly_raw = T)
  prot_coefs_raw[prot1,] <- lmm_fit_raw$coefficients
  setTxtProgressBar(pb,stepi)
  stepi <- stepi + 1
#  for (prot2 in all_prots){
#    if (prot2 %in% prot1_passed) next
#  
#   lmm_fit2 <- fit_lmm_poly3_adj_covar(d_wide, prot2, n = n_points, covariates, scale = T)
#     
#    # Euclidean distance between trajectories
#    eucl_dist <- sum(abs(as.numeric(lmm_fit1$predicted) - as.numeric(lmm_fit2$predicted)))/length(lmm_fit1$predicted)
#    
#    #Euclidean distance in the lm coefficients space
#    eucl_dist_coef <- TSdist::EuclideanDistance(lmm_fit1$coefficients, lmm_fit2$coefficients)
#    
#    prot_dist[cnt,] <- c(prot1, prot2, eucl_dist, eucl_dist_coef)
#    cnt <- cnt + 1
#  }
#  prot1_passed <- c(prot1_passed, prot1)
}

cat("\n")
close(pb)

prot_dist <- na.omit(prot_dist) %>%
  mutate(across(-c( prot1, prot2), as.numeric)) 
dist_matrix <- my_pivot_wider(prot_dist, "prot1", "prot2", "eucl")
dist_matrix[lower.tri(dist_matrix)] <- t(dist_matrix)[lower.tri(dist_matrix)]

dist_matrix_coef <- my_pivot_wider(prot_dist, "prot1", "prot2", "eucl_coef")
dist_matrix_coef[lower.tri(dist_matrix_coef)] <- t(dist_matrix_coef)[lower.tri(dist_matrix_coef)]


#convert distance to similarity metric
prot_dist$eucl_similarity <- 1 / (1 + prot_dist$eucl)
prot_dist$eucl_similarity_coef <- 1 / (1 + prot_dist$eucl_coef)

sim_matrix <- my_pivot_wider(prot_dist, "prot1", "prot2", "eucl_similarity")
sim_matrix[lower.tri(sim_matrix)] <- t(sim_matrix)[lower.tri(sim_matrix)]

sim_matrix_coef <- my_pivot_wider(prot_dist, "prot1", "prot2", "eucl_similarity_coef")
sim_matrix_coef[lower.tri(sim_matrix_coef)] <- t(sim_matrix_coef)[lower.tri(sim_matrix_coef)]


write.table(sim_matrix, file = "../results/trajectories/similarity_matrix_eucl_109signif.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(dist_matrix, file = "../results/trajectories/distance_matrix_eucl_109signif.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(prot_dist, file = "../results/trajectories/protein_distance_eucl_109signif.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(prot_trajs, file = "../results/trajectories/protein_trajectories_109signif.txt", quote = F, sep = "\t", row.names = FALSE)


# Heatmap and clustering on euclidean distance between curves
method = 'hclust_my_eucl_dist'
num_k = 8
hm <- pheatmap(dist_matrix, cutree_rows = num_k, cutree_cols = num_k, filename = paste0("../plots/trajectories/", method, "_k", num_k, "_heatmap.pdf"))
cl  = cutree(hm$tree_row, k = num_k)
plot_clusters(cl, method, num_k, out_path = paste0("../plots/trajectories/", method, "_k", num_k, "_109signif.pdf"))
write.table(as.data.frame(cl), file = paste0("../results/trajectories/", method, "_k", num_k, "_109signif.txt"), quote = F, sep = "\t")

# Get pairwise differences between i and i+1 timepoint and cluster according to them
prot_deltas <- data.frame(matrix(nrow = nrow(prot_trajs), ncol = ncol(prot_trajs) - 1))

row.names(prot_deltas) <- row.names(prot_trajs)

for (c in 1:(ncol(prot_trajs) - 1)){
  prot_deltas[,c] <- prot_trajs[,(c + 1)] - prot_trajs[,c]
}
method = 'pam_delta'
num_k = 6
#num_k = 15
pam_res <- cluster::pam(prot_deltas, num_k, diss = F)
cl <- pam_res$clustering
plot_clusters(cl, method, num_k, out_path = paste0("../plots/trajectories/", method, "_k", num_k, "_109signif.pdf"))
write.table(as.data.frame(cl), file = paste0("../results/trajectories/", method, "_k", num_k, "_109signif.txt"), quote = F, sep = "\t")

pheatmap(prot_deltas, cluster_rows = F, cluster_cols = F,filename = paste0("../plots/trajectories/", method, "_k", num_k, "_heatmap.pdf"), fontsize_row = 6)

################################################################################
# Classify using median values and paired wilcoxon test
################################################################################
source("../scripts/classify_proteins_wilcox_delta.R")
prots <- colnames(d_wide)[4:ncol(d_wide)]
prots <- signif$prot
classified <- classify_median_wilcox(d_wide, prots, wilcox_p_threshold = 0.05)
cl <- classified$Pattern
names(cl) <- classified$prot
plot_clusters(cl, out_path = "../plots/trajectories/classify_median_wilcox_109signif.pdf", add_cluster_name = T)  

classified$test_pattern <- classified$Pattern
classified$test_pattern <- gsub("-", "same",classified$test_pattern)
classified$test_pattern <- gsub("/", "up",classified$test_pattern)
classified$test_pattern <- gsub("\\", "down",classified$test_pattern,  fixed=TRUE)
classified$test_pattern <- gsub(" ", "_",classified$test_pattern)
write.table(classified, file = "../results/trajectories/classification_median_wilcox_109signif.txt", quote = F, sep = "\t", row.names = FALSE)
 

################################################################################
# LMM-based association between proteins
################################################################################

all_prots <-  colnames(d_wide)[4:ncol(d_wide)]
#prot_rmcorr <- data.frame(matrix(nrow = length(all_prots) * length(all_prots), ncol = 4))
#colnames(prot_rmcorr) <- c("prot1", "prot2", "r", "p")

prot_lmm <- data.frame(matrix(nrow = length(all_prots) * length(all_prots), ncol = 6))
colnames(prot_lmm) <- c("prot1", "prot2", "r", "p", "r_tp", "p_tp")

prot1_passed <- ""

pb = txtProgressBar(min = 0, max = length(all_prots), initial = 0) 
stepi = 0
cnt <- 1
for (prot1 in all_prots){
  #prot1_passed <- c(prot1_passed, prot1)
  setTxtProgressBar(pb,stepi)
  stepi <- stepi + 1
  #print(prot1)
  for (prot2 in all_prots){
    if (prot2 %in% prot1_passed) next
    if (prot1 == prot2) next
    #rm_corr <- get_rmcorr(d_wide_adj_covar, d_wide_adj_covar, prot1, prot2) 
    #prot_rmcorr[cnt,] <- c(prot1, prot2, rm_corr$r, rm_corr$p)
    
    lmm_res <- lmm_pheno_prot_adj_covar_noTP(d_wide, d_wide, prot1, prot2, covariates, scale = T)
    lmm_res_tp <- lmm_pheno_prot_adj_covar_withTP(d_wide, d_wide, prot1, prot2, covariates, scale = T)
    prot_lmm[cnt,] <- c(prot1, prot2, lmm_res$estimate, lmm_res$pval, lmm_res_tp$estimate, lmm_res_tp$pva)
    cnt <- cnt + 1
  }
  prot1_passed <- c(prot1_passed, prot1)
}
cat("\n")
close(pb)

prot_lmm <- na.omit(prot_lmm) %>%
  mutate(across(-c( prot1, prot2), as.numeric))  

write.table(prot_lmm, file = "../results/prot_vs_prot_lmm_adj_covar.txt", quote = F, sep = "\t", row.names = FALSE)


prot_lmm <- read.delim("../results/prot_vs_prot_lmm_adj_covar.txt", as.is = T, check.names = F, sep = "\t")
prot_lmm <- na.omit(prot_lmm) %>%
  mutate(across(-c( prot1, prot2), as.numeric))  

# fill in proteins only tested as prot1 or prot2 because of the loop check
missing_prot2 <- setdiff(unique(prot_lmm$prot1), unique(prot_lmm$prot2))
missing_prot1 <- setdiff(unique(prot_lmm$prot2), unique(prot_lmm$prot1))
mp2_res <- prot_lmm[prot_lmm$prot1 == missing_prot2, c(2,1,3,4,5,6)]
mp1_res <- prot_lmm[prot_lmm$prot2 == missing_prot1, c(2,1,3,4,5,6)]
mp1_res <- mp1_res[mp1_res$prot1 != missing_prot2,]
colnames(mp1_res) = colnames(mp2_res) = colnames(prot_lmm)

prot_lmm <- rbind(prot_lmm, mp2_res, mp1_res)

lmm_matrix <- my_pivot_wider(prot_lmm, "prot1", "prot2", "r")
lmm_matrix[lower.tri(lmm_matrix)] <- t(lmm_matrix)[lower.tri(lmm_matrix)]

lmm_matrix <- lmm_matrix[, row.names(lmm_matrix)]
diag(lmm_matrix) <- 1
lmm_matrix[lmm_matrix > 1] <- 1
pdf("../plots/lmm_prot_prot_corrplot_2.pdf", height = 25, width = 25)
corrplot(as.matrix(lmm_matrix), tl.col = 'black', order = 'hclust', addgrid.col = NA, tl.cex = 0.2)
dev.off()

c <- corrplot::corrplot(as.matrix(lmm_matrix), tl.col = 'black', order = 'hclust', addgrid.col = NA)

c_ord <- row.names(c$corr)
st <- which(c_ord == 'PRDX3')
end <- which(c_ord == 'KYAT1')

st <- which(c_ord == 'CCL16')
end <- which(c_ord == 'DKK3')
sel <- c_ord[st:end]


plot_traj_many_prots2(prot_trajs, sel, colored = F)

#
# Corrplot per TP
#

cor_1 <- cor(d_wide[d_wide$TP == 1, -c(1,2,3)], use = 'complete.obs')
cor_2 <- cor(d_wide[d_wide$TP == 2, -c(1,2,3)], use = 'complete.obs')
cor_3 <- cor(d_wide[d_wide$TP == 3, -c(1,2,3)], use = 'complete.obs')
cor_4 <- cor(d_wide[d_wide$TP == 4, -c(1,2,3)], use = 'complete.obs')

pdf("../plots/corrplot_prot_prot_per_tp.pdf")
c <- corrplot(cor_1, tl.col = 'black', tl.cex = 0.1, title = "TP = 1", addgrid.col = NA, order = 'hclust')
ord <- row.names(c$corr)
corrplot(cor_2[ord,ord], tl.col = 'black', tl.cex = 0.1, title = "TP = 2", addgrid.col = NA)
corrplot(cor_3[ord,ord], tl.col = 'black', tl.cex = 0.1, title = "TP = 3", addgrid.col = NA)
corrplot(cor_4[ord,ord], tl.col = 'black', tl.cex = 0.1, title = "TP = 4", addgrid.col = NA)
dev.off()



#
# Cluster based lmm estimate
#

method = 'pam_lmm'
num_k = 6
#num_k = 15
pam_res <- cluster::pam(lmm_matrix, num_k, diss = F)
cl <- pam_res$clustering

write.table(as.data.frame(cl), file = paste0("../results/trajectories/", method, "_k", num_k, "_all.txt"), quote = F, sep = "\t")



################################################################################
# Differentially expressed proteins between visits
################################################################################
library("limma")

run_limma<-function(joined_data, tp1, tp2) {
  df <-joined_data[joined_data$TP %in% c(tp1, tp2),]
  df$ID <- as.factor(df$ID)
  
  
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

# heatmap
counts <- with(limma_res_all[limma_res_all_sign$sign == T,c("TP1_TP2", "prot")], table(prot, TP1_TP2))
pheatmap(counts, cluster_rows = F, cluster_cols = F, color = c( "white", "dodgerblue4"))


################################################################################
# LMM protein - phenotype
################################################################################


## Hormone and lipid trajectories
lmm_res_pheno_tp <- data.frame(matrix(nrow = (ncol(pheno) -4), ncol = 2))
colnames(lmm_res_pheno_tp) <- c("pheno", "pval")
cnt <- 1

plot_list <- list()
for (ph in colnames(pheno)[4:ncol(pheno)]) {
  ph_lmm <- lmm_prot_tp_poly3_adj_covar(pheno, ph, covariates, scale = T)
  lmm_res_pheno_tp[cnt,] <- c(ph, ph_lmm)
  cnt <- cnt + 1
  
  d_subs <- pheno[,c('TP', ph)]
  colnames(d_subs)[2] <- 'pheno'
  plot_list[[ph]] <- ggplot(d_subs, aes(x = TP, y = pheno)) +
    geom_smooth(method = 'lm', formula=y ~ poly(x, 3, raw=TRUE), aes(x = TP, y = pheno)) + 
    labs(x = "Timepoint ", y = ph) +
    ggtitle(paste0(ph, "; LMM P = ", formatC(ph_lmm, digits = 3))) + 
    theme_minimal()
}
pdf("../plots/lipids_hormones_poly3.pdf", height = 15, width = 15)
grid.arrange(grobs = plot_list, ncol = 4, nrow = 4)  
dev.off()

lmm_res_pheno_tp <- na.omit(lmm_res_pheno_tp) %>%
  mutate(across(-pheno, as.numeric))  
lmm_res_pheno_tp$BH_pval <- p.adjust(lmm_res_pheno_tp$pval, method = 'BH')
write.table(lmm_res_pheno_tp, file = "../results/pheno_vs_tp_poly3_lmm_adj_age_bmi_preg_storage.txt", quote = F, sep = "\t", row.names = FALSE)



#### Proteins vs phenotypes: LMM

#lmm_res <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol = 6))
#colnames(lmm_res) <- c("prot", "pheno", "estimate", "pval", "se", "tval")

lmm_res_no_tp <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol = 6))
colnames(lmm_res_no_tp) <- c("prot", "pheno", "estimate", "pval", "se", "tval")
cnt <- 1
for (ph in colnames(pheno)[4:ncol(pheno)]) {
  cat(ph, "\n")
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    #res <- lmm_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'cubic')
    #lmm_res[cnt,] <- c(prot, ph, unlist(res))
    
    res0 <- lmm_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'none')
    lmm_res_no_tp[cnt,] <- c(prot, ph, unlist(res0))
    cnt <- cnt + 1
  }
}

lmm_res <- na.omit(lmm_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
lmm_res$BH_pval <- p.adjust(lmm_res$pval)
lmm_res <- lmm_res[order(lmm_res$pval),]
nrow(lmm_res[lmm_res$BH_pval < 0.05,])
write.table(lmm_res, file = "../results/prot_vs_pheno_withTP_lmm_adj_covar.txt", quote = F, sep = "\t", row.names = FALSE)


lmm_res_no_tp <- na.omit(lmm_res_no_tp) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
lmm_res_no_tp$BH_pval <- p.adjust(lmm_res_no_tp$pval)
lmm_res_no_tp <- lmm_res_no_tp[order(lmm_res_no_tp$pval),]
write.table(lmm_res_no_tp, file = "../results/prot_vs_pheno_noTP_lmm_adj_covar.txt", quote = F, sep = "\t", row.names = FALSE)

lmm_res <- read.delim("../results/prot_vs_pheno_withTP_lmm_adj_covar.txt", sep = "\t", as.is = T, check.names = T)
lmm_res_no_tp <- read.delim("../results/prot_vs_pheno_noTP_lmm_adj_covar.txt", sep = "\t", as.is = T, check.names = T)

# tmp <- inner_join(lmm_res, inner_join(lmm_res_linTP, lmm_res_tp, by = c("pheno", "prot")), by = c("pheno", "prot"))
# tmp$logp.x <- -log10(tmp$pval.x)
# tmp$logp.y <- -log10(tmp$pval.y)
# tmp$logp <- -log10(tmp$pval)
# 
# ggplot(tmp, aes(x = logp, y = logp.x, color = pheno)) + geom_point(alpha = 0.8)  + 
#   theme_minimal() + xlab("no TP") + ylab("linear TP") + xlim(0, 20) + ylim(0,20) + 
#   geom_vline(xintercept = min (tmp[tmp$BH_pval < 0.05, "logp"]), linetype="dashed", color = 'grey') + 
#   geom_hline(yintercept = min (tmp[tmp$BH_pval.x < 0.05, "logp.x"]), linetype="dashed", color = 'grey')


# LM per TP

lm_res <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol = 10))
colnames(lm_res) <- c("prot", "pheno", "b_1", "p_1", "b_2", "p_2", "b_3", "p_3","b_4", "p_4")
cnt <- 1
for (ph in colnames(pheno)[4:ncol(pheno)]) {
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    lm_res[cnt,] <- lm_per_tp_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T)
    cnt <- cnt + 1
  }
}

lm_res <- na.omit(lm_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 

lm_res$logp <- -log10(lm_res$pval)
lm_res <- inner_join(lm_res, lmm_res, by = c("prot", "pheno"))



# Scatter prot vs pheno
plot_list <- list()
for (i in 1:20) {
  plot_list[[i]] <- scatter_col_tp(d_wide, pheno, lmm_res[i, "prot"], lmm_res[i, "pheno"], scale = T)
}
pdf("../plots/prot_pheno_scatter_scaled_top.pdf", height = 20, width = 15)
grid.arrange(grobs = plot_list, ncol = 4, nrow = 5)  
dev.off()


# Heatmap
lmm_res_est_mat <- my_pivot_wider(lmm_res, "prot", "pheno", "estimate")


lmm_res_est_mat[lmm_res_est_mat > 1] <- 1
lmm_res_est_mat[lmm_res_est_mat < -1] <- -1

breaksList = seq(-1, 1, by = 0.01)
colorList <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))

pheatmap(as.matrix(lmm_res_est_mat), treeheight_row = 0, show_rownames = F)
pheatmap(as.matrix(lmm_res_est_mat), treeheight_row = 0, show_rownames = F, color = colorList, breaks = breaksList, filename = "../plots/pheno_vs_prot_heatmap_with_poly3TP.pdf")
dev.off()

#barplot
counts <- as.data.frame(table(lmm_res[lmm_res$BH_pval < 0.05, "pheno"]))
ggplot(counts, aes(x = reorder(Var1, -Freq), y = Freq)) + geom_bar(stat="identity", fill = my_colors[4]) + 
  theme_minimal() + xlab("Phenotype") + ylab("Number of BH significant protein associations")


#forest plot
lmm_res$lower <- lmm_res$estimate - 1.96*lmm_res$se
lmm_res$upper <- lmm_res$estimate + 1.96*lmm_res$se

lmm_res$sign <- ifelse(lmm_res$BH_pval <0.05, T, F)

sign_prots <- unique(lmm_res[lmm_res$BH_pval < 0.05, "prot"])
sign_pheno <- unique(lmm_res[lmm_res$BH_pval < 0.05, "pheno"])

tmp <- lmm_res[lmm_res$prot %in% sign_prots & lmm_res$pheno %in% sign_pheno,]
tmp$pheno <- factor(tmp$pheno, levels = c("ALT", "AST", "LH", "FSH","PROG", "X17BES", "PRL", "HDLC", "TotChol", "Trigl"))
pdf("../plots/pheno_vs_prot_forest_with_poly3TP.pdf", height = 10, width = 10)
ggplot(tmp ,aes(x = prot, y = estimate, ymin = lower, ymax = upper))+
  geom_pointrange(aes(color = sign)) +
  geom_hline(lty=2, aes(yintercept=0)) +
  facet_wrap(~pheno, ncol = 4, scales="free")+
  coord_flip() +
  scale_color_manual(values = c("grey", "black")) + 
  theme_minimal()  +
  theme(legend.position="none", axis.text.y=element_text(size=7))
dev.off()

# combined with no TP correction

tmp <- rbind(cbind(lmm_res, data.frame("adjustment" = "cubic TP")),
             cbind(lmm_res_no_tp, data.frame("adjustment" = "no TP")))

tmp$lower <- tmp$estimate - 1.96*tmp$se
tmp$upper <- tmp$estimate + 1.96*tmp$se

tmp$sign <- ifelse(tmp$BH_pval <0.05, "1", "2")

sign_prots <- unique(tmp[tmp$BH_pval < 0.05, "prot"])
sign_pheno <- unique(tmp[tmp$BH_pval < 0.05, "pheno"])

tmp <- tmp[tmp$prot %in% sign_prots & tmp$pheno %in% sign_pheno,]
tmp$pheno <- factor(tmp$pheno, levels = c("ALT", "AST", "LH", "FSH","PROG", "X17BES", "PRL", "HDLC", "TotChol", "Trigl", "HOMA_IR", "INS"))
pdf("../plots/pheno_vs_prot_forest_with_and_without_poly3TP.pdf", height = 20, width = 10)
ggplot(tmp ,aes(x = prot, y = estimate, ymin = lower, ymax = upper, group = adjustment, color = adjustment))+
  geom_linerange( position = position_dodge(width = 0.5), aes(lty = sign, color = adjustment)) +
  geom_point(size = 2.5,  stroke = 0.5, position = position_dodge(width = 0.5), aes(shape = sign))+
  geom_hline(lty=2, aes(yintercept=0)) +
  facet_wrap(~pheno, ncol = 4, scales="free")+
  coord_flip() +
  theme_minimal()  +
  scale_shape_manual(values = c(16,21)) +
  scale_color_manual(values = my_colors[c(2,3)])+
  theme(legend.position="none", axis.text.y=element_text(size=7))
dev.off()

#### Replication

repl_ice <- read.delim("../../iceland_selected.txt", check.names = F, as.is = T, sep = "\t")
repl_ice <- na.omit(repl_ice) %>%
  mutate(across(c(beta,pval), as.numeric))  

repl_ukb <- read.delim("../../ukb_selected.txt", check.names = F, as.is = T, sep = "\t")
repl_ukb <- na.omit(repl_ukb) %>%
  mutate(across(c(beta,pval), as.numeric)) 

repl_fenland <- read.delim("../../fenland_selected.txt", check.names = F, as.is = T, sep = "\t")
repl_fenland <- na.omit(repl_fenland) %>%
  mutate(across(-c(gene), as.numeric)) %>%
  pivot_longer(cols = 2:8, names_to = "phenotype", values_to = "var_explained")


tmp <- left_join(lmm_res[lmm_res$pheno %in% c("AST", "ALT", "Trigl", "TotChol", "HDLC", "LDLC", "Glucose"),], 
                 full_join(repl_ice, repl_ukb, by = c("gene", "phenotype"), suffix = c("_ice", "_ukb")), 
                 by = c("prot" = "gene", "pheno" = "phenotype"))
tmp$beta <- tmp$estimate

tmp <- left_join(tmp, repl_fenland, by = c("prot" = "gene", "pheno" = "phenotype"))
tmp[tmp == 0] = 5e-324

# fill in the absent lines
tmp[is.na(tmp$beta_ice) & tmp$pheno %in% repl_ice$phenotype & tmp$prot %in% repl_ice$gene, "pval_ice"] <- 1
tmp[is.na(tmp$beta_ukb) & tmp$pheno %in% repl_ukb$phenotype & tmp$prot %in% repl_ukb$gene, "pval_ukb"] <- 1


tmp$logp <- -log10(tmp$pval)
tmp$logp_ice <- -log10(tmp$pval_ice)
tmp$logp_ukb <- -log10(tmp$pval_ukb)

tmp$BH_signif <- ifelse(tmp$BH_pval < 0.05, T, F)


ggplot(tmp, aes(x = logp, y = logp_ice, color = BH_signif)) + 
  geom_point() + 
  theme_minimal() + 
  scale_color_manual(values = my_colors[c(3,2)]) + xlab("logp_W4H")


ggplot(tmp, aes(beta, y = beta_ice, color = BH_signif)) + 
  geom_point() + 
  geom_vline(xintercept = 0, lty = 2) + geom_hline(yintercept = 0, lty = 2) +
  theme_minimal() + scale_color_manual(values = my_colors[c(3,2)]) + xlab("beta_W4H")

summary_table <- tmp %>%
  group_by(pheno) %>%
  summarise(
    W4H_signif = sum(BH_pval < 0.05, na.rm = TRUE),
    of_signif_repl_p0.05 = sum(BH_pval < 0.05 & pval_ukb < 0.05, na.rm = TRUE),
    of_signif_repl_p0.0001 = sum(BH_pval < 0.05 & pval_ukb < 0.0001, na.rm = TRUE)
  )
View(summary_table)

write.table(tmp, file = "../results_replication_lmm_prot_vs_pheno_cubic.txt", quote = F, sep = "\t", row.names = FALSE)

ggplot(tmp, aes(beta, y = var_explained, color = BH_signif)) + 
  geom_point() + 
  geom_vline(xintercept = 0, lty = 2) + geom_hline(yintercept = 0, lty = 2) +
  theme_minimal() + scale_color_manual(values = my_colors[c(3,2)]) + xlab("beta_W4H")


#
# Interaction analysis
#


inter_res <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol = 3))
colnames(inter_res) <- c("prot", "pheno", "pval")
cnt <- 1
for (ph in colnames(pheno)[4:ncol(pheno)]) {
  cat(ph, "\n")
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res <- lmm_prot_tp_interaction_pheno_adj_covar(d_wide, pheno, prot, ph , covariates, adjust_timepoint = 'linear')
    inter_res[cnt,] <- c(prot, ph, res)
    cnt <- cnt + 1
  }
}

inter_res <- na.omit(inter_res) %>%
  mutate(across(pval, as.numeric)) 


inter_res$BH_pval <- p.adjust(inter_res$pval, method = 'BH')


################################################################################
# LMM without TP correction 
################################################################################

# tmp <- inner_join(lmm_res, lmm_res_no_tp, by = c("prot", "pheno"))
# tmp$logp.x <- -log10(tmp$pval.x)
# tmp$logp.y <- -log10(tmp$pval.y)
# 
# bh_x <- min(tmp[tmp$BH_pval.x < 0.05, "logp.x"])
# bh_y <- min(tmp[tmp$BH_pval.y < 0.05, "logp.y"])
# 
# ggplot(tmp, aes(x = logp.x, y = logp.y, color = pheno)) + 
#   geom_point() + 
#   geom_vline(xintercept = bh_x, lty = 2) + geom_hline(yintercept = bh_y, lty = 2) +
#   xlab("with cubic TP correction") + ylab("no TP correction") +
#   theme_minimal() +
#   ggtitle("-log10(P) of protein - phenotype LMM association")


################################################################################
# LMM vs rmcorr vs corr per ID
################################################################################

all_prots <- colnames(d_wide)[4:ncol(d_wide)]
cmp_lmm_rmcorr <- data.frame(matrix(nrow = length(all_prots) * (ncol(pheno) -4), ncol = 8))
colnames(cmp_lmm_rmcorr) <- c("prot", "pheno", "lmm_no_tp_beta", "lmm_no_tp_pval", "lmm_lin_tp_beta", "lmm_lin_tp_pval", "rmcorr_r", "rmcorr_pval")
cnt <- 1
for (ph in colnames(pheno)[4:ncol(pheno)]) {
  cat(ph, "\n")
  for (prot in all_prots){
    #res <- lmm_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'cubic')
    #lmm_res[cnt,] <- c(prot, ph, unlist(res))
    
    lmm_no_tp <- lmm_pheno_prot_no_adj_covar(d_wide, pheno, prot, ph, scale = T, adjust_timepoint = 'none')
    lmm_lin_tp <- lmm_pheno_prot_no_adj_covar(d_wide, pheno, prot, ph, scale = T, adjust_timepoint = 'linear')
    rmcor <- get_rmcorr(d_wide, pheno, prot, ph)
    cmp_lmm_rmcorr[cnt,] <- c(prot, ph, lmm_no_tp$estimate, lmm_no_tp$pval, lmm_lin_tp$estimate, lmm_lin_tp$pval, rmcor$r, rmcor$p)
    cnt <- cnt + 1
  }
}

cmp_lmm_rmcorr <- na.omit(cmp_lmm_rmcorr) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 

ggplot(cmp_lmm_rmcorr, aes(x = lmm_no_tp_beta, rmcorr_r)) + geom_point() + theme_minimal()
ggplot(cmp_lmm_rmcorr, aes(x = log10(lmm_no_tp_pval), log10(rmcorr_pval))) + geom_point() + theme_minimal()

################################################################################
# Phenotype - protein distances
################################################################################

all_pheno <- lmm_res_pheno_tp[lmm_res_pheno_tp$BH_pval < 0.05, "pheno"]
all_prots <- colnames(d_wide)[4:ncol(d_wide)]
#all_prots <- lmm_res_prot_tp[lmm_res_prot_tp$BH_pval < 0.05, "prot"]

n_points = 100

pheno_dist <- data.frame(matrix(nrow = length(all_prots) * length(all_pheno), ncol = 4))
colnames(pheno_dist) <- c("prot", "pheno", "eucl",  "eucl_coef")

pb = txtProgressBar(min = 0, max = length(all_prots), initial = 0) 
stepi = 0
cnt <- 1
for (prot in all_prots){
  #print(prot1)
  lmm_fit1 <- fit_lmm_poly3_adj_covar(d_wide, prot, n = n_points, covariates, scale = T)
  
  setTxtProgressBar(pb,stepi)
  stepi <- stepi + 1
  for (ph in all_pheno){
     lmm_fit2 <- fit_lmm_poly3_adj_covar(pheno, ph, n = n_points, covariates, scale = T)
       
      # Euclidean distance between trajectories
      eucl_dist <- sum(abs(as.numeric(lmm_fit1$predicted) - as.numeric(lmm_fit2$predicted)))/length(lmm_fit1$predicted)
      
      #Euclidean distance in the lm coefficients space
      eucl_dist_coef <- TSdist::EuclideanDistance(lmm_fit1$coefficients, lmm_fit2$coefficients)
      
      #area_diff <- get_AUC_difference(as.numeric(lmm_fit1$predicted), as.numeric(lmm_fit2$predicted), n_points)
      pheno_dist[cnt,] <- c(prot, ph, eucl_dist, eucl_dist_coef)
      cnt <- cnt + 1
    }
}

cat("\n")
close(pb)

pheno_dist <- na.omit(pheno_dist) %>%
  mutate(across(-c( prot, pheno), as.numeric)) 

pheno_dist <- pheno_dist %>%
  group_by(pheno) %>%
  mutate(eucl_scaled = scale_this(eucl),
         eucl_coef_scaled = scale_this(eucl_coef))

write.table(pheno_dist, file = "../results/prot_vs_pheno_noTP_distance_adj_covar_signif_pheno.txt", quote = F, sep = "\t", row.names = FALSE)

dist_matrix_ph <- my_pivot_wider(pheno_dist, "prot", "pheno", "eucl_coef")
pheatmap(dist_matrix_ph, show_rownames = F, treeheight_row = 0)

d_wide$TP <- as.numeric(d_wide$TP)



# FIX!!!! SMTH WRONG  - plotting doesn't use trajectories but makes them itself
# and adjusting per tp is not correct
#Make a dataframe with phenotypes adjusted for all covariates per visit
covariate_names <- c("Age","BMI","storage_months", "Pregnancy_category")

joined_data <- full_join(covariates, pheno, by = c("ID"), relationship = 'one-to-many')

tmp1 <- cbind(1, regress_covariates(pheno[pheno$TP == 1,], covariates))
tmp2 <- cbind(2, regress_covariates(pheno[pheno$TP == 2,], covariates))
tmp3 <- cbind(3, regress_covariates(pheno[pheno$TP == 3,], covariates))
tmp4 <- cbind(4, regress_covariates(pheno[pheno$TP == 4,], covariates))
colnames(tmp1)[1] = colnames(tmp2)[1]= colnames(tmp3)[1] = colnames(tmp4)[1] = 'TP'
pheno_adj_covar <- rbind(tmp1, tmp2, tmp3, tmp4)


pdf("../plots/prot_pheno_area_between_top.pdf", height = 15, width = 10)
#pdf("../plots/prot_pheno_eucl_dist_bottom.pdf", height = 15, width = 10)
for (ph in all_pheno){
  subs <- pheno_dist[pheno_dist$pheno == ph,]
  subs <- as.data.frame(subs[order(subs$area_diff),])
  
  plot_list <- list()
  cnt <- 1
  for (i in 1:5){
    plot_list[[cnt]] <- plot_together(d_wide, pheno, subs[i, "prot"], subs[i, "pheno"], annot = paste0(" eucl d = ", formatC(subs[i, "area_diff"], digits = 3), " scaled d = ", formatC(subs[i, "eucl_scaled"], digits = 3)))
    plot_list[[cnt + 1]] <- plot_together(d_wide, pheno, subs[i, "prot"], subs[i, "pheno"], method = 'boxplot')
    #j = nrow(subs) - i
    #plot_list[[cnt]] <- plot_together(d_wide, pheno, subs[j, "prot"], subs[j, "pheno"], annot = paste0(" eucl_coef d = ", formatC(subs[j, "eucl_coef"], digits = 3)))
    #plot_list[[cnt + 1]] <- plot_together(d_wide, pheno, subs[(j - 1), "prot"], subs[(j - 1), "pheno"], method = 'boxplot')
    
    
    cnt <- cnt + 2
  }
  grid.arrange(grobs = plot_list, ncol = 2, nrow = 5)  
}

dev.off()

plot_list <- list()
dist_metrics <- "eucl_coef"
for (ph in all_pheno) {
  subs <- as.data.frame(pheno_dist[pheno_dist$pheno == ph & pheno_dist$prot %in% all_prots,])

  top <- subs[subs[, dist_metrics] < quantile(subs[, dist_metrics], 0.05), ]
  #top <- subs[subs[, dist_metrics] < -1, ]

  #plot_list[[ph]] <- plot_traj_prots_and_pheno(d_wide, pheno, top$prot, ph, title = ph)
  plot_list[[ph]] <- plot_medians_prots_and_pheno(d_wide, pheno, top$prot, ph, title = ph)
}
pdf(paste0("../plots/trajectories/pheno_prot_medians_", dist_metrics, "_q0.05.pdf"), height = 15, width = 15)
grid.arrange(grobs = plot_list, ncol = 4, nrow = 4)  
dev.off()
