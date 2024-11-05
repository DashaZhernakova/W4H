
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/data/")
source("../scripts/traj_functions.R")

library(ggplot2)
library(rmcorr)
library(dplyr)
library(lme4)
library(grid)
library(gridExtra)

pheno <- read.delim("../../phenotypes/blood_pheno_03102024_log.txt", as.is = T, check.names = F, sep = "\t")
d_wide <- read.delim("olink_clean_CVD+INF.txt", as.is = T, check.names = F, sep = "\t")
d_long <- read.delim("olink_clean_CVD+INF_long.txt", as.is = T, check.names = F, sep = "\t")
covariates <- read.delim("../../phenotypes/phenotypes_combined_visit0.txt", sep = "\t", check.names = F, as.is = T)

ID <- gsub("_.*", "", d_wide$SampleID)
TP <- gsub(".*_", "", d_wide$SampleID)

d_wide <- cbind(ID, TP, d_wide)

pheno$Record.ID <- gsub("ID_", "", pheno$Record.ID)
colnames(pheno)[1] <- "ID"
pheno <- cbind(SampleID = paste0(pheno$ID, "_", pheno$TP), pheno)
pheno$Age <- NULL

covariates$PaNULLcovariates$Patient_id <- gsub("X","", covariates$Patient_id)
covariates <- covariates[covariates$Patient_id %in% d_wide$ID, c("Patient_id", "Age", "BMI")]
colnames(covariates)[1] <- "ID"


#
# Protein vs TP
#

lmm_res_prot_tp <- data.frame(matrix(nrow = (ncol(d_wide) -4), ncol = 2))
colnames(lmm_res_prot_tp) <- c("prot", "pval")
cnt <- 1
for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res <- lmm_prot_tp_poly3_adj_age_bmi(d_wide, prot, covariates)
    lmm_res_prot_tp[cnt,] <- c(prot, res)
    cnt <- cnt + 1
}

lmm_res_prot_tp <- na.omit(lmm_res_prot_tp) %>%
  mutate(across(-c( prot), as.numeric)) 

lmm_res_prot_tp$BH_pval <- p.adjust(lmm_res_prot_tp$pval, method = 'BH')
lmm_res_prot_tp <- lmm_res_prot_tp[order(lmm_res_prot_tp$pval),]

write.table(lmm_res_prot_tp, file = "prot_vs_tp_poly3_lmm_adj_age_bmi.txt", quote = F, sep = "\t", row.names = FALSE)
signif <- lmm_res_prot_tp[lmm_res_prot_tp$BH_pval < 0.05,]
nrow(signif)


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

pdf("../plots/radian.pdf", width = 10, height = 10)
make_radian_plot(d_wide, signif$prot)
dev.off()

#
# Cluster temporal trajectories - try on subset first
#
n_points <- 10

#all_prots <-  colnames(d_wide)[4:ncol(d_wide)]
all_prots <- sample(colnames(d_wide)[4:ncol(d_wide)], 100)

prot_dist <- data.frame(matrix(nrow = length(all_prots) * length(all_prots), ncol = 4))
colnames(prot_dist) <- c("prot1", "prot2", "eucl",  "eucl_coef")

prot_coefs <- data.frame(matrix(nrow = length(all_prots) , ncol = 4))
row.names(prot_coefs) <- all_prots

prot_coefs_raw <- data.frame(matrix(nrow = length(all_prots) , ncol = 4))
row.names(prot_coefs_raw) <- all_prots


prot_trajs <- data.frame(matrix(nrow = length(all_prots) , ncol = n_points))
row.names(prot_trajs) <- all_prots
colnames(prot_trajs) <- seq(1,4, length.out = n_points)
cnt <- 1

prot1_passed <- c()
start.time <- Sys.time()
for (prot1 in all_prots){
  gls_fit1 <- fit_gls(d_wide, prot1, n = n_points, scale = T, covariates)
  prot_coefs[prot1,] <- gls_fit1$coefficients
  prot_trajs[prot1,] <- gls_fit1$predicted
  
  gls_fit1_raw <- fit_gls(d_wide, prot1, n = n_points, scale = T, covariates, poly_raw = T)
  prot_coefs_raw[prot1,] <- gls_fit1_raw$coefficients
  
  for (prot2 in all_prots){
    if (prot2 %in% prot1_passed) next
    
    gls_fit2 <- fit_gls(d_wide, prot2, n = n_points, scale = T, covariates)
    
    # Euclidean distance between trajectories
    eucl_dist <- sum(abs(as.numeric(gls_fit1$predicted) - as.numeric(gls_fit2$predicted)))/length(gls_fit1$predicted)
    
    #Euclidean distance in the lm coefficients space
    eucl_dist_coef <- TSdist::EuclideanDistance(gls_fit1$coefficients, gls_fit2$coefficients)
    
    prot_dist[cnt,] <- c(prot1, prot2, eucl_dist, eucl_dist_coef)
    cnt <- cnt + 1
  }
  prot1_passed <- c(prot1_passed, prot1)
}
end.time <- Sys.time()

prot_dist <- na.omit(prot_dist) %>%
  mutate(across(-c( prot1, prot2), as.numeric)) 
prot_dist_wide <- my_pivot_wider(prot_dist, "prot1", "prot2", "eucl")
prot_dist_wide[lower.tri(prot_dist_wide)] <- t(prot_dist_wide)[lower.tri(prot_dist_wide)]

# prot_dist_scaled <- prot_dist %>%
#   group_by(prot1) %>%
#   mutate(eucl_scaled = scale_this(eucl), eucl_coef_scaled = scale_this(eucl_coef))
# 
# prot_dist_wide_scaled <- my_pivot_wider(prot_dist_scaled, "prot1", "prot2", "eucl_scaled")
# prot_dist_wide_scaled[lower.tri(prot_dist_wide_scaled)] <- t(prot_dist_wide_scaled)[lower.tri(prot_dist_wide_scaled)]


#convert distance to similarity metric
prot_dist$eucl_similarity <- 1 / (1 + prot_dist$eucl)

prot_sim_wide <- my_pivot_wider(prot_dist, "prot1", "prot2", "eucl_similarity")
prot_sim_wide[lower.tri(prot_sim_wide)] <- t(prot_sim_wide)[lower.tri(prot_sim_wide)]

pheatmap(prot_sim_wide)

#ggplot(prot_dist, aes(eucl, y = eucl_coef)) + geom_point()


#
# Now run on all proteins but get only orthogonal coefficients and predicted trajectories:
#

all_prots <-  colnames(d_wide)[4:ncol(d_wide)]

prot_coefs_all <- data.frame(matrix(nrow = length(all_prots) , ncol = 4))
row.names(prot_coefs_all) <- all_prots

prot_coefs_raw <- data.frame(matrix(nrow = length(all_prots) , ncol = 4))
row.names(prot_coefs_raw) <- all_prots

prot_trajs_all <- data.frame(matrix(nrow = length(all_prots) , ncol = n_points))
row.names(prot_trajs_all) <- all_prots
colnames(prot_trajs_all) <- seq(1,4, length.out = n_points)

for (prot1 in all_prots){
  gls_fit1 <- fit_gls(d_wide, prot1, n = n_points, scale = T, covariates)
  prot_coefs_all[prot1,] <- gls_fit1$coefficients
  prot_trajs_all[prot1,] <- gls_fit1$predicted
  
  gls_fit1_raw <- fit_gls(d_wide, prot1, n = n_points, scale = T, covariates, poly_raw = T)
  prot_coefs_raw[prot1,] <- gls_fit1_raw$coefficients
}

#
# LMM with TP as ordered factor - get coef
#

all_prots <-  colnames(d_wide)[4:ncol(d_wide)]

prot_coef_of <- data.frame(matrix(nrow = length(all_prots) , ncol = 3))
row.names(prot_coef_of) <- all_prots
colnames(prot_coef_of) <- c("L", "Q", "C")

for (prot1 in all_prots){
  prot_coef_of[prot1,] <- lmm_prot_tp_TP_factor_adj_age_bmi(d_wide, prot1, covariates, scale = T)
}


#
# rmcorr and lmm
#
all_prots <-  colnames(d_wide)[4:ncol(d_wide)]
prot_rmcorr <- data.frame(matrix(nrow = length(all_prots) * length(all_prots), ncol = 4))
colnames(prot_rmcorr) <- c("prot1", "prot2", "r", "p")

prot_lmm <- data.frame(matrix(nrow = length(all_prots) * length(all_prots), ncol = 4))
colnames(prot_lmm) <- c("prot1", "prot2", "r", "p")

pb = txtProgressBar(min = 0, max = length(all_prots), initial = 0) 
stepi = 0
for (prot1 in all_prots){
  #prot1_passed <- c(prot1_passed, prot1)
  setTxtProgressBar(pb,stepi)
  stepi <- stepi + 1
  #print(prot1)
  for (prot2 in all_prots){
    #if (prot2 %in% prot1_passed) next
    if (prot1 == prot2) next
    #rm_corr <- get_rmcorr(d_wide, d_wide, prot1, prot2) 
    #prot_rmcorr[cnt,] <- c(prot1, prot2, rm_corr$r, rm_corr$p)
    
    lmm_res <- lmm_pheno_prot_adj_age_bmi_v2(d_wide, d_wide, prot1, prot2, covariates, scale = T)
    prot_lmm[cnt,] <- c(prot1, prot2, lmm_res$estimate, lmm_res$pval)
    cnt <- cnt + 1
  }
}
cat("\n")
close(pb)
write.table(prot_lmm, file = "../results/prot_vs_prot_poly3_lmm_adj_age_bmi.txt", quote = F, sep = "\t", row.names = FALSE)

prot_lmm <- na.omit(prot_lmm) %>%
  mutate(across(-c( prot1, prot2), as.numeric))     
corr_matrix <- my_pivot_wider(prot_rmcorr, "prot1", "prot2", "r")
#corr_matrix[lower.tri(corr_matrix)] <- t(corr_matrix)[lower.tri(corr_matrix)]
lmm_matrix <- my_pivot_wider(prot_lmm, "prot1", "prot2", "r")
lmm_matrix <- lmm_matrix[, row.names(lmm_matrix)]
diag(lmm_matrix) <- 1
lmm_matrix[lmm_matrix > 1] <- 1
pdf("../plots/lmm_prot_prot_corrplot.pdf", height = 25, width = 25)
corrplot::corrplot(as.matrix(lmm_matrix), tl.col = 'black', order = 'hclust', addgrid.col = NA)
dev.off()
corrplot::corrplot(as.matrix(corr_matrix), tl.col = 'black', order = 'hclust')









#
# LMM protein - phenotype
#

lmm_res <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol = 4))
colnames(lmm_res) <- c("prot", "pheno", "estimate", "pval")
cnt <- 1
for (ph in colnames(pheno)[4:ncol(pheno)]) {
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res <- lmm_pheno_prot_adj_age_bmi(d_wide, pheno, prot, ph, covariates)
    lmm_res[cnt,] <- c(prot, ph, unlist(res))
    cnt <- cnt + 1
  }
}

lmm_res <- na.omit(lmm_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 

lmm_res$BH_pval <- p.adjust(lmm_res$pval)
lmm_res <- lmm_res[order(lmm_res$pval),]
