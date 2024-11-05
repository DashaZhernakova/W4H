
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

d_wide$ID <- gsub("_.*", "", d_wide$SampleID)
d_wide$TP <- gsub(".*_", "", d_wide$SampleID)

pheno$Record.ID <- gsub("ID_", "", pheno$Record.ID)
colnames(pheno)[1] <- "ID"
pheno <- cbind(SampleID = paste0(pheno$ID, "_", pheno$TP), pheno)

d_long <- read.delim("olink_clean_CVD+INF_long.txt", as.is = T, check.names = F, sep = "\t")


#
# run all distances
#

res_cor <- data.frame(matrix(ncol = 12, nrow = (ncol(pheno) - 4) * ncol(d_wide)))
cnt <- 1
for (ph in colnames(pheno)[5:ncol(pheno)]){
  print(ph)
  for (prot in colnames(d_wide)[2:(ncol(d_wide) - 2)]){
    
    # rm corr
    rm_corr <- get_rmcorr(d_wide, pheno, prot, ph)
    
    # LMM
    lmm_res <- lmm_pheno_prot(d_wide, pheno, prot, ph)
    
    # correlation of poly3 lm coefficients
    coef1 <- fit_poly3_get_coef(d_wide, prot)
    coef2 <- fit_poly3_get_coef(pheno, ph)
    poly3_corr <- cor(coef1, coef2)
    
    # DTW distance
    dtw <- get_dtw(d_wide, pheno, prot, ph)
    
    # AUC difference
    auc1 <- get_auc(d_wide, prot)
    auc2 <- get_auc(pheno, ph)
    auc_dif <- abs(auc1 - auc2)
    
    # Euclidean distance between trajectories
    traj1 <- fit_poly3(d_wide, prot, scale = T)
    traj2 <- fit_poly3(pheno, ph, scale = T)
    eucl_dist <- sum(abs(as.numeric(traj1$predicted) - as.numeric(traj2$predicted)))/nrow(traj1)
    
    # Euclidean distance in the lm coefficients space
    eucl_dist_coef <- TSdist::EuclideanDistance(coef1, coef2)
    
    
    # Frechet distance
    fre_dist <- distFrechet(traj1$TP, as.numeric(traj1$predicted), traj2$TP, as.numeric(traj2$predicted))
    
    res_cor[cnt, ] <- c(ph, prot, rm_corr$r, rm_corr$p, lmm_res$estimate, lmm_res$pval,
                        poly3_corr, dtw, auc_dif, eucl_dist, eucl_dist_coef, fre_dist)
    cnt <- cnt + 1
  }
}
colnames(res_cor) <- c("pheno", "prot", "rmcorr", "rmcorr_pval", "lmm_est", "lmm_pval",
                       "poly3_corr", "dtw", "auc_dif", "eucl", "eucl_coef", "frechet")

res_cor2 <- na.omit(res_cor) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 

dist_met = "rmcorr"
res_cor_wide <- res_cor2[,c("prot", "pheno", dist_met)] %>% 
  pivot_wider(names_from = pheno, values_from = !!sym(dist_met)) %>%
  column_to_rownames(var = "prot")

pheatmap(as.matrix(res_cor_wide), fontsize_row  = 3)

dist_mets <- c("rmcorr", "lmm_est", 
"poly3_corr", "dtw", "auc_dif", "eucl", "eucl_coef", "frechet")

for (dist_met in dist_mets){
pdf(paste0("../plots/distances/cmp_distances.", dist_met, ".smooth_poly.v2.pdf"), width = 10, height = 15)
plot_list <- list()
cnt <- 1
for (ph in unique(res_cor2$pheno)){
  if(dist_met == 'lmm_est') {
    tmp <- res_cor2[res_cor2$lmm_pval < 0.001 & res_cor2$pheno == ph,]
  } else {
    tmp <- res_cor2[ res_cor2$pheno == ph,]
  }
  for (pr in tmp[order(-tmp[,dist_met]), ][1:min(nrow(tmp),4), "prot"]){
    plot_list[[cnt]] <- plot_together(d_wide, pheno, pr, ph, formatC(tmp[tmp$prot == pr, dist_met], digits = 3), method = "poly3")
    cnt <- cnt + 1
  }
  for (pr in tmp[order(tmp[,dist_met]), ][1:min(nrow(tmp),4), "prot"]){
    plot_list[[cnt]] <- plot_together(d_wide, pheno, pr, ph, formatC(tmp[tmp$prot == pr, dist_met], digits = 3), method = "poly3")
    cnt <- cnt + 1
  }

  if (cnt == 25){
    cat(ph, "\n")
    grid.arrange(grobs = plot_list, ncol = 4, nrow = 6)  
    plot_list <- list()
    cnt <- 1
  }
}
dev.off()
}



# Interactions
res_inter <- data.frame(matrix(ncol = 3, nrow = (ncol(pheno) - 4) * ncol(d_wide)))
cnt <- 1
for (ph in colnames(pheno)[5:ncol(pheno)]){
  print(ph)
  for (prot in colnames(d_wide)[2:(ncol(d_wide) - 2)]){
    pval <- lmm_pheno_prot_inter(d_wide, pheno, prot, ph)
    res_inter[cnt,] <- c(ph, prot, pval)
    cnt <- cnt + 1
  }
}
    



#
# Cmp euclidean distances: 17BES puzzle
#
clusters <- data.frame(matrix(c('PROG', 1, "FS", 1, "LH", 1, 'FSH', 1, 'PRL', 1,
                         "X17BES", 2, "HOMA_IR", 2, "HOMA_B", 2, "INS", 2, "ALT", 2, 'AST', 2, 'Glucose', 2, 'HDLC',2 ,
                         'LDLC', 3,  'TotChol', 3, 'Trigl', 3), byrow = T, ncol = 2))
colnames(clusters) <- c('pheno', 'cluster')
res_cor3 <- left_join(res_cor2, clusters, by = "pheno")

ggplot(res_cor3, aes(x = eucl, y = eucl_coef, color = cluster, group = cluster)) + 
  geom_point(alpha = 0.3) +
  geom_smooth(method = 'lm') +
  theme_bw()

pheno_long <- pheno %>%
  pivot_longer(cols = unique(res_cor2$pheno),names_to = 'pheno', values_to = 'value') %>%
  group_by(pheno) %>%
  mutate(val_scaled = scale_this(value))
pheno_long <- left_join(pheno_long, clusters, by = "pheno")

pheno_long$val_scaled2 <- scale(pheno_long$value)
ggplot(pheno_long, aes(x = TP, y = val_scaled, group = pheno, color = cluster)) + geom_smooth(se = F) + theme_bw()


scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)

}


# to make the phenotypes comparable:
# scale the distances per phenotype

res_cor4 <- res_cor2 %>%
  group_by(pheno) %>%
  mutate(eucl_scaled = scale_this(eucl), eucl_coef_scaled = scale_this(eucl_coef)) 

ggplot(res_cor4, aes(x = eucl_scaled, y = eucl_coef_scaled, color = pheno)) + geom_point() + 
  theme_bw() + 
  xlab("Euclidean distance between curves scaled per phenotype") +
  ylab("Euclidean distance in the lm coefficient space scaled per phenotype")
