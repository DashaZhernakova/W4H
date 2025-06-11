setwd("/Users/Dasha/work/Sardinia/W4H/olink/")
source("scripts/utility_functions_for_Franceen.R")

library(dplyr)

out_basedir <- "results/pheno_batch2_prot_rm_outliers_4sd/"
d_wide <- read.delim("data/olink_clean_CVD+INF_rm_outliers_4sd.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))
covariates <- read.delim("data/covariates_age_bmi_storage_preg.txt", sep = "\t", check.names = F, as.is = T, colClasses = c(ID = "character"))
pheno <- read.delim("../phenotypes/blood_pheno_13122024_log_adj_batch_storage_rm_outiers_4sds.txt.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))

d_wide$TP <- as.numeric(d_wide$TP)

# Convert binary covariates into factors
covariates[] <- lapply(covariates, function(col) {
  if (length(unique(col)) < 3) {
    return(factor(col))
  } else {
    return(col)
  }
})
covariate_names <- c("Age","BMI","storage_months", "Pregnancy_category")

################################################################################
# Protein vs TP GAM and LMM
################################################################################

gam_res_prot_tp <- data.frame(matrix(nrow = (ncol(d_wide) -4), ncol = 7))
colnames(gam_res_prot_tp) <- c("prot", "gam_pval", "gam_edf", "gam_fval", "n", "n_samples", "lmm_pval")

cnt <- 1
for (prot in colnames(d_wide)[4:ncol(d_wide)]){
  # GAM modelling of protein trajsectories
  res_gam <- gam_prot_tp_adj_covar(d_wide, prot, covariates, scale = T, predict = F)
  # LMM modelling of protein trajsectories
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
# GAM protein vs phenotype 
################################################################################

all_pheno <- colnames(pheno)[4:ncol(pheno)]

lmm_res <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol =14))
colnames(lmm_res) <- c("prot", "pheno", 
                       paste( c("estimate", "pval", "se", "tval", "N", "N_unique"), "withTP", sep = "_"),
                       paste( c("estimate", "pval", "se", "tval", "N", "N_unique"), "noTP", sep = "_")
)

gam_res <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol = 17))
colnames(gam_res) <- c("prot", "pheno", 
                       paste(c("pval", "estimate", "se", "n", " n_samples"), "spline", sep =  "_"),
                       paste(c("pval", "estimate", "se", "n", " n_samples"), "linear", sep =  "_"),
                       paste(c("pval", "estimate", "se", "n", " n_samples"), "noTP", sep =  "_")
)

cnt <- 1
for (ph in all_pheno) {
  cat(ph, "\n")
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    # run GAM association between protein and phenotype adjusting for the visit in a NON-LINEAR way:
    res_gam <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'spline', anova_pval = F)
    # run GAM association between protein and phenotype adjusting for the visit in a LINEAR way:
    res_gam_lin <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'linear', anova_pval = F)
    # run GAM association between protein and phenotype with NO adjustment for the visit:
    res_gam_noTP <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'none', anova_pval = F)
    
    gam_res[cnt,] <- c(prot, ph, unlist(res_gam), unlist(res_gam_lin), unlist(res_gam_noTP))
    
    # run LMM association between protein and phenotype adjusting for the visit in a NON-LINEAR (cubic) way:
    res_lmm <- lmm_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'cubic')
    # run LMM association between protein and phenotype with NO adjustment for the visit:
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
