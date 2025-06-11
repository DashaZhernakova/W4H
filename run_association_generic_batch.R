args <- commandArgs(trailingOnly = TRUE)
library(tidyverse)
library(dplyr)

source("/Users/Dasha/work/Sardinia/W4H/olink//scripts/utility_functions_generic_association_analysis.R")

fname1 <- args[1]
fname2 <- args[2]
covar_fname <- args[3]
prot_list_fname <- args[4]
out_fname <- args[5]

run_lmm <- F

d1 <- read.delim(fname1, as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character", TP = "numeric"))
d2 <- read.delim(fname2, as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character", TP = "numeric"))

prot_list <- as.character(read.delim(prot_list_fname, as.is = T, check.names = F, sep = "\t", header = F))

covariates <- read.delim(covar_fname, as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character", TP = "numeric"))

cat("Running association analysis for ", length(prot_list), " proteins\n")
cat("File 1: ", fname1, "\n")
cat("File 2: ", fname2, "\n")
cat ("Covariates: ", colnames(covariates)[2:ncol(covariates)], "\n")

gam_res <- data.frame(matrix(nrow = length(prot_list) * (ncol(d2) -4), ncol = 7))
colnames(gam_res) <- c("feature1", "feature2", "pval", "edf_round", "fval", "n", " n_samples")

if (run_lmm) {
  lmm_res <- data.frame(matrix(nrow = length(prot_list) * (ncol(d2) -4), ncol = 8))
  colnames(lmm_res) <- c("feature1", "feature2", "estimate", "pval", "se", "tval", "N", "N_unique")
}

cnt <- 1
for (prot in prot_list) {
  cat(prot, "\n")
  for (ph in colnames(d2)[4:ncol(d2)]){
    res_gam <- gam_prot_pheno_adj_covar(d1, d2, prot, ph, covariates, scale = T, adjust_timepoint = 'spline', anova_pval = F)
    gam_res[cnt,] <- c(prot, ph, unlist(res_gam))
    
    if (run_lmm) {
      res_lmm <- lmm_pheno_prot_adj_covar(d1, d2, prot, ph, covariates, scale = T, adjust_timepoint = 'cubic')
      lmm_res[cnt,] <- c(prot, ph, unlist(res_lmm))
    }
    cnt <- cnt + 1
  }
}

gam_res <- na.omit(gam_res) %>%
  mutate(across(-c(feature1, feature2), as.numeric)) 
gam_res <- gam_res[order(gam_res$pval),]

write.table(gam_res, file = paste0(out_fname, ".GAM_results.txt"), quote = F, sep = "\t", row.names = FALSE)

if (run_lmm) {
  lmm_res <- na.omit(lmm_res) %>%
    mutate(across(-c(feature1, feature2), as.numeric)) 
  lmm_res <- lmm_res[order(lmm_res$pval),]
  write.table(lmm_res, file = paste0(out_fname, ".LMM_results.txt"), quote = F, sep = "\t", row.names = FALSE)
}

cat("Finished!\n")
