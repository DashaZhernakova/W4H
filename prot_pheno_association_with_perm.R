
library(ggplot2)
library(dplyr)
library(lme4)
library(grid)
library(gridExtra)

setwd('/mnt/sannaLAB-Temp/dasha/olink/data/')

# Read protein data
d_wide <- read.delim("olink_clean_CVD+INF.txt", as.is = T, check.names = F, sep = "\t")
d_wide <- cbind( gsub("_.*", "", d_wide$SampleID), gsub(".*_", "", d_wide$SampleID), d_wide)
colnames(d_wide)[c(1,2)] <- c("ID", "TP")

# keep only 1st visit protein measurements
d_wide_1 <- d_wide[d_wide$TP == 1,]
d_wide_1$TP <- NULL
d_wide_1$SampleID <- NULL


# Read phenotype data
pheno_0 <- read.delim("phenotypes_combined_visit0.txt", sep = "\t", check.names = F, as.is = T)
pheno_0$Patient_id <- gsub("X","", pheno_0$Patient_id)
pheno_0 <- pheno_0[pheno_0$Patient_id %in% d_wide$ID,]

# Keep a subset of phenotypes 
pheno_0 <- pheno_0[,c("Patient_id", "Age", "BMI", "Smoker", "Pregnancy_category", "Pill_use", "T2D_first_second", "caffe", "alcol")]

# Combine proteins and phenotypes
joined_data <- left_join(d_wide_1, pheno_0, by = c("ID" = "Patient_id"))


# run the association between pheno and prot using lm
regression_olink_pheno <- function(joined_data, prot, ph, scale = F){
  d <- na.omit(joined_data[,c("ID", prot, ph)])
  colnames(d) <- c("SampleID", "prot", "pheno")
  is_factor <- length(unique(d$pheno)) < 3
  d <- na.omit(d)
  if (is_factor){
    d$pheno <- as.factor(d$pheno)
  } 
  
  if(scale){
    d$prot <- scale(d$prot)
    if (! is_factor) d$pheno <- scale(d$pheno)
  }
  
  lm_fit <- lm(prot ~ pheno, data = d)
  coef <- summary(lm_fit)$coefficients
  return(list("pval" = coef[2,4], "est" = coef[2,1], "n" = nrow(d)))
}


# Association of all phenotypes vs all proteins
correl_res <- data.frame(matrix(nrow = (ncol(d_wide_1) -1) * (ncol(pheno_0) -1), ncol = 5))
colnames(correl_res) <- c("prot", "pheno", "pval", "estimate", "N")
cnt <- 1
for (ph in colnames(pheno_0)[2:ncol(pheno_0)]) {
  for (prot in colnames(d_wide_1)[2:ncol(d_wide_1)]){
    res <- regression_olink_pheno(joined_data, prot, ph, scale = T)
    correl_res[cnt,] <- c(prot, ph, unlist(res))
    cnt <- cnt + 1
  }
}

correl_res <- na.omit(correl_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 


#
# Permutations: 
#

set.seed(123) 

num_permutations <- 1000 #  number of permutations

# if we calculate the adjusted P per protein - phenotype pair, we can run the adjustment on the subset of nominally significant associations:
prots <- unique(correl_res[correl_res$pval < 0.01, "prot"])
correl_res_subs <- correl_res[correl_res$prot %in% prots,]

for (ph in colnames(pheno_0)[2:ncol(pheno_0)]) {
  print(ph)
  for (prot in prots){
    permuted_pvals <- c()
    for (i in 1:num_permutations) {
      # get the sample ids for which we have this protein and this phenotype measured:
      prot_subs <- na.omit(d_wide_1[,c("ID", prot)])
      pheno_subs <- na.omit(pheno_0[,c("Patient_id", ph)])
      shared_ids <- intersect(prot_subs$ID, pheno_subs$Patient_id)
      prot_subs <- prot_subs[prot_subs$ID %in% shared_ids,]
      permuted_pheno <- pheno_subs[pheno_subs$Patient_id %in% shared_ids,]
      # shuffle ids:
      permuted_pheno$Patient_id <- sample(permuted_pheno$Patient_id)
      # combine protein and permuted pheno data
      joined_perm <- left_join(prot_subs, permuted_pheno, by = c("ID" = "Patient_id"))
      
      # Re-run the linear models with permuted phenotypes
      res <- regression_olink_pheno(joined_perm, prot, ph, scale = T)
      permuted_pvals <- c(permuted_pvals, res[['pval']])
    }
    # get the P from the real analysis
    idx <- which(correl_res_subs$prot == prot & correl_res_subs$pheno == ph)
    obs_pval <- correl_res_subs[idx, "pval"]
    # calculate the adjusted P:
    correl_res_subs[idx, "padj_perm"] <- sum(permuted_pvals < obs_pval) / num_permutations
  }
}

nrow(correl_res_subs[correl_res_subs$padj_perm < 0.05,])
# [1] 51
