library(dplyr)
library(lme4)
library(lmerTest)
library(tidyverse)
library(mgcv)

#'
#' Models protein changes throughout the menstrual cycle using GAMs
#' 
#' @param d_wide protein dataframe with samples in rows and proteins in columns. First 3 columns should contain SampleID, ID (individual id) and TP (visit number 1-4)
#' @param prot name of the protein to model (column name in d_wide)
#' @param covariates dataframe with covariates to adjust for. Here covariates are not longitudinal and are measured at baseline. First column should be the ID (individual id).
#' @param scale binary, whether the protein values should be scaled
#' @param predict binary, whether the fitted values predicted by the model should also be calculated and returned
#' @param anova_pval binary, whether the significance of protein changes with visit should be calulated using anova comparison between model with and without timepoint instead of the standard GAM p-value
#
gam_prot_tp_adj_covar <- function(d_wide, prot, covariates, scale = F, predict = T, anova_pval = F){
  d_subs <- inner_join(d_wide[,c(prot, "ID", "TP")], covariates, by = c("ID"))
  colnames(d_subs)[1] <- "prot"
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs$ID <- as.factor(d_subs$ID)
  d_subs <- na.omit(d_subs)
  
  if (scale) d_subs$prot <- scale(d_subs$prot)
  
  fo_gam <- as.formula(paste("prot ~ s(TP, k = 4) + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+")))
  fo_gam_null <- as.formula(paste("prot ~ s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+")))
  
  model <- gam(fo_gam, data = d_subs,  method = 'REML')
  if (anova_pval){
    model0 <- gam(fo_gam_null, data = d_subs, method = 'REML')
    an <- anova.gam(model, model0)
    pval <- an$`Pr(>F)`[2]
  } else {
    pval <- summary(model)$s.table["s(TP)","p-value"]
  }
  edf <- summary(model)$s.table["s(TP)","edf"]
  fval <- summary(model)$s.table["s(TP)","F"]
  
  if (predict){
    covar_means <- as.data.frame(lapply(covariates[,-1], function(x) if(is.numeric(x)) mean(x, na.rm = TRUE) else as.factor(2)))
    
    new_data <- expand.grid(
      TP = seq(1, 4, length.out = 100),
      ID = unique(d_subs$ID),
      predicted = NA
    ) %>%
      bind_cols(
        covar_means[1,] 
      )
    predictions <- predict.gam(model, newdata = new_data,  exclude = "s(ID)", se.fit = T)
    new_data$predicted <- predictions$fit
    new_data$SE <- predictions$se.fit
    new_data$lower <- new_data$predicted - 1.96 * new_data$SE
    new_data$upper <- new_data$predicted + 1.96 * new_data$SE
    
    new_data2 <- unique(new_data[,c("TP", "predicted", "lower", "upper")])
    
    return(list(pval = pval,  edf = edf, fval = fval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID)), predicted = new_data2$predicted, lower = new_data2$lower, upper = new_data2$upper))
  } 
  return(list(pval = pval,  edf = edf, fval = fval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}

#'
#' Models protein changes throughout the menstrual cycle using LMMs
#'
#' @param d_wide protein dataframe with samples in rows and proteins in columns. First 3 columns should contain SampleID, ID (individual id) and TP (visit number 1-4)
#' @param prot name of the protein to model (column name in d_wide)
#' @param covariates dataframe with covariates to adjust for. Here covariates are not longitudinal and are measured at baseline. First column should be the ID (individual id).
#' @param scale binary, whether the protein values should be scaled
#'
lmm_prot_tp_poly3_adj_covar <- function(d_wide, prot, covariates, scale = F){
  d_subs <- inner_join(d_wide[,c(prot, "ID", "TP")], covariates, by = c("ID"))
  colnames(d_subs)[1] <- "prot"
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  fo_lmm <- as.formula(paste("prot ~ poly(TP,3) +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  model <- lmer(fo_lmm, data = d_subs)
  fo_lmm_base <- as.formula(paste("prot ~ ", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  model0 <- lmer(fo_lmm_base, data = d_subs)
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  
  return(pval)
}

#'
#' Models protein association with a phenotype using GAMs
#' 
#' @param d_wide protein dataframe with samples in rows and proteins in columns. First 3 columns should contain SampleID, ID (individual id) and TP (visit number 1-4)
#' @param pheno phenotype dataframe with samples in rows and proteins in columns. First 3 columns should contain SampleID, ID (individual id) and TP (visit number 1-4)
#' @param prot name of the protein to model (column name in d_wide)
#' @param ph name of the phenotype to use
#' @param covariates dataframe with covariates to adjust for. Here covariates are not longitudinal and are measured at baseline. First column should be the ID (individual id).
#' @param scale binary, whether the protein values should be scaled
#' @param adjust_timepoint how to adjust the association for the timepoint (visit). Can be one of the following: none - do not adjust for timepoint; linear - add timepoint as a linear predictor to the model, spline - add the timepoint as a non-linear predictor
#' @param adjust_pheno how to add the phenotype to the model. Can be one of the following: linear - add phenotype as a linear predictor to the model, spline - add the phenotype as a non-linear predictor
#' @param predict binary, whether the protein fitted values predicted by the model should also be calculated and returned
#' @param anova_pval binary, whether the significance of protein association with phenotype should be calculated using anova comparison between model with and without phenotype instead of the standard GAM p-value
#
gam_prot_pheno_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F, adjust_timepoint = 'spline', adjust_pheno = 'linear', anova_pval = F, predict = F){
  
  # Combine phenotype, protein and covariate data
  if(! "TP" %in% colnames(covariates) ){
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                         covariates, by = c("ID"))
  } else {
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                         covariates, by = c("ID", "TP"))
    covariates$TP = NULL
  }
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  d_subs$ID <- as.factor(d_subs$ID)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  # Create the model formulas
  if (adjust_timepoint == 'spline'){
    fo_gam <- paste("prot ~ s(pheno) + s(TP, k = 4) + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"))
    fo_gam_null <- paste("prot ~ s(TP, k = 4) + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"))
  } else if (adjust_timepoint == 'linear') {
    fo_gam <- paste("prot ~ s(pheno) + TP + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"))
    fo_gam_null <- paste("prot ~ TP + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"))
  } else if (adjust_timepoint == 'none') {
    fo_gam <- paste("prot ~ s(pheno) + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"))
    fo_gam_null <- paste("prot ~  s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"))
  } else {
    stop ("Wrong adjust_timepoint argument. Should be one of spline, linear or none.")
  }
  
  # Test for linear relation between protein and phenotype
  if (adjust_pheno != 'spline'){
    fo_gam <- gsub("s\\(pheno\\)", "pheno", fo_gam)
    
    model <- gam(as.formula(fo_gam), data = d_subs, method = 'REML')
    
    est <- summary(model)$p.table["pheno","Estimate"]
    se <- summary(model)$p.table["pheno","Std. Error"]
    pval <- summary(model)$p.table["pheno","Pr(>|t|)"]
    return(list(pval = pval,  est = est, se = se, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
  }
  
  # NON-linear relation between protein and phenotype
  model <- gam(as.formula(fo_gam), data = d_subs, method = 'REML')
  
  edf <- round(summary(model)$s.table["s(pheno)","edf"])
  fval <- summary(model)$s.table["s(pheno)","F"]
  
  if (anova_pval){
    model0 <- gam(as.formula(fo_gam_null), data = d_subs, method = 'REML')
    an <- anova.gam(model, model0)
    pval <- an$`Pr(>F)`[2]
  } else {
    pval <- summary(model)$s.table["s(pheno)","p-value"]
  }
  
  # Generate fitted values from the model. Use population averages, do not use the random effect
  if (predict){
    covar_means <- as.data.frame(lapply(covariates[,-1], function(x) if(is.numeric(x)) mean(x, na.rm = TRUE) else as.factor(2)))
    
    new_data <- expand.grid(
      TP = seq(1, 4),
      ID = unique(d_subs$ID),
      pheno = seq(min(d_subs$pheno), max(d_subs$pheno), length.out = 50),
      predicted = NA
    ) %>%
      bind_cols(
        covar_means %>%
          slice(1)   # to use the first row of covar_means
      )
    predictions <- predict(model, newdata = new_data,  exclude = "s(ID)", se.fit = T)
    new_data$predicted <- predictions$fit
    new_data$SE <- predictions$se.fit
    new_data$lower <- new_data$predicted - 1.96 * new_data$SE
    new_data$upper <- new_data$predicted + 1.96 * new_data$SE
    
    new_data2 <- unique(new_data[,c("TP", "pheno","predicted", "lower", "upper")])
    new_data2$TP <- as.factor(new_data2$TP)
    
    return(list(pval = pval,  edf = edf, fval = fval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID)), new_data = new_data2))
  } 
  
  return(list(pval = pval,  edf = edf, fval = fval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}



#'
#' Models protein association with a phenotype using LMMs
#' 
#' @param d_wide protein dataframe with samples in rows and proteins in columns. First 3 columns should contain SampleID, ID (individual id) and TP (visit number 1-4)
#' @param pheno phenotype dataframe with samples in rows and proteins in columns. First 3 columns should contain SampleID, ID (individual id) and TP (visit number 1-4)
#' @param prot name of the protein to model (column name in d_wide)
#' @param ph name of the phenotype to use
#' @param covariates dataframe with covariates to adjust for. Here covariates are not longitudinal and are measured at baseline. First column should be the ID (individual id).
#' @param scale binary, whether the protein values should be scaled
#' @param adjust_timepoint how to adjust the association for the timepoint (visit). Can be one of the following: none - do not adjust for timepoint; linear - add timepoint as a linear predictor to the model, cubic - add the timepoint as a 3rd degree polynomial
#
lmm_pheno_prot_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F, adjust_timepoint = "cubic"){
  
  # combine protein, phenotype and covariate data
  if(! "TP" %in% colnames(covariates) ){
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                         covariates, by = c("ID"))
  } else {
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                         covariates, by = c("ID", "TP"))
    covariates$TP = NULL
  }
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  # Generate the formula
  if (adjust_timepoint == 'cubic'){
    fo_lmm <- as.formula(paste("prot ~ poly(TP, 3) + pheno +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
    fo_lmm_base <- as.formula(paste("prot ~ poly(TP, 3) + ", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  } else if (adjust_timepoint == 'linear') {
    fo_lmm <- as.formula(paste("prot ~ TP + pheno +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
    fo_lmm_base <- as.formula(paste("prot ~ TP + ", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  } else if (adjust_timepoint == 'none') {
    fo_lmm <- as.formula(paste("prot ~ pheno +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
    fo_lmm_base <- as.formula(paste("prot ~ ", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  } else {
    stop ("Wrong adjust_timepoint argument. Should be one of cubic, linear or none.")
  }
  
  model <- lmer(fo_lmm, data = d_subs)
  model0 <- lmer(fo_lmm_base, data = d_subs)
  
  est <- summary(model)$coefficients["pheno", "Estimate"]
  se <- summary(model)$coefficients["pheno",2]
  tval <- summary(model)$coefficients["pheno",3]
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  
  return(list(estimate = est, pval = pval, se = se, tval = tval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}

