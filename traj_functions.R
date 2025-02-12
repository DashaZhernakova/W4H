
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
#my_colors <- c("#eda048", "#4d9bb0", "#075c62", "#86e6ca")
#my_colors <- c("#f5b026", "#e26413", "#9a4020", "#63552d")
#my_colors <- c("#097054", "#FFDE00", "#6599FF", "#FF9900")

library(tidyverse)
library(ggplot2)
library(rmcorr)
library(dplyr)
library(lme4)
library(splines)
library(dtw)
#library(nlme)

plot_together <- function(d_wide, pheno, prot, ph, annot = "", method = "poly3", scale = T){
  if (is.character(d_wide$TP)) d_wide$TP <- as.numeric(d_wide$TP)
  d_subs <- inner_join(d_wide[,c("ID", "TP", prot)], pheno[,c("ID", "TP", ph)], by = c("ID", "TP"))
  colnames(d_subs) <- c("ID", "TP", "prot", "pheno")

  d_subs$TP <- as.numeric(d_subs$TP)
  if (scale){
    d_subs$prot <- scale( d_subs$prot)
    d_subs$pheno <- scale( d_subs$pheno)
  }
  plot_title <- paste0(ph, " - ", prot)
  if (annot != "") plot_title <- paste0(plot_title, ", ", annot)
  if (method == "smooth"){
    g <- ggplot(d_subs, aes(x = TP, y = prot)) +
      geom_smooth(color = my_colors[2], aes(x = TP, y = pheno)) +  
      geom_smooth(color = my_colors[3], aes(x = TP, y = prot)) +  
      labs(x = "Timepoint ", y = "", 
           title = plot_title) +
      theme_minimal()
  } else if (method == "poly3"){
    g <- ggplot(d_subs, aes(x = TP, y = prot)) +
      geom_smooth(method = 'lm', formula=y ~ poly(x, 3, raw=TRUE), color = my_colors[2], aes(x = TP, y = pheno)) +  
      geom_smooth(method = 'lm', formula=y ~ poly(x, 3, raw=TRUE), color = my_colors[3], aes(x = TP, y = prot)) +  
      labs(x = "Timepoint ", y = "", 
           title = plot_title) +
      theme_minimal()
  } else if(method == 'boxplot'){
    d_subs_long <- d_subs[,-1] %>% pivot_longer(names_to = 'type', cols = c('prot', 'pheno'), values_to = 'value')
    colnames(d_subs_long)[3] <- "value"
    d_subs_long$TP <- as.factor(d_subs_long$TP)
    g <- ggplot(d_subs_long, aes(x = TP,  fill = type, y = value)) +
      geom_boxplot(position = position_dodge(width = 0.85), width = 0.8) + 
      labs(x = "Timepoint ", y = "", title = plot_title) +
      theme_minimal() +
      stat_summary(
        fun = median,
        geom = 'line',
        aes(group = type, color = type),
        position = position_dodge(width = 0.85),
        linewidth = 1
      ) +
      scale_fill_manual(values = c(my_colors[c(2,3)])) + 
      scale_color_manual(values = c('darkgoldenrod4', 'deepskyblue4'))
  }
  g
  
}

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


plot_traj_many_prots2 <- function(prot_trajs = NULL, prots, colored = T, signif = NULL, title = ""){
  #tmp <- as.data.frame(t(apply(prot_trajs[prots,], 1, scale)))

  tmp <- prot_trajs[prots,]
  colnames(tmp) <- colnames(prot_trajs)
  d_subs <- tmp %>%
    rownames_to_column(var = 'prot') %>%
    pivot_longer(cols = -prot, names_to = 'TP')
    
  
  d_subs$TP <- as.numeric(d_subs$TP)
  
  if (colored){
    g <- ggplot(d_subs, aes(x = TP, color = prot, y = value, group = prot)) +
      geom_line(stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F) +
      theme_minimal()
  } else {
    g <- ggplot(d_subs, aes(x = TP, y = value, group = prot)) +
      geom_line(stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F, alpha = 0.5) +
      theme_minimal() 
    if (!is.null(signif)){
      g <- g + geom_line(data = d_subs[d_subs$prot %in% signif,],aes(x = TP, y = value, group = prot), stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F,  color = 'red')
    }
    
    if (title != ""){
      g <- g + ggtitle(title) + theme(plot.title = element_text(size=10)) + theme_minimal()
    }
  }
  g
}

plot_traj_prots_and_pheno <- function(d_wide, pheno, prots, ph, title = ""){
  #tmp <- as.data.frame(t(apply(prot_trajs[prots,], 1, scale)))
  
  res_trajs <- data.frame(matrix(nrow = length(prots) + 1, ncol = 101))
  ph_fit <- fit_lmm_poly3_adj_covar(pheno, ph, n = 100, covariates, scale = T)
  res_trajs[1,] <- c("pheno", ph_fit$predicted)
  
  cnt <- 2
  for (prot in prots){
    fit <- fit_lmm_poly3_adj_covar(d_wide, prot, n = 100, covariates, scale = T)
    res_trajs[cnt,] <- c(prot, fit$predicted)
    cnt <- cnt + 1
  }
  colnames(res_trajs) <- c("feature", seq(1,4, length.out = 100))
  d_subs <- res_trajs %>%
    pivot_longer(cols = -feature,names_to = 'TP')
  
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs$value <- as.numeric(d_subs$value)
  d_subs$feature_type <- ifelse(d_subs$feature == 'pheno', "phenotype" ,"proteins")  
  g <- ggplot(d_subs, aes(x = TP, color = feature_type, y = value, group = feature)) +
      geom_line(stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F, linewidth = 0.5, alpha = 0.4) +
      geom_line(data = d_subs[d_subs$feature_type == 'phenotype',], color = my_colors[2], stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F, linewidth = 1) +
      theme_minimal() +
    scale_color_manual(values = my_colors[c(2,3)])
    
    if (title != ""){
      g <- g + ggtitle(title) + theme(plot.title = element_text(size=10)) + theme_minimal()
    }
  g
}

plot_medians_prots_and_pheno <- function(d_wide, pheno, prots, ph, title = "", scale = T){
  #tmp <- as.data.frame(t(apply(prot_trajs[prots,], 1, scale)))
  
  d_subs <- inner_join(pheno[,c("ID", "TP", ph)], d_wide[,c("ID", "TP", prots)], by = c("ID", "TP"))
  colnames(d_subs)[1:3] <- c("ID", "TP", "pheno")
  
  d_subs$TP <- as.factor(d_subs$TP)
  d_subs$ID <- NULL
  
  d_subs <- na.omit(d_subs)
  if (scale){
    d_subs[,-1] <- scale(d_subs[,-1])
  }
  
  medians <- aggregate(. ~ TP, data=d_subs, FUN=median) %>%
    pivot_longer(-TP, names_to = 'feature')
    
  
  medians$TP <- as.numeric(medians$TP)
  medians$value <- as.numeric(medians$value)
  medians$feature_type <- ifelse(medians$feature == 'pheno', "phenotype" ,"proteins")  
  g <- ggplot(medians, aes(x = TP, color = feature_type, y = value, group = feature)) + 
    geom_point() + 
    geom_line(linewidth = 0.5, alpha = 0.4) + 
    geom_line(data = medians[medians$feature_type == 'phenotype',], color = my_colors[2], linewidth = 1) +
    theme_minimal() +
    scale_color_manual(values = my_colors[c(2,3)])
  
  if (title != ""){
    g <- g + ggtitle(title) + theme(plot.title = element_text(size=10)) + theme_minimal()
  }
  g
}



get_rmcorr <- function(d_wide, pheno, prot, ph){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$ID <- as.factor(d_subs$ID)
  corr <- rmcorr(participant = ID, measure1 = prot, measure2 = pheno, dataset = d_subs)
  
  return(corr)
}


correlate_per_id <- function(d_wide, pheno, prot, ph){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  
  cor_res <- d_subs %>% 
    group_by(ID) %>%
    summarise(correl = cor(prot, pheno, method = 'spearman'))
  
  return(cor_res)
}


lmm_pheno_prot <- function(d_wide, pheno, prot, ph, scale = F){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  model <- lmer(pheno ~  prot + TP + (1|ID), data = d_subs)
  est <- summary(model)$coefficients["prot", "Estimate"]
  model0 <- lmer(pheno ~  TP + (1|ID), data = d_subs)
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  return(list(estimate = est, pval = pval))
}


lmm_pheno_prot_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F, adjust_timepoint = "cubic"){
  
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

lmm_pheno_prot_no_adj_covar <- function(d_wide, pheno, prot, ph,  scale = F, adjust_timepoint = "cubic"){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
                       
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  if (adjust_timepoint == 'cubic'){
    fo_lmm <- as.formula("prot ~ poly(TP, 3) + pheno + (1|ID)")
    fo_lmm_base <- as.formula("prot ~ poly(TP, 3) + (1|ID)")
  } else if (adjust_timepoint == 'linear') {
    fo_lmm <- as.formula("prot ~ TP + pheno + (1|ID)")
    fo_lmm_base <- as.formula("prot ~ TP +  (1|ID)")
  } else if (adjust_timepoint == 'none') {
    fo_lmm <- as.formula("prot ~ pheno + (1|ID)")
    fo_lmm_base <- as.formula("prot ~  (1|ID)")
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
  
  return(list(estimate = est, pval = pval, se = se, tval = tval))
}


lmm_prot_tp_interaction_pheno_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F, adjust_timepoint = "cubic"){
  d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                       covariates, by = c("ID"))
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  if (adjust_timepoint == 'cubic'){
    fo_lmm <- as.formula(paste("prot ~ poly(TP, 3) + pheno +  poly(TP, 3) * pheno +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
    fo_lmm_base <- as.formula(paste("prot ~ poly(TP, 3) + pheno +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  } else if (adjust_timepoint == 'linear') {
    fo_lmm <- as.formula(paste("prot ~ TP + pheno + TP * pheno + ", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
    fo_lmm_base <- as.formula(paste("prot ~ TP + pheno + ", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  } else {
    stop ("Wrong adjust_timepoint argument. Should be one of cubic or linear.")
  }
  
  model <- lmer(fo_lmm, data = d_subs)
  model0 <- lmer(fo_lmm_base, data = d_subs)
  
  #est <- summary(model)$coefficients["pheno", "Estimate"]
  #se <- summary(model)$coefficients["pheno",2]
  #tval <- summary(model)$coefficients["pheno",3]
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  
  return( pval)
}



lm_per_tp_pheno_prot_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F){
  d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                       covariates, by = c("ID"))
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  res_line <- c(prot, ph)
  for (tp in 1:4){
   fo_lm <- as.formula(paste("prot ~ pheno +", paste(colnames(covariates)[-1], collapse = "+")))
   model <- lm(fo_lm, data = d_subs[d_subs$TP == tp,])
   coefs <- summary(model)$coefficients
   res_line <- c(res_line, coefs['pheno', 1], -log10(coefs['pheno', 4]))
  }
  return(res_line)
}

lmm_pheno_prot_inter <- function(d_wide, pheno, prot, ph){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  model <- lmer(pheno ~ prot + TP + prot*TP + (1|ID), data = d_subs)
  #est <- summary(model)$coefficients["prot", "Estimate"]
  model0 <- lmer(pheno ~ prot + TP + (1|ID), data = d_subs)
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  return(pval)
}

lmm_prot_tp_poly3 <- function(d_wide, prot){
  d_subs <- d_wide[,c("ID", "TP", prot)]
  colnames(d_subs) <- c("ID", "TP", "prot")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  model <- lmer(prot ~ poly(TP,3) + (1|ID), data = d_subs)
  model0 <- lmer(prot ~ (1|ID), data = d_subs)
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  return (pval)
}

lmm_prot_tp_poly3_adj_age_bmi <- function(d_wide, prot, covariates, scale = F){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], covariates, by = c("ID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "Age", "BMI")
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  model <- lmer(prot ~ Age + BMI + poly(TP,3) + (1|ID), data = d_subs)
  model0 <- lmer(prot ~ Age + BMI + (1|ID), data = d_subs)
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  #est <- summary(model)$coefficients["prot", "Estimate"]
  return(pval)
}

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
get_ICC <- function(d_wide, prot){
  d_subs <- d_wide[,c(prot, "ID", "TP")]
  colnames(d_subs)[1] <- "prot"
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  d_subs$prot <- scale(d_subs$prot)
  
  m <- lmer(prot ~ 1 + (1|ID), data = d_subs)
  return (as.numeric(VarCorr( m )$ID))
}

get_AUC_difference <- function(traj1, traj2, n_points = 500) {
  x <- seq(1, 4, length.out = n_points)  # Replace with your actual x-values
  
  # Calculate the absolute difference between the two curves
  abs_diff <- abs(traj1 - traj2)
  
  # Calculate the area between the curves using the trapezoidal rule
  area <- sum((abs_diff[-1] + abs_diff[-length(abs_diff)]) / 2 * diff(x))
  
  return(area)
}

fit_lmm_poly3_adj_covar <- function(d_wide, prot, n = 10, covariates, scale = F, poly_raw = F){
  d_subs <- inner_join(d_wide[,c(prot, "SampleID", "ID", "TP")], covariates, by = c("ID"))
  colnames(d_subs)[1] <- "prot"
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  d_subs$ID <- as.factor(d_subs$ID)
  if (scale) d_subs[,"prot"] <- scale(d_subs[,"prot"])
  
  # make all columns with less than 3 unique values as factors  
  d_subs[] <- lapply(d_subs, function(col) {
    if (length(unique(col)) < 3) {
      return(factor(col))
    } else {
      return(col)
    }
  })
  
  if (poly_raw){
    fo_lmm <- as.formula(paste("prot ~ poly(TP, 3, raw = TRUE) +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  } else {
    fo_lmm <- as.formula(paste("prot ~ poly(TP,3) +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  }
  lmm_fit <- lmer(fo_lmm, data = d_subs)
  
  covar_means <- as.data.frame(lapply(covariates[,-1], function(x) if(is.numeric(x)) mean(x, na.rm = TRUE) else as.factor(2)))
  new_data <- cbind(data.frame(TP = seq(1,4, length.out = n),  predicted = NA), covar_means)
  new_data$predicted <- predict(lmm_fit, newdata = new_data, re.form = NA)
  
  coef <- summary(lmm_fit)$coefficients[grepl("poly\\(TP", row.names(summary(lmm_fit)$coefficients)),1]
  return(list("predicted" = new_data$predicted, "coefficients" = coef))
}

lmm_prot_tp_TP_factor_adj_age_bmi <- function(d_wide, prot, covariates, scale = F){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], covariates, by = c("ID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "Age", "BMI")
  
  d_subs$TP_of <- factor(d_subs$TP, ordered = T)
  d_subs <- na.omit(d_subs)
  if(scale) d_subs$prot <- scale(d_subs$prot)
  
  model <- lmer(prot ~ Age + BMI +  TP_of + (1|ID), data = d_subs)
  coef <- summary(model)$coefficients[c("TP_of.L", "TP_of.Q", "TP_of.C"),1]
  
  return(coef)
}

# scatter colored by visit to see the relationship at each visit
scatter_col_tp <- function(d_wide, pheno, prot, ph, scale = F, add_points = F){
  if ("SampleID" %in% colnames(d_wide)){
    d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
    colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  } else {
    d_subs <- inner_join(d_wide[,c("ID", "TP", prot)], pheno[,c("ID", "TP" ,ph)], by = c("ID", "TP"))
    colnames(d_subs) <- c("ID", "TP", "prot", "pheno")
  }
  d_subs$TP <- as.factor(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if(scale){
    d_subs$pheno <- scale(d_subs$pheno)
    d_subs$prot <- scale(d_subs$prot)
  }
  g <- ggplot(d_subs, aes(x = prot, y = pheno, colour = TP)) + 
    geom_smooth(method = 'lm', alpha = 0.2) + 
    stat_smooth(method = 'lm', se = F) +
    theme_minimal() +
    labs(x = prot, y = ph, 
         title = paste0(ph, " - ", prot))  +
    scale_color_manual(values = my_colors)
  
  if (add_points) g <- g + geom_point()
  g
}


fit_poly3 <- function(d_wide, prot, n = 7, scale = T, covariates = NULL){
  if(is.null(covariates)){
    d_subs <- d_wide[,c("ID", "TP", prot)]
    colnames(d_subs) <- c("ID", "TP", "prot")
    d_subs$TP <- as.numeric(d_subs$TP)
    d_subs <- na.omit(d_subs)
    if (scale) d_subs$prot <- scale(d_subs$prot)
  
    lm_fit <- lm(prot ~ poly(TP,3), data = d_subs)
    new_data <- data.frame(TP = seq(1,4, length.out = n), predicted = NA)
  } else {
    d_subs <- inner_join(d_wide[,c("ID", "TP", prot)], covariates, by = c("ID"))
    colnames(d_subs)[3] <- "prot"
    d_subs$TP <- as.numeric(d_subs$TP)
    d_subs <- na.omit(d_subs)
    if (scale) d_subs$prot <- scale(d_subs$prot)
    
    lm_fit <- lm(prot ~ Age + BMI + poly(TP,3), data = d_subs)
    new_data <- data.frame(TP = seq(1,4, length.out = n), Age = mean(d_subs$Age), BMI = mean(d_subs$BMI), predicted = NA)
  }
  new_data$predicted <- predict(lm_fit, newdata = new_data )
  coef <- summary(lm_fit)$coefficients[c("poly(TP, 3)1", "poly(TP, 3)2", "poly(TP, 3)3")]
  return(list("predicted" = new_data$predicted, "coefficients" = coef))
}


fit_poly3_get_coef <- function(d_wide, prot, scale = T){
  d_subs <- d_wide[,c("ID", "TP", prot)]
  colnames(d_subs) <- c("ID", "TP", "prot")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) d_subs$prot <- scale(d_subs$prot)
  
  lm_fit <- lm(prot ~ poly(TP,3), data = d_subs)
  
  return (summary(lm_fit)$coefficients[-1,1])
}

get_mean_per_tp <- function(d_wide, prot){
  d_subs <- d_wide[,c("ID", "TP", prot)]
  colnames(d_subs) <- c("ID", "TP", "prot")
  mean_prot_by_TP <- d_subs %>%
    group_by(TP) %>%
    summarize(mean_prot = mean(prot, na.rm = TRUE))
  return(mean_prot_by_TP)
}

fit_splines <- function(d_wide, prot, n = 50, scale = F){
  d_subs <- d_wide[,c("ID", "TP", prot)]
  colnames(d_subs) <- c("ID", "TP", "prot")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  spline_fit <- lm(prot ~ ns(TP, df = 3), data = d_subs)
  new_data <- data.frame(TP = seq(1,4, length.out = n), predicted = NA)
  new_data$predicted <- predict(spline_fit, newdata = new_data )
  if (scale) new_data$predicted <- scale(new_data$predicted)
  return (new_data)
}

get_dtw <- function(d_wide, pheno, prot, ph, fitting = "poly3"){
  if (fitting == "poly3"){
    traj_prot <- fit_poly3(d_wide, prot, scale = T)
    traj_ph <- fit_poly3(pheno, ph, scale = T)
  } else {
    traj_prot <- fit_splines(d_wide, prot, scale = T)
    traj_ph <- fit_splines(pheno, ph, scale = T)
  }
  dtw_dist <- dtw(traj_prot$predicted, traj_ph$predicted)$distance
  return(dtw_dist)
}

get_auc <- function(d_wide, prot, fitting = "poly3") {
  if (fitting == "poly3"){
    fit <- fit_poly3(d_wide, prot, scale = T)
  } else {
    fit <- fit_splines(d_wide, prot, scale = T)
  }
  x = fit$TP
  y = fit$predicted
  return (sum(diff(x) * (head(y,-1)+tail(y,-1)))/2)
}


make_radian_plot <- function(d_wide, prots, value = 'mean'){
  library(ggradar)
  if (value == 'mean'){
    d_subs <- d_long[d_long$Assay %in% prots,]
    mean_prot_by_TP <- d_subs %>%
      group_by(Timepoint, Assay) %>%
      summarize(value = 1 + mean(NPX, na.rm = TRUE))
  }
  
  data_wide <- mean_prot_by_TP %>%
    pivot_wider(names_from = Assay, values_from = value)
  
  data_wide$Timepoint <- as.factor(data_wide$Timepoint)

  g <- ggradar(
    data_wide,
    group.colours = my_colors,
    legend.title = "Timepoints",
    axis.label.size = 2,
    grid.label.size = 4,
    group.line.width = 1,
    group.point.size = 0,
    background.circle.colour = "white",
    gridline.mid.colour = "gray"
  )
  g
}



gls_pheno_prot <- function(d_wide, pheno, prot, ph){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  # Unstructured correlation: different correlations for every pair of timepoints
  gls_fit <- gls(pheno ~ prot, data = d_subs, correlation = corSymm(form = ~ TP | ID))
  # Autoregressive correlation: measurements taken closer in time are more highly correlated, and the correlation decreases exponentially with increasing time intervals
  #gls_fit <- gls(pheno ~ prot, data = d_subs, correlation = corAR1(form = ~ TP | ID))
  #Compound symmetry: corCompSymm
  
  pval <- summary(gls_fit)$tTable["prot", "p-value"]
  est <- summary(gls_fit)$tTable["prot", "Value"]
  return(list(estimate = est, pval = pval))
}



fit_gls <- function(d_wide, prot, n = 8, scale = F, covariates = NULL, poly_raw = FALSE){
  if (is.null(covariates)) {
    d_subs <- d_wide[,c("ID", "TP", prot)]
    colnames(d_subs) <- c("ID", "TP", "prot")
    d_subs$TP <- as.numeric(d_subs$TP)
    d_subs <- na.omit(d_subs)
    
    if (scale) d_subs$prot <- scale(d_subs$prot)
    
    # Unstructured correlation: different correlations for every pair of timepoints
    if (poly_raw) {
      gls_fit <- gls(prot ~ poly(TP, 3, raw = TRUE), data = d_subs, correlation = corSymm(form = ~ TP | ID))
    } else {
      gls_fit <- gls(prot ~ poly(TP, 3), data = d_subs, correlation = corSymm(form = ~ TP | ID))
    }  
      
      # Autoregressive correlation: measurements taken closer in time are more highly correlated, and the correlation decreases exponentially with increasing time intervals
    #gls_fit <- gls(prot ~ poly(TP,3), data = d_subs, correlation = corAR1(form = ~ TP | ID))
    new_data <- data.frame(TP = seq(1,4, length.out = n), predicted = NA)
  } else {
    d_subs <- inner_join(d_wide[,c("ID", "TP", prot)], covariates, by = c("ID"))
    colnames(d_subs)[3] <- "prot"
    d_subs$TP <- as.numeric(d_subs$TP)
    d_subs <- na.omit(d_subs)
    
    if (scale) d_subs$prot <- scale(d_subs$prot)
    if (poly_raw) {
      gls_fit <- gls(prot ~ Age + BMI + poly(TP, 3, raw = TRUE), data = d_subs, correlation = corSymm(form = ~ TP | ID))
    } else {
      gls_fit <- gls(prot ~ Age + BMI + poly(TP, 3), data = d_subs, correlation = corSymm(form = ~ TP | ID))
    }
    new_data <- data.frame(TP = seq(1,4, length.out = n), Age = mean(d_subs$Age), BMI = mean(d_subs$BMI), predicted = NA)
  }
  
  new_data$predicted <- predict(gls_fit, newdata = new_data )
  coef <- c(gls_fit$coefficients[1],tail(gls_fit$coefficients,3))
  return(list("predicted" = new_data$predicted, "coefficients" = coef))
}

my_pivot_wider <- function(d, row_names, names_from, values_from){
  d2 <- d[,c(row_names, names_from, values_from)] %>%
    pivot_wider(names_from = {{names_from}}, values_from = {{values_from}})
  d2 <- as.data.frame(d2)
  row.names(d2) <- d2[,row_names]
  d2[,row_names ] <- NULL
  return(d2)
}


# First derivative as a function of x
first_derivative <- function(coefs, x) {
  return(as.numeric(coefs[2] + 2 * coefs[3] * x + 3 * coefs[4] * x^2))
}

# Second derivative as a function of x
second_derivative <- function(coefs, x) {
  return(as.numeric(2 * coefs[3] + 6 * coefs[4] * x))
}

get_derivatives <- function(coefs, n_points = 10){
  deriv1 <- c()
  deriv2 <- c()
  for (x in seq(1,4, length.out = n_points)){
    deriv1 <- c(deriv1, first_derivative(coefs, x))
    deriv2 <- c(deriv2, second_derivative(coefs, x))
  }
  return(c(deriv1, deriv2))
  #return(list("deriv1" = deriv1, "deriv2" = deriv2))
}

get_inflection_points <- function(coefs){
  beta2 <- coefs[3]  # Coefficient for x^2
  beta3 <- coefs[4]  # Coefficient for x^3

  # Solve for the x value of the inflection point
  x_inflection <- -beta2 / (3 * beta3)
  x_inflection
}



fit_function <- function(x) {
   coefs[1] + coefs[2] * x + coefs[3] * x^2 + coefs[4] * x^3
}


get_roots_inflation_critical_points <- function(coefs){
  # 1. Find the roots of the polynomial
  roots <- polyroot(coefs)
  cat("Roots of the polynomial:", roots, "\n")
  
  a <- coefs[4]
  b = coefs[3]
  c = coefs[2]
  d = coefs[1]
  # 2. Calculate the inflection point
  inflection_point <- -b / (3 * a)
  cat("Inflection point (x-value):", inflection_point, "\n")
  
  # 3. Calculate critical points (where first derivative is zero)
  # The first derivative is: f'(x) = 3ax^2 + 2bx + c
  # Solving 3ax^2 + 2bx + c = 0
  
  # Calculate the discriminant
  discriminant <- (2 * b)^2 - 4 * 3 * a * c
  
  if (discriminant >= 0) {
    # Real solutions
    critical_points <- c((-2 * b + sqrt(discriminant)) / (6 * a),
                         (-2 * b - sqrt(discriminant)) / (6 * a))
    cat("Critical points (real x-values where f'(x) = 0):", critical_points, "\n")
  } else {
    # Complex solutions
    critical_points <- c((-2 * b + sqrt(as.complex(discriminant))) / (6 * a),
                         (-2 * b - sqrt(as.complex(discriminant))) / (6 * a))
    cat("Critical points (complex x-values where f'(x) = 0):", critical_points, "\n")
  }
}

get_simplified_coefs <- function(d_wide, prot, scale = T){
  
  d_subs <- d_wide[,c("ID", "TP", prot)]
  colnames(d_subs) <- c("ID", "TP", "prot")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) d_subs$prot <- scale(d_subs$prot)
  
  gls_fit <- gls(prot ~ poly(TP, 3), data = d_subs, correlation = corSymm(form = ~ TP | ID))
  coefs <- gls_fit$coefficients
  
  coefs <- coefs[c(2,3,4)]
  coefs[coefs > 2] <- 3
  coefs[coefs < -2] <- -3
  
  coefs[coefs > 1 & coefs < 2] <- 2
  coefs[coefs < -1 & coefs > -2] <- -2
  
  coefs[coefs < 1 & coefs > 0 ] <- 1
  coefs[coefs > -1 & coefs < 0 ] <- -1
  
  return(coefs)
}

plot_clusters <- function(cl, method = "", num_k = "", colored = F, signif = NULL, save_pdf = T, add_cluster_name = F, out_path = NA){
  plot_list = list()
  for (cluster in unique(cl)){
    title = ifelse(add_cluster_name, cluster, "")
    plot_list[[cluster]] <- plot_traj_many_prots2(prot_trajs, names(cl[cl == cluster]), colored = colored, signif, title)
  }
  
  pdf_path = ifelse(is.na(out_path), 
                    paste0("../plots/clustering_signif_v2/", method, "_k", num_k, ".pdf"),
                    out_path)
  
  
  ncols = 4
  nrows = ceiling(length(unique(cl))/4)
  
  if(save_pdf) pdf(pdf_path, width = 4*ncols, height = 4*nrows)
  grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)
  
  
  #if (length(unique(cl)) < 10){
  #  grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)  
  #} else if(length(unique(cl)) < 10) {
  #  grid.arrange(grobs = plot_list, ncol = 4, nrow = 4)  
  #} else {
  #  grid.arrange(grobs = plot_list, ncol = 5, nrow = 5)  
  #}
  if(save_pdf) dev.off()
  
}


regress_covariates_per_tp <- function(data, covariates){
  data_tp <- as.data.frame(subset(data, select = -c(TP, ID)))
  if ("SampleID" %in% colnames(data_tp)) data_tp <- subset(data_tp, select = -c(SampleID))
  row.names(data_tp) <- data$ID
  data_tp <- na.omit(data_tp)
  
  covariates <- na.omit(covariates)
  row.names(covariates) <- covariates$ID
  if (colnames(covariates)[1] == 'ID' && ncol(covariates) == 2) {
    covariates <- covariates[,-1,drop = F] 
  } else {
    covariates$ID <- NULL
  }

  ids <- intersect(row.names(covariates), row.names(data_tp))
  data_tp <- data_tp[ids,]
  covariates <- covariates[ids,,drop = F]
  
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


regress_covariates_lmm <- function(data, covar_data){
  
  if (!"SampleID" %in% colnames(covar_data)) {
    covar_data <- cbind(paste0(covar_data$ID, "_",covar_data$TP), covar_data)
    colnames(covar_data)[1] <- "SampleID"
  }
  
  d_adj <- data[,c("SampleID", "ID", "TP")]
  
  data[,"TP"] <- NULL
  covar_data[,c("ID", "TP")] <- NULL
  
  cnt <- 1
  for (ph in colnames(data)[3: (ncol(data))]){
    subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "SampleID"))
    colnames(subs)[3] <- 'pheno'
    
    fo_lmm <- as.formula(paste("pheno ~ ", paste(colnames(covar_data)[-1], collapse = "+"), "+ (1|ID)"))
    
    
    lmm_fit <- lmer(fo_lmm, data = subs)
    subs[,ph] <- subs$pheno - predict(lmm_fit, re.form = NA)
    
    d_adj <- left_join(d_adj, subs[, c("SampleID", ph)], by = "SampleID")
  }
  return(d_adj)
}
