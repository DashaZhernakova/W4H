
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
library(tidyverse)
library(ggplot2)
library(rmcorr)
library(dplyr)
library(lme4)
library(splines)
library(dtw)

plot_together <- function(d_wide, pheno, prot, ph, annot = "", method = "poly3", scale = T){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  if (scale){
    d_subs$prot <- scale( d_subs$prot)
    d_subs$pheno <- scale( d_subs$pheno)
  }
  plot_title <- paste0(ph, " - ", prot)
  if (annot != "") plot_title <- paste0(plot_title, ", ", annot)
  if (method == "smooth"){
    g <- ggplot(d_subs, aes(x = TP, y = prot)) +
      geom_smooth(color = "red", aes(x = TP, y = pheno)) +  
      geom_smooth(color = "blue", aes(x = TP, y = prot)) +  
      labs(x = "Timepoint ", y = "", 
           title = plot_title) +
      theme_minimal()
  } else if (method == "poly3"){
    g <- ggplot(d_subs, aes(x = TP, y = prot)) +
      geom_smooth(method = 'lm', formula=y ~ poly(x, 3, raw=TRUE), color = "red", aes(x = TP, y =pheno)) +  
      geom_smooth(method = 'lm', formula=y ~ poly(x, 3, raw=TRUE), color = "blue", aes(x = TP, y = prot)) +  
      labs(x = "Timepoint ", y = "", 
           title = plot_title) +
      theme_minimal()
  }
  g
  
}

plot_traj_many_prots <- function(d_wide, prots, colored = T){
  d_subs <- d_wide[,c("ID","TP", prots)] %>%
    pivot_longer(cols = {prots}) %>%
    group_by(name) %>%
    mutate(value = scale(value))
  
  d_subs$TP <- as.numeric(d_subs$TP)
  if (colored){
    g <- ggplot(d_subs, aes(x = TP, color = name, y = value, group = name)) +
      geom_smooth(method = 'lm', formula=y ~ poly(x, 3, raw=TRUE), se = F) +
      theme_bw()
  } else {
    g <- ggplot(d_subs, aes(x = TP, y = value, group = name)) +
      geom_smooth(method = 'lm', formula=y ~ poly(x, 3, raw=TRUE), se = F) +
      theme_bw()
  }
    g
}

plot_traj_many_prots2 <- function(prot_trajs, prots, colored = T, signif = NULL){
  #tmp <- as.data.frame(t(apply(prot_trajs[prots,], 1, scale)))
  tmp <- prot_trajs[prots,]
  colnames(tmp) <- colnames(prot_trajs)
  d_subs <- tmp %>%
    rownames_to_column(var = 'prot') %>%
    pivot_longer(cols = -prot, names_to = 'TP')
    
  
  d_subs$TP <- as.numeric(d_subs$TP)
  if (colored){
    g <- ggplot(d_subs, aes(x = TP, color = prot, y = value, group = prot)) +
      geom_line(stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F)+
      theme_bw() 
  } else {
    g <- ggplot(d_subs, aes(x = TP, y = value, group = prot)) +
      geom_line(stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F, alpha = 0.5) +
      theme_bw()
    if (!is.null(signif)){
      g <- g + geom_line(data = d_subs[d_subs$prot %in% signif,],aes(x = TP, y = value, group = prot), stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F,  color = 'red')
    }
  }
  g
}

get_rmcorr <- function(d_wide, pheno, prot, ph){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  
  corr <- rmcorr(participant = ID, measure1 = prot, measure2 = pheno, dataset = d_subs)
  
  return(corr)
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
  ### ??
  #model <- lmer(pheno ~  (prot|TP) + (1|ID), data = d_subs)
  model <- lmer(pheno ~  prot + TP + (1|ID), data = d_subs)
  est <- summary(model)$coefficients["prot", "Estimate"]
  model0 <- lmer(pheno ~  TP + (1|ID), data = d_subs)
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  return(list(estimate = est, pval = pval))
}

lmm_pheno_prot_adj_age_bmi <- function(d_wide, pheno, prot, ph, covariates, scale = F){
  d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                       covariates, by = c("ID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno", "Age", "BMI")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  ### ??
  #model <- lmer(pheno ~  (prot|TP) + (1|ID), data = d_subs)
  model <- lmer(pheno ~  prot + TP + Age + BMI + (1|ID), data = d_subs)
  est <- summary(model)$coefficients["prot", "Estimate"]
  model0 <- lmer(pheno ~  TP +  Age + BMI +(1|ID), data = d_subs)
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  return(list(estimate = est, pval = pval))
}

# No correction for TP!!!
lmm_pheno_prot_adj_age_bmi_v2 <- function(d_wide, pheno, prot, ph, covariates, scale = F){
  d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                       covariates, by = c("ID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno", "Age", "BMI")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  ### ??
  #model <- lmer(pheno ~  (prot|TP) + (1|ID), data = d_subs)
  model <- lmer(pheno ~  prot + TP + Age + BMI + (1|ID), data = d_subs)
  est <- summary(model)$coefficients["prot", "Estimate"]
  model0 <- lmer(pheno ~  TP +  Age + BMI +(1|ID), data = d_subs)
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  return(list(estimate = est, pval = pval))
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
scatter_col_tp <- function(d_wide, pheno, prot, ph){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.factor(d_subs$TP)
  d_subs <- na.omit(d_subs)
  g <- ggplot(d_subs, aes(x = prot, y = pheno, colour = TP)) + 
    geom_smooth(method = 'lm', alpha = 0.2) + 
    stat_smooth(method = 'lm', se = F) +
    theme_bw() +
    labs(x = prot, y = ph, 
         title = paste0(ph, " - ", prot))  +
    scale_color_manual(values = my_colors)
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
