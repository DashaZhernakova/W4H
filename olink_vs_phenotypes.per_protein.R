my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/data/")

library(ggplot2)
library(dplyr)
library(lme4)
library(grid)
library(gridExtra)
library(tibble)
library(corrplot)

d_wide <- read.delim("olink_clean_CVD+INF.txt", as.is = T, check.names = F, sep = "\t")

d_wide <- cbind( gsub("_.*", "", d_wide$SampleID), gsub(".*_", "", d_wide$SampleID), d_wide)
colnames(d_wide)[c(1,2)] <- c("ID", "TP")
d_wide$SampleID <- NULL

d_long <- read.delim("olink_clean_CVD+INF_long.txt", as.is = T, check.names = F, sep = "\t")

pheno_covar <- read.delim("../data/covariate_selection.txt", sep = "\t", check.names = F, as.is = T, colClasses = c(ID = "character"))

covar <- read.delim("../data/covariates_age_bmi_storage_preg.txt", sep = "\t", check.names = F, as.is = T, colClasses = c(ID = "character"))

################################################################################
# run a PCA, check how many PCs explain 80% of variance
################################################################################

tmp_wide <- as.data.frame(subset(d_wide[d_wide$TP == 1,], select = -c(TP,  ID)))
row.names(tmp_wide) <- d_wide[d_wide$TP == 1, "ID"]
tmp_wide <- na.omit(tmp_wide)

pca <- prcomp(tmp_wide, scale = T)
num_pcs_0.8 <- as.numeric(which(summary(pca)$importance["Cumulative Proportion",] > 0.8)[1])
num_pcs_0.8
#[1] 38

################################################################################
# Associations between proteins and covariates
################################################################################

covariate_names <- c("Age","BMI","storage_months")
joined_data <- full_join(full_join(covar, pheno_covar, by = c("ID")), d_wide, by = c("ID"), relationship = 'one-to-many')
joined_covar <- full_join(covar, pheno_covar, by = c("ID"))

# Make a dataframe with proteins adjusted for all covariates per visit
tmp1 <- cbind(1, regress_covariates(d_wide[d_wide$TP == 1,], covar))
tmp2 <- cbind(2, regress_covariates(d_wide[d_wide$TP == 2,], covar))
tmp3 <- cbind(3, regress_covariates(d_wide[d_wide$TP == 3,], covar))
tmp4 <- cbind(4, regress_covariates(d_wide[d_wide$TP == 4,], covar))
colnames(tmp1)[1] = colnames(tmp2)[1]= colnames(tmp3)[1] = colnames(tmp4)[1] = 'TP'
d_wide_adj_covar <- rbind(tmp1, tmp2, tmp3, tmp4)
joined_data_adj_covar <- full_join(subset(joined_covar, select = -c(storage_months, Age, BMI, Pregnancy_category)), d_wide_adj_covar, by = c("ID"), relationship = 'one-to-many')

# do not regress out pregnancy
tmp1 <- cbind(1, regress_covariates(d_wide[d_wide$TP == 1,], covar[,-5]))
tmp2 <- cbind(2, regress_covariates(d_wide[d_wide$TP == 2,], covar[,-5]))
tmp3 <- cbind(3, regress_covariates(d_wide[d_wide$TP == 3,], covar[,-5]))
tmp4 <- cbind(4, regress_covariates(d_wide[d_wide$TP == 4,], covar[,-5]))
colnames(tmp1)[1] = colnames(tmp2)[1]= colnames(tmp3)[1] = colnames(tmp4)[1] = 'TP'
d_wide_adj_covar_no_preg <- rbind(tmp1, tmp2, tmp3, tmp4)
joined_data_adj_covar_no_preg <- full_join(subset(joined_covar, select = -c(storage_months, Age, BMI)), d_wide_adj_covar_no_preg, by = c("ID"), relationship = 'one-to-many')


############## Linear mixed models ###############

res_prot_cov_lmm <- data.frame(matrix(nrow = (ncol(pheno_covar) - 1) * (ncol(d_wide) - 2), ncol = 5))
colnames(res_prot_cov_lmm) <- c("covariate", "prot", "pval", "est", "n")
cnt <- 1
for (ph in c(colnames(pheno_covar)[-1], "Pregnancy_category")){
  for (prot in colnames(d_wide)[-c(1,2)]){
    res <- lmm_olink_pheno_adj_covar(joined_data, prot, ph, covariate_names)
    res_prot_cov_lmm[cnt,] <- c(ph, prot, unlist(res))
    cnt <- cnt + 1
  }
}
res_prot_cov_lmm <- na.omit(res_prot_cov_lmm) %>%
  mutate(across(-c(covariate, prot, ), as.numeric)) 

num_tests <- num_pcs_0.8 * length(unique(res_prot_cov_lmm$covariate))
res_prot_cov_lmm$bonf_sign <- ifelse(res_prot_cov_lmm$pval < 0.05/num_tests, T, F)
res_prot_cov_lmm$BH_qval <- p.adjust(res_prot_cov_lmm$pval)

write.table(res_prot_cov_lmm, file = "../results/prot_vs_covar/proteins_vs_covariates_lmm.txt",  quote = F, sep = "\t", row.names = FALSE)

############## Linear models ###############

bonf_threshold <- 0.05 / (num_pcs_0.8 * (ncol(pheno_covar) - 1))


res_prot_cov_lm <- data.frame(matrix(nrow = (ncol(d_wide) - 2), ncol = 17))
colnames(res_prot_cov_lm) <- c("TP","covariate", "prot", 
                               "pval_1", "est_1", "n_1", 
                               "pval_2", "est_2", "n_2",
                               "pval_3", "est_3", "n_3",
                               "pval_4", "est_4", "n_4",
                               "num_nominally_sign", "num_bonf_sign")
cnt <- 1
for (ph in c(colnames(pheno_covar)[-1], "Pregnancy_category")){
  for (prot in colnames(d_wide)[-c(1,2)]){
    tmp_all_tp <- c()
    for (tp in 1:4){
      res <- lm_olink_pheno_adj_covar(joined_data, tp, prot, ph, covariate_names)
      tmp_all_tp <- c(tmp_all_tp, unlist(res))
    }
    res_prot_cov_lm[cnt,] <- c(tp, ph, prot, tmp_all_tp, 
                               sum(tmp_all_tp[c(1,4,7,10)] < 0.05),
                               sum(tmp_all_tp[c(1,4,7,10)] < bonf_threshold))
    
    cnt <- cnt + 1
  }
}
res_prot_cov_lm <- na.omit(res_prot_cov_lm) %>%
  mutate(across(-c(covariate, prot, ), as.numeric)) 
res_prot_cov_lm$min_pval <- apply(res_prot_cov_lm[,paste("pval", seq(1:4), sep = "_")], 1, min)
res_prot_cov_lm <- res_prot_cov_lm[order(res_prot_cov_lm$pval),]

write.table(res_prot_cov_lm, file = "../results/prot_vs_covar/proteins_vs_covariates_lm.txt",  quote = F, sep = "\t", row.names = FALSE)


############## Specific phenotypes ############### 
ph = 'Pregnancy_category'
ph = "Pill_use"
subs <- res_prot_cov_lmm[res_prot_cov_lmm$covariate == ph,]
subs$BH_qval <- p.adjust(subs$pval, method = 'BH')
subs$bonf_sign <- ifelse(subs$pval < 0.05/num_pcs_0.8, T, F)

############## Plotting ############### 

### boxplots protein - covariates
tmp <- full_join(res_prot_cov_lm, res_prot_cov_lmm, by = c("covariate", "prot"), relationship = "many-to-one", suffix = c(".lm", ".lmm"))

lm_sign <- res_prot_cov_lm[res_prot_cov_lm$num_bonf_sign > 0 & res_prot_cov_lm$num_nominally_sign > 1,]
lmm_sign <- res_prot_cov_lmm[res_prot_cov_lmm$bonf_sign == T,]

plot_list <- list()
for (i in 1:nrow(lm_sign)){
  prot = lm_sign[i, 'prot']
  cov = lm_sign[i, 'covariate']
  plot_list[[i]] <- make_pheno_prot_plot(joined_data_adj_covar_no_preg, prot, cov, res_lm = res_prot_cov_lm, ylabel = paste0(prot, " adjusted"))
}

plot_list[[i + 1]] <- make_pheno_prot_plot(joined_data_adj_covar_no_preg, 'MZB1', 'Pregnancy_category', res_lm = res_prot_cov_lm, ylabel = paste0("MZB1", " adjusted"))

pdf('../plots/prot_vs_covar/prot_covar_boxplots.pdf', width = 15, height = 15)
grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)  
dev.off()


# corrplot for lmm
res_prot_cov_lmm$logp <- -log10(res_prot_cov_lmm$pval)
prots <- res_prot_cov_lmm[res_prot_cov_lmm$pval < 0.005,"prot"]
subs <- res_prot_cov_lmm[res_prot_cov_lmm$prot %in% prots,]
subs[subs$est > 1,"est"] <- 1
subs[subs$est < -1,"est"] <- -1
lmm_ests <- my_pivot_wider(subs, "prot", "covariate", "est")
lmm_logp <- my_pivot_wider(subs, "prot", "covariate", "logp")
lmm_pval <- my_pivot_wider(subs, "prot", "covariate", "pval")

p_th <- 0.05 / (num_pcs_0.8 * length(unique(res_prot_cov_lmm$covariate)))
color_limits <- range(subs$est)


pdf("../plots/prot_vs_covar/prot_vs_covariates_lmm_p0.005.pdf")
corrplot(as.matrix(lmm_ests), col.lim = color_limits, is.corr = F, sig.level = c(p_th, 0.01, 0.05), 
         p.mat = as.matrix(lmm_pval), insig = "label_sig", tl.col = 'black', pch.cex = 1)

dev.off()

################################################################################
#### Diet
################################################################################


covariate_names <- c("Age","BMI","storage_months", "Pregnancy_category")

diet <- read.delim("FFQ_visit0_clean.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))

joined_data <- full_join(full_join(covar, diet, by = c("ID")), d_wide, by = c("ID"), relationship = 'one-to-many')
joined_data_adj_covar <- full_join(diet, d_wide_adj_covar, by = c("ID"), relationship = 'one-to-many')

############## Linear mixed models ###############

res_prot_diet_lmm <- data.frame(matrix(nrow = (ncol(diet) - 1) * (ncol(d_wide) - 2), ncol = 5))
colnames(res_prot_diet_lmm) <- c("covariate", "prot", "pval", "est", "n")
cnt <- 1
for (ph in colnames(diet)[-1]){
  print(ph)
  for (prot in colnames(d_wide)[-c(1,2)]){
    res <- lmm_olink_pheno_adj_covar(joined_data, prot, ph, covariate_names)
    res_prot_diet_lmm[cnt,] <- c(ph, prot, unlist(res))
    cnt <- cnt + 1
  }
}
res_prot_diet_lmm <- na.omit(res_prot_diet_lmm) %>%
  mutate(across(-c(covariate, prot, ), as.numeric)) 

res_prot_diet_lmm <- res_prot_diet_lmm[res_prot_diet_lmm$n > 100,]
# number of dietary products: 20
#num_tests <- num_pcs_0.8 * length(unique(res_prot_diet_lm$covariate))
num_tests <- num_pcs_0.8 * 20

res_prot_diet_lmm$bonf_sign <- ifelse(res_prot_diet_lmm$pval < 0.05/num_tests, T, F)

write.table(res_prot_diet_lmm, file = "../results/prot_vs_covar/proteins_vs_diet_lmm.txt", quote = F, sep = "\t", row.names = FALSE)

############## Linear mixed models ###############

# number of dietary products: 20
bonf_threshold <- 0.05 / (num_pcs_0.8 * 20)


res_prot_diet_lm <- data.frame(matrix(nrow = (ncol(d_wide) - 2), ncol = 16))
colnames(res_prot_diet_lm) <- c("covariate", "prot", 
                               "pval_1", "est_1", "n_1", 
                               "pval_2", "est_2", "n_2",
                               "pval_3", "est_3", "n_3",
                               "pval_4", "est_4", "n_4",
                               "num_nominally_sign", "num_bonf_sign")
cnt <- 1
for (ph in colnames(diet)[-1]){
  for (prot in colnames(d_wide)[-c(1,2)]){
    tmp_all_tp <- c()
    for (tp in 1:4){
      res <- lm_olink_pheno_adj_covar(joined_data, tp, prot, ph, covariate_names)
      tmp_all_tp <- c(tmp_all_tp, unlist(res))
    }
    res_prot_diet_lm[cnt,] <- c(ph, prot, tmp_all_tp, 
                               sum(tmp_all_tp[c(1,4,7,10)] < 0.05),
                               sum(tmp_all_tp[c(1,4,7,10)] < bonf_threshold))
    
    cnt <- cnt + 1
  }
}
res_prot_diet_lm <- na.omit(res_prot_diet_lm) %>%
  mutate(across(-c(covariate, prot, ), as.numeric)) 

res_prot_diet_lm$min_pval <- apply(res_prot_diet_lm[,paste("pval", seq(1:4), sep = "_")], 1, min)

#View(res_prot_diet_lm[,c(1,2,ncol(res_prot_diet_lm) - 2, ncol(res_prot_diet_lm) - 1, ncol(res_prot_diet_lm) )])

write.table(res_prot_diet_lm, file = "../results/prot_vs_covar/proteins_vs_diet_lm.txt", quote = F, sep = "\t", row.names = FALSE)

############## Plotting ###############

lm_sign <- res_prot_diet_lm[res_prot_diet_lm$num_bonf_sign > 0 & res_prot_diet_lm$num_nominally_sign > 1,]
lmm_sign <- res_prot_diet_lmm[res_prot_diet_lmm$bonf_sign == T,]

plot_list <- list()

cnt <- 1
for(i in 1:nrow(lm_sign)){
  plot_list[[cnt]] <- make_pheno_prot_plot(joined_data_adj_covar, lm_sign[i, "prot"], lm_sign[i, "covariate"], res_lm = lm_sign, ylabel = paste0(prot, " adjusted"))
  cnt <- cnt + 1
}

for (i in 1:nrow(lmm_sign)){
  plot_list[[cnt]] <- make_pheno_prot_plot(joined_data_adj_covar, lmm_sign[i, "prot"], lmm_sign[i, "covariate"], res_lm = res_prot_diet_lm, ylabel = paste0(prot, " adjusted"))
  cnt <- cnt + 1
}

pdf('../plots/prot_vs_covar/prot_diet_boxplots.pdf', width = 15, height = 15)
grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)  
dev.off()

################################################################################
# Longitudinal covariates
################################################################################

pheno_tp <- read.delim("../../phenotypes/phenotypes_combined_per_tp.txt", sep = "\t", check.names = F, as.is = T)
pheno_tp$Patient_id <- gsub("X","", pheno_tp$Patient_id)
pheno_tp$Code <- gsub("X","", pheno_tp$Code)
pheno_tp <- pheno_tp[pheno_tp$Patient_id %in% d_wide$ID,]
pheno_tp$Visit_number <- as.numeric(pheno_tp$Visit_number)

d_wide$TP <- as.numeric(d_wide$TP)

pheno_tp$Medication_thyroid <- NULL
pheno_tp$Medication_NSAID <- NULL
pheno_tp$Medication_paracetamol <- NULL
pheno_tp <- pheno_tp[,!grepl("GR\\b", colnames(pheno_tp))]

colnames(pheno_tp) <- gsub("Visit_number", "TP", colnames(pheno_tp))
colnames(pheno_tp) <- gsub("Code", "SampleID", colnames(pheno_tp))
colnames(pheno_tp) <- gsub("Patient_id", "ID", colnames(pheno_tp))

# Set WHR of 079 to NA because it is a strong outlier
pheno_tp[pheno_tp$ID == '079', 'WHR'] <- NA

# log-transform dietary phenotypes
#pheno_tp <- pheno_tp %>%
#  mutate(across(matches("_3gg$"), ~ log(. + 1)))

joined_data <- full_join(full_join(covar, pheno_tp, by = c("ID")), d_wide, by = c("ID", "TP"))
joined_data_adj_covar <- full_join(full_join(covar, pheno_tp, by = c("ID")), d_wide_adj_covar, by = c("ID", "TP"))


covariate_names <- c("Age","BMI","storage_months", "Pregnancy_category")

correl_res_tp <- data.frame(matrix(nrow = (ncol(d_wide) -3) * (ncol(pheno_tp) -3), ncol = 5))
colnames(correl_res_tp) <- c("prot", "pheno",  "pval", "estimate", "N")
cnt <- 1
for (ph in colnames(pheno_tp)[4:ncol(pheno_tp)]) {
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res <- lmm_olink_pheno_adj_covar(joined_data, prot, ph, covariate_names)
    correl_res_tp[cnt,] <- c(prot, ph, unlist(res))
    cnt <- cnt + 1
  }
}

correl_res_tp <- na.omit(correl_res_tp) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 

num_tests <- num_pcs_0.8 * length(unique(correl_res_tp$pheno))
correl_res_tp$bonf_sign <- ifelse(correl_res_tp$pval < 0.05/num_tests, T, F)


write.table(correl_res_tp, file = "../results/prot_vs_covar/proteins_vs_long_covariates_lmm.txt", quote = F, sep = "\t", row.names = FALSE)

############## Plotting ###############

plot_pheno_prot_by_tp(d_wide_adj_covar, pheno_tp[pheno_tp$ColestGR_3gg > 0,], "CCL17", "ColestGR_3gg")

subs <- correl_res_tp[correl_res_tp$prot %in% unique(correl_res_tp[correl_res_tp$pval < 0.001, "prot"]), ]
subs[subs$estimate > 1,"estimate"] <- 1
subs[subs$estimate < -1,"estimate"] <- -1

lmm_pval <- my_pivot_wider(subs , row_names = "prot", names_from =  "pheno", values_from = 'pval')
lmm_est <- my_pivot_wider(subs , row_names = "prot", names_from =  "pheno", values_from = 'estimate')



p_th <- 0.05 / (num_pcs_0.8 * length(unique(correl_res_tp$pheno)))

color_lims_diet <- range(subs[subs$pheno != 'WHR',"estimate"])
diet_covs <- colnames(lmm_est)[grepl("_3gg", colnames(lmm_est))]

pdf("../plots/prot_vs_covar/prot_vs_longitudinal_covar_corrplot_diet.pdf")
corrplot(as.matrix(lmm_est[, colnames(lmm_est) %in% diet_covs]), is.corr = F, sig.level = c(p_th, 0.01, 0.05), 
         p.mat = as.matrix(lmm_pval[, colnames(lmm_pval) %in% diet_covs]), insig = "label_sig", tl.col = 'black', pch.cex = 1)
dev.off()

pdf("../plots/prot_vs_covar/prot_vs_longitudinal_covar_scatter_plots.pdf")
plot_pheno_prot_by_tp(d_wide_adj_covar, pheno_tp, 'THPO', 'FibreGR_3gg')

#
# Functions
#

lmm_olink_pheno_adj_covar <- function(joined_data, prot, ph, covariate_names, scale = F){
  d <- na.omit(joined_data[,c("ID", "TP", prot, ph, covariate_names)])
  colnames(d) <- c("ID", "TP", "prot", "pheno", covariate_names)
  is_factor <- length(unique(d$pheno)) < 3
  if ( length(unique(d$pheno)) < 2) return (list("pval" = NA, "est" = NA, "n" = NA))
  d <- na.omit(d)
  if (is_factor){
    d$pheno <- as.factor(d$pheno)
  } 
  
  if(scale){
    d$prot <- scale(d$prot)
    if (! is_factor) d$pheno <- scale(d$pheno)
  }
  fo_lmm <- as.formula(paste(" prot ~ pheno + TP +", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
  lmm_fit <- lmer(fo_lmm, data = d)
  est <-  summary(lmm_fit)$coefficients[2,1]
  fo_lmm_base <- as.formula(paste(" prot ~ TP + ", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
  lmm_fit_base <- lmer(fo_lmm_base, data = d)
  an <- suppressMessages(anova(lmm_fit, lmm_fit_base))
  pval <- an$`Pr(>Chisq)`[2]
  return(list("pval" = pval, "est" = est, "n" = nrow(d)))
}

lm_olink_pheno_adj_covar <- function(joined_data, tp, prot, ph, covariate_names, scale = F){
  d <- na.omit(joined_data[joined_data$TP == tp,c("ID", prot, ph, covariate_names)])
  colnames(d) <- c("ID", "prot", "pheno", covariate_names)
  is_factor <- length(unique(d$pheno)) < 3
  if ( length(unique(d$pheno)) < 2) return (list("pval" = NA, "est" = NA, "n" = NA))
  d <- na.omit(d)
  if (is_factor){
    d$pheno <- as.factor(d$pheno)
  } 
  
  if(scale){
    d$prot <- scale(d$prot)
    if (! is_factor) d$pheno <- scale(d$pheno)
  }
  fo <- as.formula(paste(" prot ~ pheno +", paste(covariate_names, collapse = "+")))
  lm_fit <- lm(fo, data = d)
  coef <- summary(lm_fit)$coefficients
  return(list("pval" = coef[2,4], "est" = coef[2,1], "n" = nrow(d)))
}




make_pheno_prot_plot <- function(joined_data, prot, ph, res_lm = NULL, ylabel = prot){
  if (! "TP" %in% colnames(joined_data)){
    d <- na.omit(joined_data[,c("ID", prot, ph)])
    colnames(d) <- c("SampleID", "prot", "pheno")
    is_factor <- length(unique(d$pheno)) < 3
    
    if(! is_factor){
      
      p <- ggplot(d, aes(x = pheno, y = prot)) + 
        geom_point() +
        geom_smooth(method = 'lm') +
        theme_bw() +
        xlab(ph) + ylab(ylabel)
    } else {
      
      d$pheno <- as.factor(d$pheno)
      
      p <- ggplot(d, aes(x = pheno, y = prot, group = pheno)) + 
        geom_boxplot() + 
        geom_jitter(width = 0.1) +
        theme_bw() +
        xlab(ph) + ylab(ylabel)
      
    }
  } else {
    d <- na.omit(joined_data[,c("TP","ID", prot, ph)])
    colnames(d) <- c("TP", "SampleID", "prot", "pheno")
    is_factor <- length(unique(d$pheno)) < 3
    
    if(! is_factor){
      
      p <- ggplot(d, aes(x = pheno, y = prot)) + 
        geom_point() +
        geom_smooth(method = 'lm') +
        theme_bw() +
        xlab(ph) + ylab(ylabel) +
        facet_wrap(~TP)
    } else {
      
      d$pheno <- as.factor(d$pheno)
      p <- ggplot(d, aes(x = pheno, y = prot, group = pheno)) + 
        geom_boxplot() + 
        geom_jitter(width = 0.1) +
        theme_bw() +
        xlab(ph) + ylab(ylabel) + 
        facet_wrap(~TP)
    }
      if (!is.null(res_lm)){
        
        if (ncol(res_lm) < 10){ # long format
          pvals <- res_lm[res_lm$prot == prot & res_lm$covariate == ph,]
          colnames(pvals) <- gsub("covariate","pheno",  colnames(pvals))
        } else { # wide format
          pvals <- res_lm[res_lm$prot == prot & res_lm$covariate == ph, c("covariate", "prot", paste("pval", seq(1:4), sep = "_"))] %>%
            pivot_longer(cols = paste("pval", seq(1:4), sep = "_"), names_to = "TP", values_to = "pval")
          colnames(pvals) <- gsub("covariate","pheno",  colnames(pvals))
          pvals$TP <- gsub("pval_", "", pvals$TP)
        }
        
        xpos = ifelse(is_factor, max(levels(d$pheno)), max(d$pheno))
        p <- p + geom_text(data = pvals, 
                    aes(x = xpos, y = max(d$prot) * 1.05, label = paste("lm p =", formatC(pval, digits = 3))),
                    color = "black",
                    size = 3,
                    hjust = 1)
      
    }
  }
  return(p)
}

plot_pheno_prot_by_tp <- function(d_wide, pheno, prot, ph){
  d_subs <- inner_join(d_wide[,c("ID", "TP", prot)], pheno[, c("ID", "TP", ph)], by = c("ID", "TP"))
  colnames(d_subs) <- c( "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.factor(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  p <- ggplot(d_subs, aes(x = pheno, y = prot, color = TP, group = TP)) + 
    geom_point() + 
    geom_smooth(method = 'lm') + 
    xlab(ph) + ylab(prot) +
    theme_bw() + scale_color_manual(values = my_colors)
  return(p)
}

regress_covariates <- function(data, covariates){
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
    covariates <- covariates[ids,]
    
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
