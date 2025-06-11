
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
my_colors <- c("#F65C00", "#eddb6d", "#006DB3", "#00A072")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/")

library(ggplot2)
library(dplyr)
library(patchwork)
library(lubridate)
library(purrr)

out_basedir <- "results/pheno_batch2_prot_rm_outliers_4sd/"

d_wide <- read.delim("data/olink_clean_CVD+INF_rm_outliers_4sd.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))
d_wide$TP <- as.numeric(d_wide$TP)

################################################################################
# PCA on proteins with missing data (nipals)
################################################################################

num_pcs_80 <- c()
pca_per_tp <- data.frame()
for (tp in 1:4){
  tmp_wide <- as.data.frame(subset(d_wide[d_wide$TP == tp,], select = -c(TP, SampleID, ID)))
  row.names(tmp_wide) <- d_wide[d_wide$TP == tp, "ID"]
  #tmp_wide <- na.omit(tmp_wide)
  pca <- pcaMethods::pca(tmp_wide, method = 'nipals', nPcs = 50, center = T, scale = 'uv')
  cumulative_variance <- cumsum(pca@R2)
  num_pcs_80 <- c(num_pcs_80, which(cumulative_variance >= 0.80)[1])

  pca10 <- as.data.frame(pca@scores)[,1:10] %>%
  rownames_to_column(var = 'ID') 
  pca_per_tp <- rbind(pca_per_tp, data.frame(TP = tp, pca10))
}
max(as.numeric(num_pcs_80))

write.table(pca_per_tp, file = paste0(out_basedir, "olink_clean_CVD+INF_rm_outliers_4sd.PCA.txt"), quote = F, sep = "\t", row.names = FALSE)
################################################################################
# Technical covariates
################################################################################

#
# plates, well position. 
#
plate_pos <- read.delim("data/plate_well_position.txt", as.is = T, check.names = F, sep = "\t")
plate_pos$PlateID <- as.factor(plate_pos$PlateID)
plate_pos$Well_pos1 <- as.factor(plate_pos$Well_pos1)
plate_pos$Well_pos2 <- as.factor(plate_pos$Well_pos2)
plate_pos$WellID <- NULL
plate_pos <- plate_pos %>%
  separate_wider_delim(cols = 'SampleID', delim = "_", names = c("ID", "TP"))

tech_correl <-run_kruskal_test_each_TP(pca_per_tp, plate_pos)

bonf_cutoff <- 0.05/nrow(tech_correl)
# [1] 0.0004166667
tech_correl$Bonferroni_sign <- ifelse(tech_correl$KW_test_pval < bonf_cutoff, T, F)

plate_pos$TP <- as.numeric(plate_pos$TP)
tech_merged <- full_join(plate_pos, pca_per_tp, by = c("ID", "TP"))

pdf(paste0(out_basedir, "plots/correlations_with_covariates/well_pos_vs_PCs.pdf"))
ggplot(tech_merged, aes(x = Well_pos1, y = PC1, group = Well_pos1)) + geom_boxplot() + theme_bw() + facet_wrap(~TP)
dev.off()
write.table(tech_correl, file = paste0(out_basedir, "correlations_with_covariates/plate_vs_PCs_KW.txt"), quote = F, sep = "\t", row.names = FALSE)

#
# Season and storage time.
#


collect_date <- read.delim("data/blood_collection_date.txt", as.is = T, check.names = F, sep = "\t", colClasses = c("character", "character"))
collect_date$datavisitaodierna <- as.Date(collect_date$datavisitaodierna, format = "%d/%m/%Y")

shipment_date <- as.Date("24/09/2024", format = "%d/%m/%Y")
collect_date$storage_months<- interval(collect_date$datavisitaodierna, shipment_date) %/% months(1)
collect_date$storage_quarters <- round(collect_date$storage_months / 4)

getSeason <- function(DATES) {
  WS <- as.Date("2012-12-15", format = "%Y-%m-%d") # Winter Solstice
  SE <- as.Date("2012-3-15",  format = "%Y-%m-%d") # Spring Equinox
  SS <- as.Date("2012-6-15",  format = "%Y-%m-%d") # Summer Solstice
  FE <- as.Date("2012-9-15",  format = "%Y-%m-%d") # Fall Equinox
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2012-%m-%d"))
  
  ifelse (d >= WS | d < SE, "Winter",
          ifelse (d >= SE & d < SS, "Spring",
                  ifelse (d >= SS & d < FE, "Summer", "Autumn")))
}
collect_date$season <- getSeason(collect_date$datavisitaodierna)

collect_date$datavisitaodierna <- NULL

colnames(collect_date)[1] <- "ID"
collect_date$season <- factor(collect_date$season, levels = c("Winter", "Spring", "Summer", "Autumn"))

season_kw <-  run_kruskal_test_each_TP(pca_per_tp, collect_date[,c("ID", "season", "storage_quarters")])
season_lm <-  run_lm_each_TP(pca_per_tp, collect_date[,c("ID", "storage_months")])

pls <- plot_PCA_boxplot(pca_per_tp, collect_date[,c("ID", "storage_months")], 'PC1', 'storage_months')
pdf(paste0(out_basedir, "correlations_with_covariates/plots/PC_vs_storage_months.pdf"), width = 10, height = 5)
pls[[1]] + pls[[2]]
dev.off()

pls <- plot_PCA_boxplot(pca_per_tp, collect_date[,c("ID", "season", "storage_quarters")], 'PC1', 'season')
pdf(paste0(out_basedir, "correlations_with_covariates/plots/PC_vs_season.pdf"), width = 10, height = 5)
pls[[1]] + pls[[2]]
dev.off()



################################################################################
# Main covariates
################################################################################


pheno_0 <- read.delim("../phenotypes/categorical_phenotypes_combined_visit0.txt", sep = "\t", check.names = F, as.is = T)
pheno_0$ID <- gsub("X","", pheno_0$Patient_id)
pheno_0$Patient_id <- NULL
pheno_0 <- pheno_0[pheno_0$ID %in% d_wide$ID,]

factor_counts <- lapply(pheno_0, function(column) {
  if (length(unique(column)) < 5) {
    return(table(column, useNA = "ifany"))
  } else {
    return(NULL)
  }
})

for (col_name in names(factor_counts)) {
  cat("Factor levels for", col_name, ":\n")
  print(factor_counts[[col_name]])
  cat("\n")
}

# Remove factor phenotypes with < 10 observations per level
pheno_0$COVID_vaccine <- NULL
pheno_0$T1D_first_second <- NULL
pheno_0$Antibiotic_use_30_days <- NULL
pheno_0$caffe_frequenza <- NULL
pheno_0$lactose_intolerance <- NULL

# Kruskal Wallis test on categorical covariates
covar_kw_res <- run_kruskal_test_each_TP(pca_per_tp, pheno_0)
colnames(covar_kw_res)[4] <- "pval"

# Linear regression on Age and BMI
pheno_cont_0 <- read.delim("../phenotypes/continuous_phenotypes_combined_visit0.txt", sep = "\t", check.names = F, as.is = T)
pheno_cont_0$ID <- gsub("X","", pheno_cont_0$ID)
pheno_cont_0 <- pheno_cont_0[pheno_cont_0$ID %in% d_wide$ID,]
covar_lm_res <- run_lm_each_TP(pca_per_tp, pheno_cont_0)

combined_covar_res <- rbind(covar_kw_res, subset(covar_lm_res, select = -c(beta)))

combined_covar_res_5pcs <- combined_covar_res[combined_covar_res$PC %in% c("PC1", "PC2", "PC3", "PC4", "PC5"),]
combined_covar_res_5pcs <- combined_covar_res_5pcs[! combined_covar_res_5pcs$covariate %in% c("Age_category", "BMI_ranges", "ID"),]
num_tests_per_tp <- nrow(combined_covar_res_5pcs[combined_covar_res_5pcs$TP == '1',])
combined_covar_res_5pcs$bonf_sign <- ifelse(combined_covar_res_5pcs$pval < 0.05/num_tests_per_tp, T, F)

combined_covar_res_5pcs[combined_covar_res_5pcs$bonf_sign == T,]

combined_covar_res_5pcs$logp <- -log10(combined_covar_res_5pcs$pval)
combined_covar_res_5pcs$pval <- as.numeric(combined_covar_res_5pcs$pval)

write.table(combined_covar_res_5pcs, file = paste0(out_basedir, 'correlations_with_covariates/PC_vs_covariates.txt'), quote = F, sep = "\t", row.names = FALSE)

pdf(paste0(out_basedir, 'correlations_with_covariates/plots/corrplots_PC_vs_covariates.pdf'), width = 10, height  = 12)
color_limits <- range(combined_covar_res_5pcs$logp)
par(mfrow=c(2,2))
for (tp in 1:4){
  
  logp_mat <- my_pivot_wider(combined_covar_res_5pcs[combined_covar_res_5pcs$TP == tp,c("PC", "covariate", "logp")], row_names = 'PC', names_from = 'covariate', values_from = 'logp')
  pval_mat <- my_pivot_wider(combined_covar_res_5pcs[combined_covar_res_5pcs$TP == tp,c("PC", "covariate", "pval")], row_names = 'PC', names_from = 'covariate', values_from = 'pval')
  corrplot(as.matrix(logp_mat), col.lim = color_limits, is.corr = F, sig.level = c(0.05/num_tests_per_tp, 0.05), 
           p.mat = as.matrix(pval_mat), insig = "label_sig", tl.col = 'black', mar = c(1, 1, 1, 1), tl.cex = 2)
}
dev.off()


pheno_0 <- subset(pheno_0, select = -c(Age_category, BMI_ranges))
pheno_comb <- full_join(pheno_cont_0, pheno_0, by = c("ID"))



################################################################################
# Functions
################################################################################


run_kruskal_test_each_TP <- function(pca_per_tp, covariates, num_pcs = 10){
  kw_res <- data.frame(matrix(ncol = 4))
  colnames(kw_res) <- c("TP", "PC", "covariate", "KW_test_pval")
  cnt <- 1
  for (tp in 1:4){
    
    pca10 <- pca_per_tp[pca_per_tp$TP == tp,]
    pca10$TP <- NULL
    
    if('TP' %in% colnames(covariates)){
      tmp_covariates <- covariates[covariates$TP == tp,]
      tmp_covariates$TP = NULL
    } else {
      tmp_covariates <- covariates
    }
    
    merged <- inner_join(tmp_covariates, pca10, by = c("ID"))
    
    for (pc in colnames(pca10)[2:ncol(pca10)]){
      for (cov in colnames(tmp_covariates)[2:ncol(tmp_covariates)]){
        subs <- merged[,c(pc,cov) ]
        colnames(subs) <- c("PC", "cov")
        kw <- kruskal.test(PC ~ cov, data=subs)$p.value
        
        kw_res[cnt,] <- c(tp, pc, cov, kw)
        cnt <- cnt + 1
      }
    }
  }
  
  kw_res <- na.omit(kw_res) %>%
    mutate(across(-c(PC, covariate), as.numeric))   
  
  kw_res
}

run_lm_each_TP <- function(pca_per_tp, covariates, num_pcs = 10){
  lm_res <- data.frame(matrix(ncol = 5))
  colnames(lm_res) <- c("TP", "PC", "covariate", "beta", "pval")
  cnt <- 1
  for (tp in 1:4){
    pca10 <- pca_per_tp[pca_per_tp$TP == tp,]
    pca10$TP <- NULL
    
    if('TP' %in% colnames(covariates)){
      tmp_covariates <- covariates[covariates$TP == tp,]
      tmp_covariates$TP = NULL
    } else {
      tmp_covariates <- covariates
    }
    
    merged <- inner_join(tmp_covariates, pca10, by = c("ID"))
    
    for (pc in colnames(pca10)[2:ncol(pca10)]){
      for (cov in colnames(tmp_covariates)[2:ncol(tmp_covariates)]){
        #cat(tp, pc, cov, "\n")
        subs <- merged[,c(pc,cov) ]
        colnames(subs) <- c("PC", "cov")
        lm_fit <- lm(PC ~ cov, data = subs)
        lm_res[cnt,] <- c(tp, pc, cov, summary(lm_fit)$coefficients['cov', 1], summary(lm_fit)$coefficients['cov', 4])
        cnt <- cnt + 1
      }
    }
  }
  
  lm_res <- na.omit(lm_res) %>%
    mutate(across(-c(PC, covariate), as.numeric))   
  
  lm_res
}

plot_PCA_boxplot <- function(pca_per_tp, covariates, pc, cov){
  
  pca_all_tps <- data.frame()
  for (tp in 1:4){
    pca10 <- pca_per_tp[pca_per_tp$TP == tp,]
    pca10$TP <- NULL
    
    if('TP' %in% colnames(covariates)){
      tmp_covariates <- covariates[covariates$TP == tp,]
      tmp_covariates$TP = NULL
    } else {
      tmp_covariates <- covariates
    }
    
    merged <- inner_join(tmp_covariates, pca10, by = c("ID"))
    pca_all_tps <- rbind(pca_all_tps, cbind(tp, merged))
  }
  

  subs <- pca_all_tps[c("PC1", "PC2", pc, cov, "tp")]
  colnames(subs) <- c("PC1", "PC2", "PC", "cov", "TP")
  if (length(unique(subs$cov)) < 5){
    p1 <- ggplot(subs, aes(x = cov, y = PC, group = cov)) + geom_boxplot() + theme_bw() + xlab(cov) + ylab(pc) + facet_wrap(~TP)
    p2 <- ggplot(pca_all_tps, aes(x = PC1, y = PC2, color = season)) + geom_point() + scale_color_manual(values = my_colors) + theme_bw()
  } else {
    p1 <- ggplot(subs, aes(x = cov, y = PC)) + geom_point() + geom_smooth(method = 'lm')+ theme_bw() + xlab(cov) + ylab(pc) + facet_wrap(~TP)
    p2 <- ggplot(subs, aes(x = PC1, y = PC2, color = cov)) + geom_point() + theme_bw() + facet_wrap(~TP)
  }
  return(list(p1,p2))
}

regression_olink_pheno <- function(joined_data, prot, ph, scale = F, kw = F){
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
  b <- NA
  pval <- NA
  if (! kw || ! is_factor){
    lm_fit <- lm(prot ~ pheno, data = d)
    pval <- summary(lm_fit)$coefficients[2,4]
    b <- summary(lm_fit)$coefficients[2,1]
    
  } else if (is_factor) {
    pval <- kruskal.test(prot ~ pheno, data=d)$p.value
    b <- NA
  }
  return(list("pval" = pval, "est" = b, "n" = nrow(d)))
}



