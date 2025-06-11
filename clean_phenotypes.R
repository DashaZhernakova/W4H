library(dplyr)
library(lme4)
library(lubridate)
library(lmerTest)
library(patchwork)

setwd("/Users/Dasha/work/Sardinia/W4H/olink/data/")

d <- read.delim("../../phenotypes/W4H_Trieste_Dati_volontarie_allsamples_noduplicates13122024.csv", 
                sep = ",", as.is = T, check.names = F, colClasses = c(sesso = "character"))

d$data.accettazione. <- gsub("2024", "24", d$data.accettazione.)
d$storage_months<- interval(as.Date(d$data.accettazione, format = "%d/%m/%Y"), as.Date(d$data.analisi.profilo.ormonale, format = "%d/%m/%Y")) %/% months(1)

d <- d[,c("Record.ID", "Numero.Visita","Glucose_mg_dL", "TotChol_mg.dL", "HDLC_mg_dL", "LDLC_mg_dL", "Trigl_mg_dL", "AST_U_L", "ALT_U_L", "INS_mIU_L", "HOMA_IR", "HOMA_B", "PROG_ng_mL",  "LH_mlU_mL", "FSH_mlU_ml", "X17BES_pg_mL", "PRL_ng_mL", "Batch_hormones_measurement", "storage_months")]
colnames(d) <- gsub("_[mnpUu].*", "", colnames(d))
colnames(d) <- gsub("Numero.Visita", "TP", colnames(d))

d <- d[rowSums(!is.na(d[, 3:ncol(d)])) > 0, ]

table(d[,c("TP", "Batch_hormones")])


#
# Batch effects
#
all_hormones <- c("FSH", "INS","HOMA_B", "HOMA_IR", "LH", "PRL", "PROG", "X17BES")

batch_effects <- data.frame(matrix(nrow = ncol(d) - 3, ncol = 2))
colnames(batch_effects) <- c("pheno", "batch_effect_pval")
cnt <- 1
for (ph in all_hormones){
  subs <- d[,c("Record.ID", "TP", ph, "Batch_hormones")]
  colnames(subs)[3] <- 'pheno'
  lmm_fit <- lmer(pheno ~ TP + Batch_hormones + (1|Record.ID), data = subs)
  lmm_fit_0 <- lmer(pheno ~ TP + (1|Record.ID), data = subs)
  an <- anova(lmm_fit, lmm_fit_0)  
  pval <- an$`Pr(>Chisq)`[2] 
  batch_effects[cnt,] <- c(ph, pval)
  cnt <- cnt + 1
}
batch_effects$batch_effect_pval <- as.numeric(batch_effects$batch_effect_pval)

d <- cbind(d[,c("Record.ID", "TP", "Batch_hormones", "storage")], d[,3: (ncol(d) - 2)])

plot_list <- list()
for (ph in all_hormones){
  ph_name = sym(ph)
  plot_list[[ph]] <- ggplot(d, aes(x = Batch_hormones, y = !!ph_name, group = Batch_hormones)) + facet_wrap(~TP) + geom_boxplot() + theme_bw()
}
pdf("../plots/hormones_batch_effect_boxplots.pdf", height = 20, width = 10)
grid.arrange(grobs = plot_list, ncol = 2, nrow = 4)  
dev.off()

plot_list <- list()
for (ph in all_hormones){
  ph_name = sym(ph)
  plot_list[[ph]] <- ggplot(d, aes( x = !!ph_name, group = Batch_hormones, fill = Batch_hormones)) +  geom_density(alpha = 0.5) + theme_bw() + facet_wrap(~TP)
}
pdf("../plots/hormones_batch_effect_distributions.pdf", height = 20, width = 15)
grid.arrange(grobs = plot_list, ncol = 2, nrow = 4)  
dev.off()

#
# Distributions
#

plot_list = list()
for (ph in colnames(d)[5: ncol(d)]){
  plot_list[[ph]] <- ggplot(d, aes(x = !!ensym(ph))) + geom_density() + theme_bw()
  plot_list[[paste0(ph, "_log")]] <- ggplot(d, aes(x = log( !!ensym(ph) + 1) )) + geom_density() + theme_bw()
}
pdf("../plots/lipids_hormones_distributions_batch2.pdf", height = 20, width = 20)
grid.arrange(grobs = plot_list, ncol = 6, nrow = 5)  
dev.off()


#
# Log-transform certain phenotypes
#
pheno_to_log <- c("Trigl", "AST", "ALT", "FSH", "INS","HOMA_B", "HOMA_IR", "LH", "PRL", "PROG", "X17BES" )
d_log <- d
d_log[, pheno_to_log] <- log(d[, pheno_to_log])


#
# Correct for batch effect by matching the means of the affected hormones and regress out storage time
#
all_hormones <- c("FSH", "INS","HOMA_B", "HOMA_IR", "LH", "PRL", "PROG", "X17BES")

d_log$Record.ID <- paste0(gsub("ID_", "", d_log$Record.ID), "_", d_log$TP)
d_log <- cbind(d_log$Record.ID, gsub("_.*","",d_log$Record.ID), d_log[,2:ncol(d_log)])
colnames(d_log)[1:2] <- c("SampleID", "ID") 

# get the effect of storage time

for (ph in all_hormones){
  fo <- as.formula(paste(ph, "~ storage + (1|ID)"))
  lmm_fit <- lmer(fo, data = d_log)
  pval <- summary(lmm_fit)$coefficients["storage", 5]
  
  cat(ph, pval, "\n")
}


# regress storage time
all_hormones <- c("FSH", "INS","HOMA_B", "HOMA_IR", "LH", "PRL", "PROG", "X17BES")
d_log_adj <- regress_covariates_lmm(data = d_log[, c("SampleID", "ID", "TP", all_hormones)], covar_data = d_log[,c("SampleID", "storage"),drop =F])
d_log_adj$ID <- NULL
d_log_adj$TP <- NULL
d_log_adj <- full_join(d_log[,! colnames(d_log) %in% all_hormones], d_log_adj, by = "SampleID")

# correct for batch effect
pheno_to_adj_batch <- c("FSH", "INS","HOMA_B", "HOMA_IR", "LH", "PRL")
d_matched_means <- d_log_adj[,! colnames(d_log_adj) %in% pheno_to_adj_batch]


for (ph in pheno_to_adj_batch){
  d_matched_means <- full_join(d_matched_means, make_equal_means(d_log_adj, ph), by = "SampleID")
}
d_matched_means$Batch_hormones <- NULL
d_matched_means$storage <- NULL


d_log$Batch_hormones <- NULL
d_log$storage <- NULL
d_log_rm_outliers4 <- remove_outliers_dataframe(d_log,  sd_cutoff = 4)$cleaned_data

d_matched_means4 <- remove_outliers_dataframe(d_matched_means,  sd_cutoff = 5)$cleaned_data


#
# Write final tables 
#
write.table(d, file = "../../phenotypes/blood_pheno_13122024.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(d_log, file = "../../phenotypes/blood_pheno_13122024_log.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(d_matched_means, file = "../../phenotypes/blood_pheno_13122024_log_adj_batch_storage.txt", quote = F, sep = "\t", row.names = FALSE)

write.table(d_log_rm_outliers4, file = "../../phenotypes/blood_pheno_13122024_log_rm_outiers_4sds.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(d_matched_means4, file = "../../phenotypes/blood_pheno_13122024_log_adj_batch_storage_rm_outiers_4sds.txt.txt", quote = F, sep = "\t", row.names = FALSE)


make_equal_means <- function(pheno_with_batch, ph){
  b1 <- pheno_with_batch[pheno_with_batch$Batch_hormones == 'first', c("SampleID",  "Batch_hormones", ph)]
  b2 <- pheno_with_batch[pheno_with_batch$Batch_hormones == 'second', c("SampleID",  "Batch_hormones", ph)]
  
  m1 <- mean(b1[,ph], na.rm = T)
  m2 <- mean(b2[,ph], na.rm = T)
  mean_dif <- m2 - m1
  
  if (m1 > m2) stop("Mean of batch 1 is larger than the mean of batch2")
  
  b2[,ph] <- b2[,ph] - mean_dif
  
  adjusted <- rbind(b1, b2)
  
  pheno_name = sym(ph)
  p1 <- ggplot(pheno_with_batch, aes (x = !!pheno_name, group = Batch_hormones, fill = Batch_hormones)) + 
    geom_density(alpha = 0.3)+ geom_vline(xintercept = m1) + geom_vline(xintercept = m2) + theme_minimal()
  
  p2 <- ggplot(adjusted, aes (x = !!pheno_name, group = Batch_hormones, fill = Batch_hormones)) + 
    geom_density(alpha = 0.3)+ geom_vline(xintercept = m1)  + theme_minimal() 
  
  print(p1 + p2)
  adjusted$Batch_hormones <- NULL
  
  return (adjusted)
}




regress_covariates_lmm <- function(data, covar_data, covars_longitudinal = T){
  
  if (!"SampleID" %in% colnames(covar_data) & covars_longitudinal) {
    covar_data <- cbind(paste0(covar_data$ID, "_",covar_data$TP), covar_data)
    colnames(covar_data)[1] <- "SampleID"
  }
  
  d_adj <- data[,c("SampleID", "ID", "TP")]
  
  data[,"TP"] <- NULL
  covar_data[,"TP"] <- NULL
  
  cnt <- 1
  for (ph in colnames(data)[3: (ncol(data))]){
    if (covars_longitudinal){
      covar_data$ID <- NULL
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "SampleID"))
    } else {
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "ID"))
    }
    colnames(subs)[3] <- 'pheno'
    
    fo_lmm <- as.formula(paste("pheno ~ ", paste(colnames(covar_data)[-1], collapse = "+"), "+ (1|ID)"))
    
    
    lmm_fit <- lmer(fo_lmm, data = subs)
    subs[,ph] <- subs$pheno - predict(lmm_fit, re.form = NA)
    
    d_adj <- left_join(d_adj, subs[, c("SampleID", ph)], by = "SampleID")
  }
  return(d_adj)
}


remove_outliers_per_feature <- function(d, sd_cutoff = 4) {
  zscore <- scale(d)  # Compute z-scores
  d[abs(zscore) > sd_cutoff] <- NA  # Replace outliers with NA
  return(d)
}


remove_outliers_dataframe <- function(df, sd_cutoff = 4) {
  # Create a copy of the data frame to store the outlier mask
  outlier_mask <- df %>%
    mutate(across(.cols = -c(SampleID, ID, TP), .fns = ~abs(scale(.)) > sd_cutoff))
  
  # Apply the outlier removal
  df_cleaned <- df %>%
    mutate(across(.cols = -c(SampleID, ID, TP), .fns = ~remove_outliers_per_feature(., sd_cutoff)))
  
  # Return the cleaned data and the outlier mask
  list(cleaned_data = df_cleaned, outlier_mask = outlier_mask)
}

