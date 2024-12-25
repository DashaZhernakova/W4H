library(dplyr)
library(lme4)

d <- read.delim("../../phenotypes/W4H_Trieste_Dati_volontarie_allsamples_noduplicates13122024.csv", 
                sep = ",", as.is = T, check.names = F, colClasses = c(sesso = "character"))


d <- d[,c("Record.ID", "Numero.Visita","Glucose_mg_dL", "TotChol_mg.dL", "HDLC_mg_dL", "LDLC_mg_dL", "Trigl_mg_dL", "AST_U_L", "ALT_U_L", "INS_mIU_L", "HOMA_IR", "HOMA_B", "PROG_ng_mL",  "LH_mlU_mL", "FSH_mlU_ml", "X17BES_pg_mL", "PRL_ng_mL", "Batch_hormones_measurement")]
colnames(d) <- gsub("_[mnpUu].*", "", colnames(d))
colnames(d) <- gsub("Numero.Visita", "TP", colnames(d))

d <- d[rowSums(!is.na(d[, 3:ncol(d)])) > 0, ]

table(d[,c("TP", "Batch_hormones")])


#
# Batch effects
#

batch_effects <- data.frame(matrix(nrow = ncol(d) - 3, ncol = 2))
colnames(batch_effects) <- c("pheno", "batch_effect_pval")
cnt <- 1
for (ph in colnames(d)[4: ncol(d)]){
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

d <- cbind(d[,c("Record.ID", "TP", "Batch_hormones")], d[,3: (ncol(d) - 1)])
#
# Distributions
#

plot_list = list()
for (ph in colnames(d)[4: ncol(d)]){
  plot_list[[ph]] <- ggplot(d, aes(x = !!ensym(ph))) + geom_density() + theme_bw()
  plot_list[[paste0(ph, "_log")]] <- ggplot(d, aes(x = log( !!ensym(ph) + 1) )) + geom_density() + theme_bw()
}
pdf("../plots/lipids_hormones_distributions_batch2.pdf", height = 20, width = 20)
grid.arrange(grobs = plot_list, ncol = 6, nrow = 5)  
dev.off()


#
# Log-transform all pheno
#

d_log <- d
d_log[, 4:ncol(d_log)] <- log(d[, 4:ncol(d)])



#
# Regress out batch effect from log-transformed data
#

d_adj <- d[,c("Record.ID", "TP")]

cnt <- 1
for (ph in colnames(d)[4: (ncol(d))]){
  print(ph)
  subs <- na.omit(d_log[,c("Record.ID", "TP", ph, "Batch_hormones")])
  colnames(subs)[3] <- 'pheno'
  lmm_fit <- lmer(pheno ~ Batch_hormones + (1|Record.ID), data = subs)
  subs[,ph] <- subs$pheno - predict(lmm_fit, re.form = NA)
  
  d_adj <- left_join(d_adj, subs[, c("Record.ID", "TP", ph)], by = c("Record.ID", "TP"))
}


#
# Write final tables 
#
write.table(d, file = "../../phenotypes/blood_pheno_13122024.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(d_log, file = "../../phenotypes/blood_pheno_13122024_log.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(d_adj, file = "../../phenotypes/blood_pheno_13122024_log_adj_batch.txt", quote = F, sep = "\t", row.names = FALSE)
