pheno_old <- read.delim("../../phenotypes/blood_pheno_03102024_log.txt", as.is = T, check.names = F, sep = "\t")

pheno <- read.delim("../../phenotypes/blood_pheno_13122024_log.txt", as.is = T, check.names = F, sep = "\t")

pheno$Record.ID <- gsub("ID_", "", pheno$Record.ID)
colnames(pheno) <- gsub("Record.ID","ID",colnames(pheno))
pheno <- cbind(SampleID = paste0(pheno$ID, "_", pheno$TP), pheno)
pheno$Age <- NULL

covariates2 <- left_join(pheno[,c("ID", "TP", "Batch_hormones")], covariates, by = "ID")
covariates2 <- cbind(paste0(covariates2$ID, "_",covariates2$TP), covariates2)
colnames(covariates2)[1] <- "SampleID"

pheno$Batch_hormones <- NULL

storage <- na.omit(unique(covariates2[,c("ID", "storage_months")]))
tmp1 <- cbind(1, regress_covariates(d_wide[d_wide$TP == 1,], storage))
tmp2 <- cbind(2, regress_covariates(d_wide[d_wide$TP == 2,], storage))
tmp3 <- cbind(3, regress_covariates(d_wide[d_wide$TP == 3,], storage))
tmp4 <- cbind(4, regress_covariates(d_wide[d_wide$TP == 4,], storage))
colnames(tmp1)[1] = colnames(tmp2)[1]= colnames(tmp3)[1] = colnames(tmp4)[1] = 'TP'
d_wide_adj_storage <- rbind(tmp1, tmp2, tmp3, tmp4)
d_wide_adj_storage$SampleID <- paste0(d_wide_adj_storage$ID, "_", d_wide_adj_storage$TP)

d_wide_adj_storage2 <- regress_covariates_lmm(d_wide, unique(covariates2[,c("SampleID", "storage_months")]))


covariate_data <- subset(covariates2, select = -c(storage_months, TP, ID))



pheno_old <- read.delim("../../phenotypes/blood_pheno_03102024_log.txt", as.is = T, check.names = F, sep = "\t")


pheno_old$Record.ID <- gsub("ID_", "", pheno_old$Record.ID)
colnames(pheno_old) <- gsub("Record.ID","ID",colnames(pheno_old))
pheno_old <- cbind(SampleID = paste0(pheno_old$ID, "_", pheno_old$TP), pheno_old)
pheno_old$Age <- NULL



####### Old phenotypes; original direction prot ~ pheno with protein precorrected for storage time

lmm_res_orig_direction <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol =8))
colnames(lmm_res_orig_direction) <- c("prot", "pheno", "estimate", "pval", "se", "tval", "N", "N_unique")

covariate_data_tmp <- subset(covariate_data, select = -c(Batch_hormones))
covariate_data_tmp$SampleID <- gsub("_.*","", covariate_data_tmp$SampleID)
colnames(covariate_data_tmp) <- gsub("SampleID", "ID", colnames(covariate_data_tmp))
covariate_data_tmp <- unique(covariate_data_tmp)

cnt <- 1
for (ph in colnames(pheno_old)[4:ncol(pheno_old)]) {
  cat(ph, "\n")
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res <- lmm_pheno_prot_adj_covar(d_wide_adj_storage2, pheno_old, prot, ph, covariate_data_tmp, scale = T, adjust_timepoint = 'cubic')
    lmm_res_orig_direction[cnt,] <- c(prot, ph, unlist(res))
    
    cnt <- cnt + 1
  }
}

lmm_res_orig_direction <- na.omit(lmm_res_orig_direction) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
lmm_res_orig_direction$BH_pval <- p.adjust(lmm_res_orig_direction$pval)
lmm_res_orig_direction <- lmm_res_orig_direction[order(lmm_res_orig_direction$pval),]
nrow(lmm_res_orig_direction[lmm_res_orig_direction$BH_pval < 0.05,])


tmp <- full_join(lmm_res, lmm_res_orig_direction, by = c("pheno", "prot"))
tmp$logp.x <- -log10(tmp$pval.x)
tmp$logp.y <- -log10(tmp$pval.y)

ggplot (tmp, aes (x = estimate.x, y = estimate.y)) + geom_point() + 
  geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("LMM estimate old approach (38 significant)") + 
  ylab("LMM estimate precorrected for storage (36 significant)")


ggplot (tmp, aes (x = logp.x, y = logp.y)) + geom_point() + 
  geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("LMM pval old approach (38 significant)") + 
  ylab("LMM pval precorrected for storage (36 significant)")



####### Old phenotypes; reversed direction pheno ~ prot with protein precorrected for storage time



lmm_res_reversed_old_pheno <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno_old) -4), ncol = 8))
colnames(lmm_res_reversed_old_pheno) <- c("prot", "pheno", "estimate", "pval", "se", "tval", "N", "N_unique")
cnt <- 1
for (ph in colnames(pheno_old)[4:ncol(pheno_old)]) {
  cat(ph, "\n")
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res <- run_lmm_pheno_prot_reverse(d_wide_adj_storage2, pheno_old, prot, ph, subset(covariate_data, select = -c(Batch_hormones)), scale = T)
    lmm_res_reversed_old_pheno[cnt,] <- c(prot, ph, unlist(res))
    
    cnt <- cnt + 1
  }
}

lmm_res_reversed_old_pheno <- na.omit(lmm_res_reversed_old_pheno) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
lmm_res_reversed_old_pheno$BH_pval <- p.adjust(lmm_res_reversed_old_pheno$pval)
lmm_res_reversed_old_pheno <- lmm_res_reversed_old_pheno[order(lmm_res_reversed_old_pheno$pval),]
nrow(lmm_res_reversed_old_pheno[lmm_res_reversed_old_pheno$BH_pval < 0.05,])

tmp <- full_join(lmm_res_orig_direction, lmm_res_reversed_old_pheno, by = c("pheno", "prot"))
tmp$logp.x <- -log10(tmp$pval.x)
tmp$logp.y <- -log10(tmp$pval.y)

ggplot (tmp, aes (x = estimate.x, y = estimate.y)) + geom_point() + 
  geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("LMM estimate prot_adj ~ pheno_old (34 significant)") + 
  ylab("LMM estimate pheno_old ~ prot_adj (25 significant)")


ggplot (tmp, aes (x = logp.x, y = logp.y)) + geom_point() + 
  geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("LMM pval prot_adj ~ pheno_old (34 significant)") + 
  ylab("LMM pval pheno_old ~ prot_adj (25 significant)")



###### New direction pheno ~ prot, new phenotype data, proteins precorrected for storage time

lmm_res_reversed <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol = 8))
colnames(lmm_res_reversed) <- c("prot", "pheno", "estimate", "pval", "se", "tval", "N", "N_unique")
cnt <- 1
for (ph in colnames(pheno)[4:ncol(pheno)]) {
  cat(ph, "\n")
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res <- run_lmm_pheno_prot_reverse(d_wide_adj_storage2, pheno, prot, ph, covariate_data, scale = T)
    lmm_res_reversed[cnt,] <- c(prot, ph, unlist(res))
    
    cnt <- cnt + 1
  }
}

lmm_res_reversed <- na.omit(lmm_res_reversed) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
lmm_res_reversed$BH_pval <- p.adjust(lmm_res_reversed$pval)
lmm_res_reversed <- lmm_res_reversed[order(lmm_res_reversed$pval),]
nrow(lmm_res_reversed[lmm_res_reversed$BH_pval < 0.05,])


tmp <- full_join(lmm_res_reversed_old_pheno, lmm_res_reversed, by = c("pheno", "prot"))
tmp$logp.x <- -log10(tmp$pval.x)
tmp$logp.y <- -log10(tmp$pval.y)

ggplot (tmp, aes (x = estimate.x, y = estimate.y)) + geom_point() + 
  geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("LMM estimate old pheno pheno ~ prot (25 significant)") + 
  ylab("LMM estimate new pheno pheno ~ prot (27 significant)")

ggplot (tmp, aes (x = logp.x, y = logp.y)) + geom_point() + 
  geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("LMM pval old pheno pheno ~ prot (25 significant)") + 
  ylab("LMM pval new pheno pheno ~ prot (27 significant)")


###### New phenotypes, old direction, preadjusted for batch by Fabio
pheno_adj_batch <- read.delim("../../phenotypes/residuals_hormones.csv", as.is = T, sep = ",", check.names = F)
ids <- gsub("X", "",pheno_adj_batch$Code)
pheno_adj_batch <- cbind(ids, gsub("_.*", "", ids), gsub(".*_", "", ids), pheno_adj_batch)
pheno_adj_batch$Code <- NULL
pheno_adj_batch$Bacth <- NULL
pheno_adj_batch$Var.4 <- NULL
colnames(pheno_adj_batch)[1:3] <- c("SampleID", "ID", "TP")

x1 <- pheno
x2 <- pheno_adj_batch

row.names(x1) <- pheno$SampleID
row.names(x2) <- pheno_adj_batch$SampleID

plot(x1[ids, "PROG"], x2[ids, "PROG"], pch = 16, col = x1$TP, xlab = "old data", ylab = "new residuals")



run_lmm_pheno_prot_reverse <- function(prot_data, pheno_data, prot, ph, covariate_data, scale = T){

  d_subs <- inner_join(inner_join(prot_data[,c("SampleID", "ID", "TP", prot)], pheno_data[,c("SampleID" ,ph)], by = c("SampleID")),
                       covariate_data, by = c("SampleID"))
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)


  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }

  fo_lmm <- as.formula(paste("pheno ~ poly(TP, 3) + prot +", paste(colnames(covariate_data)[-1], collapse = "+"), "+ (1|ID)"))
  fo_lmm_base <- as.formula(paste("pheno ~ poly(TP, 3) + ", paste(colnames(covariate_data)[-1], collapse = "+"), "+ (1|ID)"))

  model <- lmer(fo_lmm, data = d_subs)
  model0 <- lmer(fo_lmm_base, data = d_subs)

  est <- summary(model)$coefficients["prot", "Estimate"]
  se <- summary(model)$coefficients["prot",2]
  tval <- summary(model)$coefficients["prot",3]
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  
  return(list(estimate = est, pval = pval, se = se, tval = tval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}




#################################################################


run_both_dir <- function(prot_data, pheno_data, prot, ph){
d_subs <- inner_join(inner_join(prot_data[,c("SampleID", "ID", "TP", prot)], pheno_data[,c("SampleID" ,ph)], by = c("SampleID")),
                     covariate_data, by = c("SampleID"))
colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
d_subs$TP <- as.numeric(d_subs$TP)
d_subs <- na.omit(d_subs)

d_subs$prot <- scale(d_subs$prot)
d_subs$pheno <- scale(d_subs$pheno)

lmm_fit1 <- lmer(prot ~ pheno +(1|ID), data = d_subs)
lmm_fit10 <- lmer(prot ~ (1|ID), data = d_subs)
est1 <- summary(lmm_fit1)$coefficients["pheno", "Estimate"]
var1 <- unlist(summary(lmm_fit1)$varcor)

an1 <- suppressMessages(anova(lmm_fit1, lmm_fit10))
pval1 <- an1$`Pr(>Chisq)`[2]


lmm_fit2 <- lmer(pheno ~ prot + (1|ID), data = d_subs)
lmm_fit20 <- lmer(pheno ~  (1|ID), data = d_subs)
est2 <- summary(lmm_fit2)$coefficients["prot", "Estimate"]
var2 <- unlist(summary(lmm_fit2)$varcor)
an2 <- suppressMessages(anova(lmm_fit2, lmm_fit20))
pval2 <- an2$`Pr(>Chisq)`[2]

return(c(prot, ph, est1, var1, pval1, est2, var2, pval2))

}


cmp <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol = 8))
colnames(cmp) <- c("prot", "pheno", "beta1", "var1", "pval1", "beta2", "var2", "pval2")
cnt <- 1
for (ph in colnames(pheno)[4:ncol(pheno)]) {
  cat(ph, "\n")
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res <- run_both_dir(d_wide_adj_storage2, pheno, prot, ph)
    cmp[cnt,] <- res
    
    cnt <- cnt + 1
  }
}

cmp <- na.omit(cmp) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 


ggplot(cmp,aes(x = beta1, y = beta2, color = var1)) + geom_point(alpha = 0.5) + geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("prot ~ pheno") + ylab("pheno ~ prot") + 
  guides(color=guide_legend(title="Variance of PROT \nexplained by ID")) +  scale_color_viridis_c()
ggplot(cmp,aes(x = beta1, y = beta2, color = var2)) + geom_point(alpha = 0.5) + geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("prot ~ pheno") + ylab("pheno ~ prot") + 
  guides(color=guide_legend(title="Variance of PHENO \nexplained by ID")) +  scale_color_viridis_c()



ggplot(cmp,aes(x = -log10(pval1), y = -log10(pval2), color = var1)) + geom_point(alpha = 0.5) + geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("prot ~ pheno") + ylab("pheno ~ prot") + guides(color=guide_legend(title="Variance of PROT \nexplained by ID")) +
  xlim(0,10) + ylim(0,10)+  scale_color_viridis_c()

ggplot(cmp,aes(x = -log10(pval1), y = -log10(pval2), color = var2)) + geom_point(alpha = 0.5) + geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("prot ~ pheno") + ylab("pheno ~ prot") + guides(color=guide_legend(title="Variance of PHENO \nexplained by ID")) +
  xlim(0,10) + ylim(0,10) +  scale_color_viridis_c()


ggplot(cmp,aes(x = -log10(pval1), y = -log10(pval2), color = pheno)) + geom_point(alpha = 0.5) + geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("prot ~ pheno") + ylab("pheno ~ prot") + guides(color=guide_legend(title="Variance of PHENO \nexplained by ID")) +
  xlim(0,10) + ylim(0,10) +  scale_color_viridis_d()


ggplot(cmp, aes(x = pheno, y = var2))  + geom_point() +theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylab("how much variance is explained by ID\nin LMM pheno ~ prot + TP + (1|ID)")

ggplot(cmp[cmp$prot %in% sample(unique(cmp$prot), 10),], aes(x = prot, y = var2, color = pheno))  + geom_point() +theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylab("how much variance is explained by ID\nin LMM pheno ~ prot + TP + (1|ID)")


### Compare how much variation is explained by ID for hormones and for proteins
icc_all <- data.frame()
cnt <- 1
for (ph in colnames(pheno)[4:ncol(pheno)]) {
  res <- get_ICC(pheno, ph)
  icc_all <- rbind(icc_all, c(ph, res, "pheno"))
}
for (prot in colnames(d_wide)[4:ncol(d_wide)]){
  res <- get_ICC(d_wide, prot)
  icc_all <- rbind(icc_all, c(prot, res, "prot"))
}
colnames(icc_all) <- c("feature", "icc", "type")
icc_all <- na.omit(icc_all) %>%
  mutate(across(-c(feature, type), as.numeric)) 


ggplot(icc_all, aes(x = feature, y = icc, color = type, group = type)) + geom_boxplot()




######
pheno_to_adj <- c("FSH", "INS","HOMA_B", "HOMA_IR", "LH", "PRL")

pheno_with_batch <- full_join(covariates2[,c("SampleID", "Batch_hormones")], pheno, by = "SampleID")
pheno_matched_means <- pheno_with_batch[,! colnames(pheno_with_batch) %in% pheno_to_adj]

for (ph in pheno_to_adj){
  pheno_matched_means <- full_join(pheno_matched_means, make_equal_means(pheno_with_batch, ph), by = "SampleID")
  
}
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


lmm_res_means <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol =8))
colnames(lmm_res_means) <- c("prot", "pheno", "estimate", "pval", "se", "tval", "N", "N_unique")

cnt <- 1
for (ph in colnames(pheno_old)[4:ncol(pheno_old)]) {
  cat(ph, "\n")
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res <- lmm_pheno_prot_adj_covar(d_wide, pheno_matched_means, prot, ph, covariates, scale = T, adjust_timepoint = 'cubic')
    lmm_res_means[cnt,] <- c(prot, ph, unlist(res))
    
    cnt <- cnt + 1
  }
}

lmm_res_means <- na.omit(lmm_res_means) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
lmm_res_means$BH_pval <- p.adjust(lmm_res_means$pval)

nrow(lmm_res_means[lmm_res_means$BH_pval < 0.05,])

tmp <- full_join(lmm_res, lmm_res_means, by = c("pheno", "prot"))
tmp$logp.x <- -log10(tmp$pval.x)
tmp$logp.y <- -log10(tmp$pval.y)

ggplot (tmp, aes (x = estimate.x, y = estimate.y)) + geom_point() + 
  geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("LMM estimate old pheno") + 
  ylab("LMM estimate new pheno means matched")


ggplot (tmp, aes (x = logp.x, y = logp.y)) + geom_point() + 
  geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("LMM pval old pheno") + 
  ylab("LMM pval new pheno means matched")



###### Pheno adjusted by Fabio

pheno_adj_batch <- read.delim("../../phenotypes/residuals_hormones.csv", as.is = T, sep = ",", check.names = F)
ids <- gsub("X", "",pheno_adj_batch$Code)
pheno_adj_batch <- cbind(ids,  pheno_adj_batch)
pheno_adj_batch$Code <- NULL
pheno_adj_batch$Bacth <- NULL
pheno_adj_batch$Var.2 <- NULL
colnames(pheno_adj_batch)[1] <- "SampleID"

pheno_adj <- full_join(pheno[,! colnames(pheno) %in% pheno_to_adj], pheno_adj_batch[, colnames(pheno_adj_batch) %in% c("SampleID", pheno_to_adj)], by = "SampleID")

lmm_res_RF <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol =8))
colnames(lmm_res_RF) <- c("prot", "pheno", "estimate", "pval", "se", "tval", "N", "N_unique")

cnt <- 1
for (ph in colnames(pheno_adj)[4:ncol(pheno_adj)]) {
  cat(ph, "\n")
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res <- lmm_pheno_prot_adj_covar(d_wide, pheno_adj, prot, ph, covariates, scale = T, adjust_timepoint = 'cubic')
    lmm_res_RF[cnt,] <- c(prot, ph, unlist(res))
    
    cnt <- cnt + 1
  }
}

lmm_res_RF <- na.omit(lmm_res_RF) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
lmm_res_RF$BH_pval <- p.adjust(lmm_res_RF$pval)



nrow(lmm_res_RF[lmm_res_RF$BH_pval < 0.05,])

tmp <- full_join(lmm_res_means, lmm_res_RF, by = c("pheno", "prot"))
tmp$logp.x <- -log10(tmp$pval.x)
tmp$logp.y <- -log10(tmp$pval.y)

ggplot (tmp, aes (x = estimate.x, y = estimate.y)) + geom_point() + 
  geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("LMM estimate new pheno means matched") + 
  ylab("LMM estimate new pheno adjusted using RF")


ggplot (tmp, aes (x = logp.x, y = logp.y)) + geom_point() + 
  geom_abline(slope = 1, lty = 2) + 
  theme_minimal() + xlab("LMM pval new pheno means matched") + 
  ylab("LMM pval new pheno adjusted using RF")

