
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/data/")

library(ggplot2)
library(dplyr)
library(lme4)
library(grid)
library(gridExtra)

d_wide <- read.delim("olink_clean_CVD+INF.txt", as.is = T, check.names = F, sep = "\t")

d_wide <- cbind( gsub("_.*", "", d_wide$SampleID), gsub(".*_", "", d_wide$SampleID), d_wide)
colnames(d_wide)[c(1,2)] <- c("ID", "TP")

d_long <- read.delim("olink_clean_CVD+INF_long.txt", as.is = T, check.names = F, sep = "\t")


#
# Technical covariates: plates and well position
#

colnames(d_wide) <- gsub("Cardiometabolic","CVD", colnames(d_wide)) 

tmp_wide2 <- as.data.frame(d_wide[,-c(1,2,3)])
row.names(tmp_wide2) <- d_wide$SampleID
tmp_wide2 <- na.omit(tmp_wide2)

pca <- prcomp(tmp_wide2, scale = T)
pca10 <- as.data.frame(pca$x[,1:10]) %>%
  rownames_to_column(var = 'SampleID') %>%
  separate_wider_delim(cols = "SampleID", names = c("ID", "TP"), delim = "_")

plate_pos <- read.delim("plate_well_position.txt", as.is = T, check.names = F, sep = "\t")
plate_pos$PlateID <- as.factor(plate_pos$PlateID)
plate_pos$Well_pos1 <- as.factor(plate_pos$Well_pos1)
plate_pos$Well_pos2 <- as.factor(plate_pos$Well_pos2)

plate_pos <- plate_pos %>%
  separate_wider_delim(cols = 'SampleID', delim = "_", names = c("ID", "TP"))

tech_merged <- inner_join(plate_pos, pca10, by = c("ID", "TP"))

tech_correl <- data.frame(tp = character(), 
                          prot = character(),
                          ANOVA_well_pos1 = numeric(),
                          ANOVA_well_pos2 = numeric(),
                          ANOVA_plateID = numeric())
cnt <- 1
for (tp in 1:4){
  for (pc in colnames(pca10)[3:ncol(pca10)]){
    
    subs <- tech_merged[tech_merged$TP == tp, ]
    colnames(subs) <- gsub(pc, "PC", colnames(subs))
    aov_well1 <- summary(aov(PC ~ Well_pos1, data=subs))[[1]][["Pr(>F)"]][1]
    aov_well2 <- summary(aov(PC ~ Well_pos2, data=subs))[[1]][["Pr(>F)"]][1]
    aov_plate <- summary(aov(PC ~ PlateID, data=subs))[[1]][["Pr(>F)"]][1]
    tech_correl[cnt,] <- c(tp, pc, aov_well1, aov_well2, aov_plate)
    cnt <- cnt + 1
  }
}

tech_correl <- na.omit(tech_correl) %>%
  mutate(across(-c( prot), as.numeric))   


ggplot(tech_merged, aes(x = Well_pos1, y = PC1, group = Well_pos1)) + geom_boxplot() + theme_bw() + facet_wrap(~TP)


#
# Phenotypes
#

pheno_0 <- read.delim("../../phenotypes/phenotypes_combined_visit0.txt", sep = "\t", check.names = F, as.is = T)
pheno_tp <- read.delim("../../phenotypes/phenotypes_combined_per_tp.txt", sep = "\t", check.names = F, as.is = T)
pheno_0$Patient_id <- gsub("X","", pheno_0$Patient_id)
pheno_0 <- pheno_0[pheno_0$Patient_id %in% d_wide$ID,]



# First we focus on the most prevalent phenotypes and check if they affect PCA of proteins
pheno_0 <- pheno_0[,c("Patient_id", "Age", "BMI", "Smoker", "Pregnancy_category", "Pill_use", "caffe", "alcol")]

d_wide_1 <- d_wide[d_wide$TP == 1,]
d_wide_1$TP <- NULL
d_wide_1$SampleID <- NULL

tmp_d_wide_1 <- d_wide_1
row.names(tmp_d_wide_1) <- d_wide_1$ID
tmp_d_wide_1$ID <-NULL

pca <-  prcomp(na.omit(tmp_d_wide_1), scale = T)
pca10 <- as.data.frame(pca$x[,1:10]) %>%
  rownames_to_column(var = 'ID')

joined_data_pca <- left_join(pca10, pheno_0, by = c("ID" = "Patient_id"))

correl_pc_res <- data.frame(matrix(nrow = (ncol(pca10) - 1) * (ncol(pheno_0) -1), ncol = 5))
colnames(correl_pc_res) <- c("PC", "pheno", "pval", "estimate", "N")
cnt <- 1
for (ph in colnames(pheno_0)[2:ncol(pheno_0)]) {
  for (pc in colnames(pca10)[2:ncol(pca10)]){
    res <- regression_olink_pheno(joined_data_pca, pc, ph, scale = T)
    correl_pc_res[cnt,] <- c(pc, ph, unlist(res))
    cnt <- cnt + 1
  }
}

correl_pc_res <- na.omit(correl_pc_res) %>%
  mutate(across(-c(pheno, PC), as.numeric)) 


correl_pc_res_subs <- run_permutation_adjustment(pheno_0, pca10, correl_pc_res, p_cutoff = 0.01)
correl_pc_res_subs$qval <- qvalue(p = run_permutation_adjustment$padj_perm)$qvalues
write.table(correl_pc_res_subs, file = '../results/PC_vs_pheno_assoc.txt',quote = F, sep = "\t", row.names = FALSE)


#
# Plotting
#
pca_to_plot <- pca10

ph='caffe'
pca_to_plot$pheno <- tmp_pheno[row.names(pca_to_plot),ph]
pca_to_plot$pheno <- as.factor(pca_to_plot$pheno)

ggplot(pca_to_plot, aes(x = PC1, y = PC2, col = pheno)) + 
  geom_point() +
  theme_bw() +
  scale_color_discrete("Do you drink coffee?", labels=c("Yes","No"))

ggplot(pca_to_plot, aes(x = pheno, y = PC1, group = pheno)) + 
  geom_boxplot() +
  theme_bw() + xlab("Do you drink coffee?") + scale_x_discrete(labels=c("Yes", "No"))


ph=sym('Pregnancy_category')

ggplot(joined_data_pca, aes(x = !!ph, y = PC3, group = !!ph)) + 
  geom_boxplot() +
  theme_bw() 

ph=sym('BMI')

ggplot(joined_data_pca, aes(x = PC2, y = PC3, col = !!ph)) + 
  geom_point() +
  theme_bw() 





# Read the pheno data again and look for associations between phenotypes and proteins
# Count phenotype values remove some of the phenotypes with low counts
pheno_0 <- read.delim("../../phenotypes/phenotypes_combined_visit0.txt", sep = "\t", check.names = F, as.is = T)
pheno_tp <- read.delim("../../phenotypes/phenotypes_combined_per_tp.txt", sep = "\t", check.names = F, as.is = T)
pheno_0$Patient_id <- gsub("X","", pheno_0$Patient_id)
pheno_0 <- pheno_0[pheno_0$Patient_id %in% d_wide$ID,]

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
pheno_0$Ex_smoker <- NULL
pheno_0$COVID_vaccine <- NULL
pheno_0$T1D_first_second <- NULL
pheno_0$Antibiotic_use_30_days <- NULL
pheno_0$caffe_frequenza <- NULL
pheno_0$lactose_intolerance <- NULL
pheno_0$Passive_smoking <- NULL

pheno_0$Weight <- NULL # remove weight, because BMI is already there
pheno_0$Height <- NULL

pheno_0$CHOL_first_second <- NULL
#pheno_0$COVID_diagnosis <- NULL
#pheno_0$alcol_frequenza <- NULL

d_wide_1 <- d_wide[d_wide$TP == 1,]
d_wide_1$TP <- NULL
d_wide_1$SampleID <- NULL

joined_data <- left_join(d_wide_1, pheno_0, by = c("ID" = "Patient_id"))

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

correl_res <- correl_res[order(correl_res$pval),]
correl_res_subs <- run_permutation_adjustment(pheno_0, d_wide_1, correl_res, p_cutoff = 0.01)
correl_res_subs$qval <- qvalue(p = correl_res_subs$padj_perm)$qvalues
write.table(correl_res_subs, file = "../results/correl_prot_pheno_baseline_p0.01.v2.txt", quote = F, sep = "\t", row.names = FALSE)

#
# plotting
#
prot = 'WASF1'
ph = 'caffe'
make_pheno_prot_plot(joined_data, prot, ph)


pdf("../plots/prot-pheno_corrplots_scaled.pdf")
par(mfrow = c(1,2))
prot_subs <- unique(correl_res[correl_res$pval < 1e-3, "prot"])
correl_res_wide_subs <- correl_res[correl_res$prot %in% prot_subs, c("prot", "pheno", "logp")] %>%
  pivot_wider(names_from = pheno, values_from = logp)
correl_res_wide_subs <- as.data.frame(correl_res_wide_subs)
row.names(correl_res_wide_subs) <- correl_res_wide_subs$prot
correl_res_wide_subs$prot <- NULL
corrplot::corrplot(as.matrix(correl_res_wide_subs), is.corr = FALSE, tl.col = 'black')

correl_res_wide_subs <- correl_res[correl_res$prot %in% prot_subs, c("prot", "pheno", "estimate")] %>%
  pivot_wider(names_from = pheno, values_from = estimate)
correl_res_wide_subs <- as.data.frame(correl_res_wide_subs)
row.names(correl_res_wide_subs) <- correl_res_wide_subs$prot
correl_res_wide_subs$prot <- NULL
corrplot::corrplot(as.matrix(correl_res_wide_subs), is.corr = FALSE, tl.col = 'black')

dev.off()


#
# Check the correlation of lipids and hormones with these phenotypes
#

lipids <- read.delim("../../phenotypes/blood_pheno_03102024_log.txt", as.is = T, check.names = F, sep = "\t")
lipids$Record.ID <- gsub("ID_", "", lipids$Record.ID)
colnames(lipids)[1] <- "ID"
lipids <- cbind(SampleID = paste0(lipids$ID, "_", lipids$TP), lipids)
lipids$Age <- NULL
lipids <- lipids[lipids$TP == 1,]
lipids$SampleID = NULL
lipids$TP = NULL

joined_data <- left_join(lipids, pheno_0, by = c("ID" = 'Patient_id'))

correl_res_lip <- data.frame(matrix(nrow = (ncol(lipids) -1) * (ncol(pheno_0) -1), ncol = 5))
colnames(correl_res_lip) <- c("prot", "pheno", "pval", "estimate", "N")
cnt <- 1
for (ph in colnames(pheno_0)[2:ncol(pheno_0)]) {
  for (prot in colnames(lipids)[2:ncol(lipids)]){
    res <- regression_olink_pheno(joined_data, prot, ph, scale = T)
    correl_res_lip[cnt,] <- c(prot, ph, unlist(res))
    cnt <- cnt + 1
  }
}

correl_res_lip <- na.omit(correl_res_lip) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 

correl_res_lip$BH_qval <- p.adjust(correl_res_lip$pval, method = 'BH')


# Association adjusted for age and BMI
joined_data <- left_join(d_wide_1, pheno_0, by = c("ID" = "Patient_id"))

correl_res <- data.frame(matrix(nrow = (ncol(d_wide_1) -1) * (ncol(pheno_0) -1), ncol = 5))
colnames(correl_res) <- c("prot", "pheno", "pval", "estimate", "N")
cnt <- 1
for (ph in colnames(pheno_0)[2:ncol(pheno_0)]) {
  if (! ph %in% c("Age", "BMI")) {
  for (prot in colnames(d_wide_1)[2:ncol(d_wide_1)]){
    res <- regression_olink_pheno_adj_age_bmi(joined_data, prot, ph, scale = T)
    correl_res[cnt,] <- c(prot, ph, unlist(res))
    cnt <- cnt + 1
  }
  }
}

correl_res <- na.omit(correl_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 

correl_res <- correl_res[order(correl_res$pval),]
correl_res$logp <- -1*log10(correl_res$pval)

correl_res_perm <- run_permutation_adjustment(pheno_0, d_wide_1, correl_res, p_cutoff = 1, num_permutations = 100)
correl_res_subs$qval <- qvalue(p = correl_res_subs$padj_perm)$qvalues
write.table(correl_res_subs, file = "../results/correl_prot_pheno_baseline_adj_age_bmi.v2.txt", quote = F, sep = "\t", row.names = FALSE)



#
# PERMANOVA
#
library(vegan)
tmp_prot <- joined_data[,colnames(d_wide_1)]
row.names(tmp_prot) <- tmp_prot$ID
tmp_prot$ID <- NULL

eucl_dist <- dist(tmp_prot, method = "euclidean", 
                      diag = TRUE, upper = TRUE)

tmp_pheno <- joined_data[,colnames(pheno_0)[2:ncol(pheno_0)]]
row.names(tmp_pheno) <- joined_data$ID

adonis_pvals <- c()
for (ph in colnames(tmp_pheno)){
  non_na <- which(!is.na(tmp_pheno[,ph]))
  eucl_dist <- dist(tmp_prot[non_na,], method = "euclidean", 
                    diag = TRUE, upper = TRUE)
  ad <- adonis2(eucl_dist ~ tmp_pheno[non_na,ph],permutations=10000)
  adonis_p <-ad$`Pr(>F)`[1]
  adonis_pvals <- c(adonis_pvals, adonis_p)
  if (adonis_p < 0.05) cat(ph, adonis_p, "\n")
}

# coffee is nominally significant



#
# Longitudinal
#

pheno_tp <- read.delim("../../phenotypes/phenotypes_combined_per_tp.txt", sep = "\t", check.names = F, as.is = T)
pheno_tp$Patient_id <- gsub("X","", pheno_tp$Patient_id)
pheno_tp$Code <- gsub("X","", pheno_tp$Code)
pheno_tp <- pheno_tp[pheno_tp$Patient_id %in% d_wide$ID,]

factor_counts <- lapply(pheno_tp, function(column) {
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

pheno_tp$Medication_thyroid <- NULL
pheno_tp$Medication_NSAID <- NULL
pheno_tp$Medication_paracetamol <- NULL
pheno_tp <- pheno_tp[,!grepl("_3gg", colnames(pheno_tp))]


colnames(pheno_tp) <- gsub("Visit_number", "TP", colnames(pheno_tp))
colnames(pheno_tp) <- gsub("Code", "SampleID", colnames(pheno_tp))


lmm_pheno_prot(d_wide, pheno_tp, "PROK1", "CarboidratiGR_3gg")

correl_res_tp <- data.frame(matrix(nrow = (ncol(d_wide) -3) * (ncol(pheno_tp) -3), ncol = 5))
colnames(correl_res_tp) <- c("prot", "pheno",  "estimate", "pval", "N")
cnt <- 1
for (ph in colnames(pheno_tp)[4:ncol(pheno_tp)]) {
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res <- lmm_pheno_prot(d_wide, pheno_tp, prot, ph)
    correl_res_tp[cnt,] <- c(prot, ph, unlist(res))
    cnt <- cnt + 1
  }
}

correl_res_tp <- na.omit(correl_res_tp) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 

correl_res_tp$BH_qval <- p.adjust(correl_res_tp$pval, method = "BH")
correl_res_tp$logp <- -1*log10(correl_res_tp$pval)

correl_res_tp <- correl_res_tp[!correl_res_tp$pheno %in% c("Medication_NSAID", "Medication_paracetamol"),]

pdf("../plots/prot-pheno_corrplots_longitudinal.pdf")
par(mfrow = c(1,2))
prot_subs <- unique(correl_res_tp[correl_res_tp$pval < 0.001, "prot"])
correl_res_wide_subs <- correl_res_tp[correl_res_tp$prot %in% prot_subs, c("prot", "pheno", "logp")] %>%
  pivot_wider(names_from = pheno, values_from = logp)
correl_res_wide_subs <- as.data.frame(correl_res_wide_subs)
row.names(correl_res_wide_subs) <- correl_res_wide_subs$prot
correl_res_wide_subs$prot <- NULL
corrplot::corrplot(as.matrix(correl_res_wide_subs), is.corr = FALSE, tl.col = 'black')

correl_res_wide_subs <- correl_res_tp[correl_res_tp$prot %in% prot_subs, c("prot", "pheno", "estimate")] %>%
  pivot_wider(names_from = pheno, values_from = estimate)
correl_res_wide_subs <- as.data.frame(correl_res_wide_subs)
row.names(correl_res_wide_subs) <- correl_res_wide_subs$prot
correl_res_wide_subs$prot <- NULL
corrplot::corrplot(as.matrix(correl_res_wide_subs), is.corr = FALSE, tl.col = 'black')

dev.off()

plot_pheno_prot_by_tp(d_wide, pheno_tp, 'HK2', 'CaloriesGR')


#
# Functions
#


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


regression_olink_pheno_adj_age_bmi <- function(joined_data, prot, ph, scale = F){
  d <- na.omit(joined_data[,c("ID", "Age", "BMI", prot, ph)])
  colnames(d) <- c("SampleID", "Age", "BMI", "prot", "pheno")
  is_factor <- length(unique(d$pheno)) < 3
  d <- na.omit(d)
  if (is_factor){
    d$pheno <- as.factor(d$pheno)
  } 
  
  if(scale){
    d$prot <- scale(d$prot)
    if (! is_factor) d$pheno <- scale(d$pheno)
  }
  
  lm_fit <- lm(prot ~ pheno + Age + BMI, data = d)
  coef <- summary(lm_fit)$coefficients
  return(list("pval" = coef[2,4], "est" = coef[2,1], "n" = nrow(d)))
}


make_pheno_prot_plot <- function(joined_data, prot, ph){
  d <- na.omit(joined_data[,c("ID", prot, ph)])
  colnames(d) <- c("SampleID", "prot", "pheno")
  is_factor <- length(unique(d$pheno)) < 3
  
  if(! is_factor){

    p <- ggplot(d, aes(x = pheno, y = prot)) + 
      geom_point() +
      geom_smooth(method = 'lm') +
      theme_bw() +
      xlab(ph) + ylab(prot)
  } else {
    
    d$pheno <- as.factor(d$pheno)

    p <- ggplot(d, aes(x = pheno, y = prot, group = pheno)) + 
      geom_boxplot() + 
      geom_jitter(width = 0.1) +
      theme_bw() +
      xlab(ph) + ylab(prot)
      
  }
  return(p)
}



lmm_pheno_prot <- function(d_wide, pheno, prot, ph, scaled = F){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  is_factor <- length(unique(d$pheno)) < 3
  
  if (scaled){
    d_subs$prot <- scale(d_subs$prot)
    if (!is_factor) d_subs$pheno <- scale(d_subs$pheno)
  }
  model <- lmer(prot ~  pheno + TP + (1|ID), data = d_subs)
  est <- summary(model)$coefficients[1, "Estimate"]
  model0 <- lmer(prot ~  TP + (1|ID), data = d_subs)
  an <- anova(model, model0)
  pval <- an$`Pr(>Chisq)`[2]
  return(list(estimate = est, pval = pval, n = nrow(d_subs)))
}


plot_pheno_prot_by_tp <- function(d_wide, pheno, prot, ph, log_transform = T){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.factor(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (log_transform){
    d_subs$pheno <- log(d_subs$pheno)
  }
  
  p <- ggplot(d_subs, aes(x = pheno, y = prot, color = TP, group = TP)) + 
    geom_point() + 
    geom_smooth(method = 'lm') + 
    xlab(ph) + ylab(prot) +
    theme_bw() + scale_color_manual(values = my_colors)
  return(p)
}


run_permutation_adjustment <- function(pheno_0, prot_df, correl_res, p_cutoff = 0.01,  num_permutations = 1000, adjust_age_sex = FALSE){
  
  #set.seed(123) 
  
  # if we calculate the adjusted P per protein - phenotype pair, we can run the adjustment on the subset of nominally significant associations:
  correl_res_subs <- correl_res[correl_res$pval < p_cutoff, ]
  
  pb = txtProgressBar(min = 0, max = nrow(correl_res_subs), initial = 0) 
  stepi = 0
    
  for (l in 1:nrow(correl_res_subs)) {rc
    ph = correl_res_subs[l,'pheno']
    prot = correl_res_subs[l,'prot']
    permuted_pvals <- c()
    
    setTxtProgressBar(pb,stepi)
    stepi <- stepi + 1
    for (i in 1:num_permutations) {
        # get the sample ids for which we have this protein and this phenotype measured:
        prot_subs <- na.omit(prot_df[,c("ID", prot)])
        pheno_subs <- na.omit(pheno_0[,c("Patient_id", ph)])
        shared_ids <- intersect(prot_subs$ID, pheno_subs$Patient_id)
        prot_subs <- prot_subs[prot_subs$ID %in% shared_ids,]
        permuted_pheno <- pheno_subs[pheno_subs$Patient_id %in% shared_ids,]
        # shuffle ids:
        permuted_pheno$Patient_id <- sample(permuted_pheno$Patient_id)
        # combine protein and permuted pheno data
        joined_perm <- left_join(prot_subs, permuted_pheno, by = c("ID" = "Patient_id"))
        
        # Re-run the linear models with permuted phenotypes
        if (adjust_age_sex) {
          res <- regression_olink_pheno_adj_age_bmi(joined_perm, prot, ph, scale = T)
        } else {
          res <- regression_olink_pheno(joined_perm, prot, ph, scale = T)
        }
       permuted_pvals <- c(permuted_pvals, res[['pval']])
      }
      # get the P from the real analysis
      obs_pval <- correl_res_subs[l, "pval"]
      # calculate the adjusted P:
      correl_res_subs[l, "padj_perm"] <- sum(permuted_pvals < obs_pval) / num_permutations
  }
  cat ("\n")
  close(pb)
  return(correl_res_subs)
}




