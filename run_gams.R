library(mgcv)

out_basedir <- "results/pheno_batch2_prot_rm_outliers_5sd/"
d_wide <- read.delim("data/olink_clean_CVD+INF_rm_outliers_5sd.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))
covariates <- read.delim("data/covariates_age_bmi_storage_preg.txt", sep = "\t", check.names = F, as.is = T, colClasses = c(ID = "character"))

d_wide$TP <- as.numeric(d_wide$TP)

covariates[] <- lapply(covariates, function(col) {
  if (length(unique(col)) < 3) {
    return(factor(col))
  } else {
    return(col)
  }
})



################################################################################
# Protein vs TP LMM
################################################################################
all_prots <- colnames(d_wide)[4:ncol(d_wide)]

gam_res_prot_tp <- data.frame(matrix(nrow = length(all_prots), ncol = 3))
colnames(gam_res_prot_tp) <- c("prot", "pval", "pval_linearity")
cnt <- 1

n_points = 100
prot_trajs <- data.frame(matrix(nrow = length(all_prots) , ncol = n_points))
row.names(prot_trajs) <- all_prots
colnames(prot_trajs) <- seq(1,4, length.out = n_points)


for (prot in all_prots){
  res <- gam_prot_tp_adj_covar(d_wide, prot, covariates)
  gam_res_prot_tp[cnt,] <- c(prot, res$pval, res$p_nonlinear)
  prot_trajs[prot,] <- res$predicted
  cnt <- cnt + 1
}

gam_res_prot_tp <- gam_res_prot_tp %>%
  mutate(across(-c( prot), as.numeric)) 

gam_res_prot_tp$BH_pval <- p.adjust(gam_res_prot_tp$pval, method = 'BH')
gam_res_prot_tp <- gam_res_prot_tp[order(gam_res_prot_tp$pval),]

write.table(gam_res_prot_tp, file = paste0(out_basedir, "prot_vs_tp_gam_adj_age_bmi_preg_storage.txt"), quote = F, sep = "\t", row.names = FALSE)
nrow(gam_res_prot_tp[gam_res_prot_tp$BH_pval < 0.05,])

write.table(prot_trajs, file = paste0(out_basedir, "trajectories/prot_vs_tp_gam_adj_age_bmi_preg_storage.txt"), quote = F, sep = "\t", row.names = FALSE)



lmm_res_prot_tp <- read.delim(paste0(out_basedir, "prot_vs_tp_poly3_lmm_adj_age_bmi_preg_storage.txt"), sep = "\t", as.is = T, check.names = F)

tmp <- full_join(lmm_res_prot_tp, gam_res_prot_tp, by = "prot")

tmp$logp.x <- -log10(tmp$pval.x)
tmp$logp.y <- -log10(tmp$pval.y)

bh.x <- min(tmp[tmp$BH_pval.x < 0.05, "logp.x"])
bh.y <- min(tmp[tmp$BH_pval.y < 0.05, "logp.y"], na.rm = T)

tmp[tmp$logp.y == 'Inf', "logp.y"] <- 10
ggplot(tmp, aes(logp.x, logp.y)) + geom_point() + theme_minimal() + 
  geom_hline(yintercept = bh.y) + geom_vline(xintercept = bh.x) +
  xlab("LMM") + ylab("GAM")


signif_prots <- lmm_res_prot_tp[lmm_res_prot_tp$BH_pval < 0.05, "prot"]
prot_trajs_lmm <- read.delim(paste0(out_basedir, "trajectories/protein_trajectories_109signif.txt"), as.is = T, sep = "\t", check.names = F, row.names =1)

pdf(paste0(out_basedir,"trajectories/gam_vs_lmm.pdf"), height = 25, width = 20)
plot_list <- list()
cnt <- 1
for (prot in signif_prots){
  subs <- as.data.frame(t(rbind(colnames(prot_trajs), prot_trajs_lmm[prot,], prot_trajs[prot,])))
  colnames(subs) <- c("TP", "prot1", "prot2")
  subs <- subs %>% mutate_all(as.numeric)
  plot_list[[prot]] <- ggplot(subs) + geom_line(aes(TP, prot1), color = 'blue') + geom_line(aes(TP, prot2), color = 'red') + theme_minimal() + ylab(prot) + ggtitle("GAM - red, LMM -blue")
  cnt <- cnt + 1
  if (cnt == 21){
    grid.arrange(grobs = plot_list, ncol = 4, nrow = 5)
    plot_list <- list()
    cnt <- 1
  }
}
dev.off()


prot_trajs_gam_all <- prot_trajs
signif_prots <- gam_res_prot_tp[gam_res_prot_tp$BH_pval < 0.05, "prot"]

prot_trajs_gam <- prot_trajs_gam_all[signif_prots,]
# GAM-based deltas
prot_deltas <- data.frame(matrix(nrow = nrow(prot_trajs_gam), ncol = ncol(prot_trajs_gam) - 1))

row.names(prot_deltas) <- row.names(prot_trajs_gam)

for (c in 1:(ncol(prot_trajs_gam) - 1)){
  prot_deltas[,c] <- prot_trajs_gam[,(c + 1)] - prot_trajs_gam[,c]
}
method = 'pam_delta_GAM'
num_k = 6
#num_k = 15
pam_res <- cluster::pam(prot_deltas, num_k, diss = F)
cl <- pam_res$clustering
plot_clusters(cl, method, num_k, out_path = paste0(out_basedir,"trajectories/", method, "_k", num_k, "_131signif.pdf"), prot_trajs = prot_trajs_gam)
write.table(as.data.frame(cl), file = paste0(out_basedir, "trajectories/", method, "_k", num_k, "_131signif.txt"), quote = F, sep = "\t")

method = 'pam_traj_GAM'
num_k = 8
#num_k = 15
pam_res <- cluster::pam(prot_trajs_gam, num_k, diss = F)
cl <- pam_res$clustering
plot_clusters(cl, method, num_k, out_path = paste0(out_basedir,"trajectories/", method, "_k", num_k, "_131signif.pdf"), prot_trajs = prot_trajs_gam)
write.table(as.data.frame(cl), file = paste0(out_basedir, "trajectories/", method, "_k", num_k, "_131signif.txt"), quote = F, sep = "\t")







################################################################################
# GAM protein - phenotype
################################################################################

gam_res <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol = 8))
colnames(gam_res) <- c("prot", "pheno", "edf", "fval", "pval",  "p_nonlinear","N", "N_unique")

gam_res_no_tp <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno) -4), ncol = 8))
colnames(gam_res_no_tp) <- c("prot", "pheno", "edf", "fval", "pval", "p_nonlinear",  "N", "N_unique")

cnt <- 1
for (ph in colnames(pheno)[4:ncol(pheno)]) {
  cat(ph, "\n")
  for (prot in colnames(d_wide)[4:ncol(d_wide)]){
    res <- gam_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'none')
    gam_res_no_tp[cnt,] <- c(prot, ph, unlist(res))
    
    res2 <- gam_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'spline')
    gam_res[cnt,] <- c(prot, ph, unlist(res2))
    
     cnt <- cnt + 1
  }
}

gam_res <- gam_res %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
gam_res$BH_pval <- p.adjust(gam_res$pval)
gam_res <- gam_res[order(gam_res$pval),]
nrow(gam_res[gam_res$BH_pval < 0.05,])
write.table(gam_res, file = "../results/pheno_prot/prot_vs_pheno_splineTP_gam_adj_covar_no_outliers.txt", quote = F, sep = "\t", row.names = FALSE)

gam_res_no_tp <- gam_res_no_tp %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
gam_res_no_tp$BH_pval <- p.adjust(gam_res_no_tp$pval)
gam_res_no_tp <- gam_res_no_tp[order(gam_res_no_tp$pval),]
nrow(gam_res_no_tp[gam_res_no_tp$BH_pval < 0.05,])
write.table(gam_res_no_tp, file = "../results/pheno_prot/prot_vs_pheno_noTP_gam_adj_covar_no_outliers.txt", quote = F, sep = "\t", row.names = FALSE)


lmm_res <- read.delim("../results/pheno_prot/prot_vs_pheno_withTP_lmm_adj_covar.txt", sep = "\t", as.is = T, check.names = T)
lmm_res_no_tp <- read.delim("../results/pheno_prot/prot_vs_pheno_noTP_lmm_adj_covar.txt", sep = "\t", as.is = T, check.names = T)

tmp <- full_join(lmm_res, gam_res, by = c("prot", "pheno"))
tmp$logp.x <- -log10(tmp$pval.x)
tmp$logp.y <- -log10(tmp$pval.y)

bh.x <- min(tmp[tmp$BH_pval.x < 0.05, "logp.x"])
bh.y <- min(tmp[tmp$BH_pval.y < 0.05, "logp.y"], na.rm = T)

tmp[tmp$logp.y == 'Inf', "logp.y"] <- 20
ggplot(tmp, aes(logp.x, logp.y)) + geom_point() + theme_minimal() + 
  geom_hline(yintercept = bh.y) + geom_vline(xintercept = bh.x) +
  xlab("LMM") + ylab("GAM")


#### Functions #####

gam_prot_tp_adj_covar <- function(d_wide, prot, covariates, scale = F, rm_outliers = F){
  d_subs <- inner_join(d_wide[,c(prot, "ID", "TP")], covariates, by = c("ID"))
  colnames(d_subs)[1] <- "prot"
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs$ID <- as.factor(d_subs$ID)
  d_subs <- na.omit(d_subs)
  
  if (rm_outliers) d_subs <- remove_outliers_zscore(d_subs, "prot")
  
  fo_gam <- as.formula(paste("prot ~ s(TP, k = 4) + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+")))
  fo_gam_linear <- as.formula(paste("prot ~ TP + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+")))
  
  model <- gam(fo_gam, data = d_subs, method = 'ML')
  pval <- summary(model)$s.table["s(TP)","p-value"]
  
  model0 <- gam(fo_gam_linear, data = d_subs, method = 'ML')
  an <- anova(model, model0, test = 'LRT')
  p_nonlinear <- an$`Pr(>Chi)`[2]
    
  covar_means <- as.data.frame(lapply(covariates[,-1], function(x) if(is.numeric(x)) mean(x, na.rm = TRUE) else as.factor(2)))
  
  new_data <- expand.grid(
    TP = seq(1, 4, length.out = 100),
    ID = unique(d_subs$ID),
    predicted = NA
  ) %>%
    bind_cols(
      covar_means %>%
        slice(1)   # Ensure you use the first row of covar_means
    )
  new_data$predicted <- predict(model, newdata = new_data,  exclude = "s(ID)")
  new_data2 <- unique(new_data[,c("TP", "predicted")])
  
  return(list(pval = pval,  p_nonlinear = p_nonlinear, predicted = new_data2$predicted))
}

gam_pheno_prot_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F, adjust_timepoint = "spline"){
  
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
  
  d_subs <- remove_outliers_zscore(d_subs, "prot")
  
  if (adjust_timepoint == 'spline'){
    fo_lmm_txt <- paste("prot ~ s(pheno) + s(TP, k = 4) + s(ID,  bs = 're') +", paste(colnames(covariates)[-1], collapse = "+"))
  } else if (adjust_timepoint == 'linear') {
    fo_lmm_txt <- paste("prot ~ s(pheno) + TP + s(ID,  bs = 're') +", paste(colnames(covariates)[-1], collapse = "+"))
  } else if (adjust_timepoint == 'none') {
    fo_lmm_txt <- paste("prot ~ s(pheno) + s(ID,  bs = 're') +", paste(colnames(covariates)[-1], collapse = "+"))
  } else {
    stop ("Wrong adjust_timepoint argument. Should be one of spline, linear or none.")
  }
  
  fo_lmm_linear <- gsub("s\\(pheno\\)","pheno", fo_lmm_txt)
  model <- gam(as.formula(fo_lmm_txt), data = d_subs, method = 'ML')
  model0 <- gam(as.formula(fo_lmm_linear), data = d_subs, method = 'ML')
  
  coefs <- summary(model)$s.table
  
  an <- anova(model, model0, test = 'LRT')
  p_nonlinear <- an$`Pr(>Chi)`[2]
  
  #covar_means <- as.data.frame(lapply(covariates[,-1], function(x) if(is.numeric(x)) mean(x, na.rm = TRUE) else as.factor(2)))
  #new_data <- rbind(
  #  cbind(data.frame(TP = 1,  pheno = seq(min(d_subs$pheno), max(d_subs$pheno), length.out = 100), predicted = NA), covar_means),
  #  cbind(data.frame(TP = 2,  pheno = seq(min(d_subs$pheno), max(d_subs$pheno), length.out = 100), predicted = NA), covar_means),
  #  cbind(data.frame(TP = 3,  pheno = seq(min(d_subs$pheno), max(d_subs$pheno), length.out = 100), predicted = NA), covar_means),
  #  cbind(data.frame(TP = 4,  pheno = seq(min(d_subs$pheno), max(d_subs$pheno), length.out = 100), predicted = NA), covar_means))
  
  #new_data$predicted <- predict(model, newdata = new_data, exclude = "s(ID, bs = 're')")
  
  return(list(edf = coefs["s(pheno)", "edf"], fval = coefs["s(pheno)", "F"], pval = coefs["s(pheno)", "p-value"], p_nonlinear = p_nonlinear, 
              n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}


remove_outliers_zscore <- function(d, col = "prot", threshold = 3){
  d$zscore <- scale(d[, col])
  data_no_outliers <- d[abs(d$zscore) <= threshold, ]
  return(data_no_outliers)
}
