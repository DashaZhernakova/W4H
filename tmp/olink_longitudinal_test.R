library(ggplot2)
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/data/")

d <- read.delim("olink_clean_CVD+INF_long.txt", as.is = T, check.names = F, sep = "\t")
pheno <- read.delim("pheno_subset.fmt.txt", check.names = F, as.is = T, sep = "\t")
pheno$id <- sprintf("%03d", pheno$id)

prots <- sample(unique(d$Assay), 100)
ggplot(d[d$Assay %in% prots,], aes(x = Timepoint, y = NPX, group  = Assay)) +  geom_smooth(size = 0.3,se = F) + theme_bw()


#
# SplinectomeR
#
library(splinectomeR)
prots <- sample(unique(d$Assay), 300)
plot_list <- list()
trends <- data.frame(matrix(nrow = 1000, ncol = length(prots)))
colnames(trends) <- prots
for (prot in prots){
  trend_per_prot <- trendyspliner(data = d[d$Assay == prot,], cases = 'ID', xvar = 'Timepoint', yvar = 'NPX',
                               perms = 99, mean_center = F)
  trends[,prot] <- trend_per_prot$imputed_curve[,"var1"]
  #pval <- trend_per_prot$pval
  #if (pval < 0.001) {
  #  plot_list[[prot]] <- trendyspliner.plot.perms(trend_per_prot, xlabel = 'time', ylabel = 'weight') + ggtitle(paste0(prot, ", pval = ", pval))
  #}
}

row.names(trends) <- trend_per_prot$imputed_curve[,"x"]
trends4 <- as.data.frame(t(trends[c(1.15, 2.001006, 3.001074, 3.857500),]))

trends <- as.data.frame(t(trends))
trends4 <- trends[,colnames(trends) %in% c(1.15, 2.001006, 3.001074, 3.857500)]
cldSDQ <- cld(trends4,timeInData=1:4,time=c(1,2,3,4))
kml(cldSDQ,nbRedrawing=3,toPlot="both")
  


#
# rm correlations
#
library(rmcorr)
library(dplyr)
library(lme4)
pheno <- read.delim("../../phenotypes/blood_pheno_03102024.txt", as.is = T, check.names = F, sep = "\t")
d_wide <- read.delim("olink_clean_CVD+INF.txt", as.is = T, check.names = F, sep = "\t")

d_wide$ID <- gsub("_.*", "", d_wide$SampleID)
d_wide$TP <- gsub(".*_", "", d_wide$SampleID)

pheno$Record.ID <- gsub("ID_", "", pheno$Record.ID)
pheno$SampleID <- paste0(pheno$Record.ID, "_", pheno$TP)

prot <- "REN"
ph <- "LDLC"

res_cor <- data.frame(matrix(ncol = 4, nrow = ncol(pheno) * ncol(d_wide)))
cnt <- 1
for (ph in colnames(pheno)[5:15]){
  print(ph)
  for (prot in colnames(d_wide)[2:(ncol(d_wide) - 2)]){
  
    d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
    colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")

    corr <- rmcorr(participant = ID, measure1 = prot, measure2 = pheno, dataset = d_subs)
    res_cor[cnt, ] <- c(ph, prot, corr$r, corr$p)
    cnt <- cnt + 1
  }
}
colnames(res_cor) <- c("pheno", "prot", "r", "pval")
res_cor$r <- as.numeric(res_cor$r)
res_cor$pval <- as.numeric(res_cor$pval)
res_cor <- res_cor[order(res_cor$pval),]


plot_together <- function(d_wide, pheno, prot, ph){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  g <- ggplot(d_subs, aes(x = TP, y = scale(prot))) +
    geom_smooth(color = "red", aes(x = TP, y = scale(pheno))) +  
    geom_smooth(color = "blue", aes(x = TP, y = scale(prot))) +  
    labs(x = "Timepoint ", y = prot, 
         title = paste0(ph, " - ", prot)) +
    theme_minimal()
  
  g
  
}

get_rmcorr <- function(d_wide, pheno, prot, ph){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  
  corr <- rmcorr(participant = ID, measure1 = prot, measure2 = pheno, dataset = d_subs)
  
  return(corr)
}

lmm_pheno_prot <- function(d_wide, pheno, prot, ph){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  model <- lmer(pheno ~ prot + TP + (1|ID), data = d_subs)
  est <- summary(model)$coefficients["prot", "Estimate"]
  model0 <- lmer(pheno ~  TP + (1|ID), data = d_subs)
  an <- anova(model, model0)
  pval <- an$`Pr(>Chisq)`[2]
  return(list(estimate = est, pval = pval))
}

res_cor$lmer_p <- NA
for (i in 1:nrow(res_cor)){
  res_cor[i, "lmer_p"] <- lmm_pheno_prot(d_wide, pheno, res_cor[i, "prot"], res_cor[i, "pheno"])
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

scatter_col_tp(d_wide, pheno, "RETN", "PROG")
scatter_col_tp(d_wide, pheno, "TNFSF10", "PROG")



#
# Cluster trajectories
#
library(splines)
prot <- "PROK1"

fit_poly3 <- function(d_wide, prot, n = 7, scale = F){
  d_subs <- d_wide[,c("ID", "TP", prot)]
  colnames(d_subs) <- c("ID", "TP", "prot")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  lm_fit <- lm(prot ~ poly(TP,3), data = d_subs)
  new_data <- data.frame(TP = seq(1,4, length.out = n), predicted = NA)
  new_data$predicted <- predict(lm_fit, newdata = new_data )
  if (scale) new_data$predicted <- scale(new_data$predicted)
  return (new_data)
}

fit_poly3_get_coef <- function(d_wide, prot){
  d_subs <- d_wide[,c("ID", "TP", prot)]
  colnames(d_subs) <- c("ID", "TP", "prot")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  d_subs$prot <- scale(d_subs$prot)
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

get_dtw <- function(d_wide, pheno, prot, ph){
  spline_prot <- fit_splines(d_wide, prot)
  spline_ph <- fit_splines(pheno, ph)
  
  dtw_dist <- dtw(scale(spline_prot$predicted), scale(spline_ph$predicted))$distance
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

ggplot(d_subs, aes(x = TP, y = prot)) +
  geom_point() +  # Plot the observed data
  geom_line(data = new_data, aes(x = TP, y = predicted), color = "blue") +  # Add the fitted line
  labs(title = "Fitted Polynomial Model", x = "Timepoint", y = "Protein levels")

# 0. Make the trajectories
prots <- sample(unique(d$Assay), 50)
n = 7
trajectories <- data.frame(matrix(nrow = length(prots), ncol = n))
cnt <- 1
for (prot in prots){
  trajectories[cnt,] <- fit_poly3(d_wide, prot, n)$predicted
  cnt <- cnt + 1
}
colnames(trajectories) <- seq(1,4, length.out = n)
trajectories <- scale(trajectories)

# 1. KLM
cld_data <- cld(trajectories,timeInData=1:7,time=seq(1,4, length.out = n))
### 2. Building "optimal" clusteration (with only 3 redrawings)
kml(cld_data,nbRedrawing=5,toPlot="both")
X11(type = "Xlib")
kml::choice(cld_data)

# 2. DTW
library(dtw)

dtw_results <- data.frame()

for (ph in colnames(pheno)[5:15]){
  print(ph)
  for (prot in colnames(d_wide)[2:(ncol(d_wide) - 2)]){
    d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
    colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
    d_subs$TP <- as.numeric(d_subs$TP)
    d_subs <- na.omit(d_subs)
    dtw_dist <- dtw(scale(d_subs$prot), scale(d_subs$pheno))$distance
    dtw_results <- rbind(dtw_results, data.frame(Protein = prot, Phenotype = ph, DTW_Distance = dtw_dist))
  }
}

# Lower DTW distance indicates more similar trajectories
dtw_results <- dtw_results %>% arrange(DTW_Distance)
plot_together(d_wide, pheno, "TP53INP1", "PROG")


# k-means

# on poly coef
fit_results <- data.frame(matrix(nrow = ncol(d_wide) - 2, ncol = 3))
cnt <- 1
for (prot in colnames(d_wide)[2:(ncol(d_wide) - 2)]){
  co <- fit_poly3_get_coef(d_wide, prot)
  fit_results[cnt,] <-  co
  row.names(fit_results)[cnt] <- prot
  cnt <- cnt + 1
}
fit_results <- na.omit(fit_results)

num_k = 5
km <- kmeans(fit_results, num_k, nstart = 50)
cl <- km$cluster

clusters <- as.data.frame(cl)

#ss <- silhouette(km$cluster, dist(cor(t(fit_results))))
#mean_sil_score <- aggregate(ss[,3]~ss[,1], FUN=mean)



#
#
library(gridExtra)
pheno$TP = as.numeric(pheno$TP)
plot_list <- list()
cnt <- 1
for (ph in colnames(pheno)[5:ncol(pheno)]){
  ph = ensym(ph)
  g <- ggplot(pheno, aes(y = !!ph, x = TP)) + geom_smooth() + theme_bw()
  g2 <- ggplot(pheno, aes(x = !!ph)) + geom_density() + theme_bw()
  plot_list[[cnt]] <- g
  plot_list[[cnt+1]] <- g2
  cnt <- cnt + 2
 
}
pdf("../plots/pheno_distr_log.pdf", width = 15, height = 25)
grid.arrange(grobs = plot_list, ncol = 4, nrow = 7)  
dev.off()


pheno <- read.delim("../../phenotypes/blood_pheno_03102024.txt", as.is = T, check.names = F, sep = "\t")
pheno <- pheno %>% select(-one_of(c("sesso", "TSH", "FT4", "TST")))
# Log-transform PROG, HOMA-IR, HOMA_B, INS_mIU_L, ALT_U_L, AST_U_L, Trigl
ph_to_log <- c("PROG", "HOMA_IR", "HOMA_B", "INS", "ALT", "AST", "Trigl", "LH", "FSH", "X17BES", "PRL")
pheno2 <- pheno
for (ph in ph_to_log){
  pheno2[,ph] <- log(pheno[,ph] + 1)
}
write.table(pheno2, file = "../../phenotypes/blood_pheno_03102024_log.txt", row.names = F, quote = F, sep = "\t")
