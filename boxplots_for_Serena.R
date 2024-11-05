my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/data/")
source("../scripts/traj_functions.R")

library(ggplot2)
library("readxl")
library(dplyr)
library(lme4)
library(grid)
library(gridExtra)
library(RColorBrewer)

d_wide <- read.delim("olink_clean_CVD+INF.txt", as.is = T, check.names = F, sep = "\t")
d_long <- read.delim("olink_clean_CVD+INF_long.txt", as.is = T, check.names = F, sep = "\t")

ID <- gsub("_.*", "", d_wide$SampleID)
TP <- gsub(".*_", "", d_wide$SampleID)

d_wide <- cbind(ID, TP, d_wide)

olink_panels <- as.data.frame(read_excel("olink_protein_names.xlsx"))
inf_prots <- intersect(colnames(d_wide), olink_panels[olink_panels$Panel == 'INF', "Assay"])
d_wide_inf <- d_wide[,c("ID", "TP", "SampleID",  inf_prots)]

lmm_res <- data.frame(matrix(nrow = (ncol(d_wide_inf) -4), ncol = 2))
colnames(lmm_res) <- c("prot", "pval")
cnt <- 1
for (prot in colnames(d_wide_inf)[4:ncol(d_wide_inf)]){
  res <- lmm_prot_tp_poly3(d_wide_inf, prot)
  lmm_res[cnt,] <- c(prot, res)
  cnt <- cnt + 1
}

lmm_res <- na.omit(lmm_res) %>%
  mutate(across(-c( prot), as.numeric)) 

lmm_res$BH_pval <- p.adjust(lmm_res$pval, method = 'BH')
lmm_res <- lmm_res[order(lmm_res$pval),]

signif <- lmm_res[lmm_res$BH_pval < 0.05,]
nrow(signif)


colors <- c("F" = "#1e9f77", "O" = "#d95e00", "eL" = "#7570b3", "lL" = "#e7298a")

d_wide_inf_long <- d_wide_inf %>%
  pivot_longer(cols = 4:ncol(d_wide_inf), names_to = 'prot', )

d_wide_inf_long$Phase <- NA
d_wide_inf_long[d_wide_inf_long$TP == 1, "Phase"] <- "F"
d_wide_inf_long[d_wide_inf_long$TP == 2, "Phase"] <- "O"
d_wide_inf_long[d_wide_inf_long$TP == 3, "Phase"] <- "eL"
d_wide_inf_long[d_wide_inf_long$TP == 4, "Phase"] <- "lL"

d_wide_inf$Phase <- factor(d_wide_inf$Phase, ordered = T, levels = c('F', 'O', 'eL', 'lL'))

subs <- d_wide_inf_long[d_wide_inf_long$prot %in% head(lmm_res, 12)$prot,]
subs$TP <- as.factor(subs$TP)


pdf("../plots/INF_prot_vs_TP_LMM.pdf")
ggplot(subs, aes(x = TP, y = value, fill = TP, group = TP)) + 
  geom_boxplot() +  # Customize outliers
  scale_x_discrete(labels=c("F", "O", "eL", "lL")) +
  scale_fill_brewer(palette= 'Dark2', labels=c("Follicular", "Ovulation", "early Luteal", "late Luteal")) +  # Apply custom colors
  theme_classic() + labs(fill="Time points") + theme(strip.background =element_rect(fill="lightgrey")) +
  facet_wrap(~prot, scales ="free")

dev.off()
