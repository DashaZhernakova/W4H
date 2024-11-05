library(dplyr)
library(ggplot2)

library(tidyr)
library(pheatmap)
library(patchwork)
library(tibble)
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")

setwd("/Users/Dasha/work/Sardinia/W4H/olink/data/")
PANEL = "INF"
d <- read.delim("INF_Q-14695_NPX_2024-09-28.txt", as.is = T, check.names = F, sep = "\t")

PANEL = "CVD"
d <- read.delim("CVD_Q-08150_NPX_2024-09-28.txt", as.is = T, check.names = F, sep = "\t")

d <- d[,-1]
# Extension controls are used to generate NPX values, remove
d <- d[d$AssayType != "ext_ctrl",]
#ggplot(d, aes(x = SampleID, y = NPX, col = SampleType)) + geom_point()

# Remove all technical entries
d_clean <- d[d$SampleType == "SAMPLE" & d$AssayType == "assay",]
d_clean$Timepoint <- gsub(".*_", "", d_clean$SampleID)
d_clean$ID <- gsub("_.*", "", d_clean$SampleID)

cat("Removing all technical samples and assays.\n")
cat("Total number of samples:", length(unique(d_clean$SampleID)), "\n\tNumber of individuals:", length(unique(d_clean$ID)), "\n\tNumber of proteins: ", length(unique(d_clean$Assay)), "\n")

cat("Proteins that completely failed QC:\n")
print(unique(d_clean[d_clean$Normalization == "EXCLUDED","Assay"]))

d_clean <- d_clean[d_clean$Normalization != "EXCLUDED",]


# Plate effects
prots <- sample(unique(d_clean$Assay), 6)
p <- ggplot(d_clean[d_clean$Assay %in% prots,], aes(x = Assay, y = NPX, color = PlateID)) +
  geom_boxplot() +
  scale_color_manual(values = my_colors) + 
  theme_bw()

p2 <- ggplot(d_clean[d_clean$Assay %in% prots,], aes(x = Assay, y = PCNormalizedNPX, color = PlateID)) +
  geom_boxplot() +
  scale_color_manual(values = my_colors) + 
  theme_bw()

pdf(paste0("../plots/", PANEL, ".plate_effects1.pdf"), width = 10, height = 5)
print(p + p2)
dev.off()




#Some proteins' measurements failed for specific samples

count_per_sample <- d_clean %>%
  group_by(ID,Timepoint, SampleQC) %>%
  summarise(sample_count = n_distinct(Assay)) %>%
  ungroup()

p <- ggplot(count_per_sample[count_per_sample$SampleQC == "PASS",], aes(x = ID, col = Timepoint, y = sample_count)) + 
  geom_point() +
  labs(x = "Sample", y = "Number of proteins ", fill = "Timepoint") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1.5, hjust=1, size = 0.5)) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  ggtitle("Number of proteins that pass QC per sample")

id_tmp <- count_per_sample[count_per_sample$SampleQC == "PASS" & count_per_sample$sample_count < 200, "ID"]
p <- p + annotate("text",x=Inf,y=-Inf,hjust=1,vjust=-0.5, color = "red", label = paste0("sample id: ", id_tmp))
pdf(paste0("../plots/", PANEL, ".tps_per_sample.pdf"))
print(p)
dev.off()

#The heatmap shows the number of timepoints measured for each sample and protein. some samples don't have all timepoints

count_non_na <- d_clean %>%
  group_by(Assay, ID) %>%
  summarise(non_na_count = sum(!is.na(NPX)), .groups = 'drop') %>%
  pivot_wider(names_from = ID, values_from = non_na_count, values_fill = 0)

count_non_na <- as.data.frame(count_non_na)
row.names(count_non_na) <- count_non_na$Assay 
count_non_na<- count_non_na[,-1]

#p <- pheatmap(as.matrix(count_non_na), cluster_rows = F, cluster_cols = F, fontsize = 6, color = c("#d1e9f9", "#80bae0", "#4a9ab0","#3b6988"), breaks = c(0,1,2,3,4))

#pdf("/Users/Dasha/work/Sardinia/W4H/olink//plots/heatmap_num_tp.pdf", width = 8, height = 15)
pheatmap(as.matrix(count_non_na), cluster_rows = F, cluster_cols = F, fontsize = 6, color = c("#d1e9f9", "#80bae0", "#4a9ab0","#3b6988"), breaks = c(0,1,2,3,4))
#dev.off()

# Boxplots for all samples

pdf(paste0("/Users/Dasha/work/Sardinia/W4H/olink//plots/", PANEL, ".boxplots_horiz.pdf"), width = 5, height = 35)
ggplot(d_clean,aes(y = SampleID, x = NPX, fill = Timepoint)) + 
  geom_boxplot() +
  scale_fill_manual(values = my_colors) + 
  theme_bw()
dev.off()

pdf(paste0("/Users/Dasha/work/Sardinia/W4H/olink//plots/", PANEL, ".boxplots_prot.pdf"), width = 5, height = 35)
ggplot(d_clean,aes(y = Assay, x = NPX)) + 
  geom_boxplot() +
  theme_bw()
dev.off()



#protein distributions look like.
#Doesn't depend on the timepoint. However some proteins seem to have a bi-modal distribution

p1 <- ggplot(d_clean, aes(x = NPX, fill = Timepoint, color = Timepoint)) +
  geom_histogram(alpha = 0.4, position = "identity", bins = 70) +
  scale_fill_manual(values = my_colors) + 
  theme_bw() +
  ggtitle("All samples and proteins combined together")


prots <- sample(unique(d_clean$Assay), 6)
p2 <- ggplot(d_clean[d_clean$Assay %in% prots,], aes(x = NPX)) +
  geom_histogram( bins = 70) +
  scale_fill_manual(values = my_colors) + 
  theme_bw() +
  facet_grid(cols = vars(Assay)) +
  ggtitle("For random proteins set 1")

prots <- sample(unique(d_clean$Assay), 6)
p3 <- ggplot(d_clean[d_clean$Assay %in% prots,], aes(x = NPX)) +
  geom_histogram( bins = 70) +
  scale_fill_manual(values = my_colors) + 
  theme_bw() +
  facet_grid(cols = vars(Assay)) +
  ggtitle("For random proteins set 2")

prots <- sample(unique(d_clean$Assay), 6)
p4 <- ggplot(d_clean[d_clean$Assay %in% prots,], aes(x = NPX)) +
  geom_histogram( bins = 70) +
  scale_fill_manual(values = my_colors) + 
  theme_bw() +
  facet_grid(cols = vars(Assay)) +
  ggtitle("For random proteins set 3")

prots <- sample(unique(d_clean$Assay), 6)
p6 <- ggplot(d_clean[d_clean$Assay %in% prots,], aes(x = NPX)) +
  geom_histogram( bins = 70) +
  scale_fill_manual(values = my_colors) + 
  theme_bw() +
  facet_grid(cols = vars(Assay)) +
  ggtitle("For random proteins set 4")

if (PANEL == "INF") {
  prots <- c("CST7", "GBP2", "IL4", "MICB_MICA", "PNLIPRP2")
} else {
  prots <- c("TNNI3", "TCL1B", "SERPINB5", "PM20D1", "CHIT1")
}
p5 <- ggplot(d_clean[d_clean$Assay %in% prots,], aes(x = NPX)) +
  geom_histogram( bins = 70) +
  scale_fill_manual(values = my_colors) + 
  theme_bw() +
  facet_grid(cols = vars(Assay)) +
  ggtitle("Bi-modal distribution")

pdf(paste0("../plots/", PANEL, ".distributions1.pdf"), width = 4, height = 4)
(p1)
dev.off()

pdf(paste0("../plots/", PANEL, ".distributions.pdf"), width = 15, height = 15)
print(p2 / p3 / p4 / p6 / p5)
dev.off()




# PCA
d_wide <- d_clean[,c("Assay","SampleID", "NPX")] %>%
  pivot_wider(names_from = Assay, values_from = NPX, values_fill = NA) %>%
  drop_na()

d_wide <- as.data.frame(d_wide)
row.names(d_wide) <- d_wide$SampleID
d_wide <- d_wide[, -1]

pca <- prcomp(d_wide, scale = T)
pca_to_plot <- as.data.frame(pca$x[,1:4])

pca_to_plot$Timepoint <- gsub(".*_", "", row.names(pca_to_plot))


g1 <- ggplot(pca_to_plot, aes(x = PC1, y = PC2, col = Timepoint)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = my_colors) +
  theme(legend.position="none")

g2 <- ggplot(pca_to_plot, aes(x = PC2, y = PC3, col = Timepoint)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = my_colors) + 
  theme(legend.position="none")

g3 <- ggplot(pca_to_plot, aes(x = PC3, y = PC4, col = Timepoint)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = my_colors)

pdf(paste0("/Users/Dasha/work/Sardinia/W4H/olink//plots/", PANEL, "_pca.pdf"), width = 10, height = 5)
print(g1 + g2 + g3)
dev.off()


# Association of PCs with phenotypes
pheno <- read.delim("/Users/Dasha/work/Sardinia/W4H/olink/data/pheno_subset.fmt.txt", check.names = F, as.is = T, sep = "\t")
pheno$id <- sprintf("%03d", pheno$id)
pca_to_plot <- pca_to_plot %>%
  rownames_to_column( var = "SampleID")
pca_to_plot$SampleID <- gsub("_.*", "", pca_to_plot$SampleID)
pca_to_plot_pheno <- left_join(pca_to_plot[pca_to_plot$Timepoint == "1",], pheno, by = c("SampleID" = "id"))
pca_to_plot_pheno <- pca_to_plot_pheno[pca_to_plot_pheno$Timepoint == "1",]

ggplot(pca_to_plot_pheno, aes(x = PC1, y = PC2, col = age )) + 
  geom_point() +
  theme_bw() 


# Association of proteins with phenotypes

d_wide_1 <- d_wide[grepl(".*_1", row.names(d_wide)),]
row.names(d_wide_1) <- gsub("_1", "", row.names(d_wide_1))
d_wide_1 <- d_wide_1 %>%
  rownames_to_column( var = "SampleID")

joined_data <- inner_join(d_wide_1, pheno, by = c("SampleID" = "id"))

regression_olink_pheno <- function(joined_data, prot, ph){
  d <- na.omit(joined_data[,c("SampleID", prot, ph)])
  colnames(d) <- c("SampleID", "prot", "pheno")
  is_factor <- length(unique(d$pheno)) < 3
  
  if (is_factor){
    d$pheno <- as.factor(d$pheno)
  } 
  
  lm_fit <- lm(prot ~ pheno, data = d)
  coef <- summary(lm_fit)$coefficients
  return(list("pval" = coef[2,4], "est" = coef[2,1], "n" = nrow(d)))
}

correl_res <- data.frame(matrix(nrow = (ncol(d_wide_1) -1) * (ncol(pheno) -1), ncol = 5))
colnames(correl_res) <- c("prot", "pheno", "pval", "r", "N")
cnt <- 1
for (prot in colnames(d_wide_1)[-1]){
  for (ph in colnames(pheno)[-1]){

    #res <- correlate_olink_pheno(joined_data, prot, ph)
    res <- regression_olink_pheno(joined_data, prot, ph)
    correl_res[cnt,] <- c(prot, ph, unlist(res))
    cnt <- cnt + 1
  }
}
correl_res$pval <- as.numeric(correl_res$pval)
correl_res <- correl_res[order(correl_res$pval),]
correl_res$p_adjBH <- p.adjust(correl_res$pval, method = "BH")
write.table(correl_res, file = paste0(PANEL, "_lm_phenotypes.txt"), quote = F, col.names = NA, sep = "\t")


make_pheno_prot_plot <- function(joined_data, correl_res_line){
  
  prot <- correl_res_line$prot
  ph <- correl_res_line$pheno
  
  d <- na.omit(joined_data[,c("SampleID", prot, ph)])
  is_factor <- length(unique(d[,ph])) < 3
  
  
  if(! is_factor){
    
    
    prot <- ensym(prot)
    ph <- ensym(ph)
    p <- ggplot(d, aes(x = !!ph, y = !!prot)) + 
      geom_point() +
      geom_smooth(method = 'lm') +
      theme_bw() + 
      ggtitle(paste0("P = ", formatC(as.numeric(correl_res_line$pval), format = 'g', digits = 2), "; r = ", formatC(as.numeric(correl_res_line$r), format = 'f', digits = 2) )) 
    
  } else {
    
    d[,ph] <- as.factor(d[,ph])
    prot <- ensym(prot)
    ph <- ensym(ph)
    p <- ggplot(d, aes(x = !!ph, y = !!prot, group = !!ph)) + 
      geom_boxplot() + 
      theme_bw() + 
      ggtitle(paste0("P = ", formatC(as.numeric(correl_res_line$pval), format = 'g', digits = 2), "; r = ", formatC(as.numeric(correl_res_line$r), format = 'f', digits = 2) )) 
    
  }
  return(p)
}

plot_list <- list()
cnt <- 1
for (i in 1:nrow(correl_res[correl_res$pval < 0.001,]) ){
  plot_list[[cnt]] <- make_pheno_prot_plot(joined_data, correl_res[i, ])
  cnt <- cnt + 1
}

library(gridExtra)
pdf(paste0("../plots/", PANEL, "_corr_with_pheno.pdf"), width = 20, height = 20)
grid.arrange(grobs = plot_list, ncol = 5, nrow = 4)  # 5 columns, 3 rows for 13 plots
dev.off()


##### LMM
library(lme4)

lmm_res <- data.frame(matrix(nrow = length(unique(d_clean$Assay)), ncol = 2))
row.names(lmm_res) <- unique(d_clean$Assay)
colnames(lmm_res) <- c("lmer_pval", "BH_qval")
for (prot in unique(d_clean$Assay)){
  d_clean$Timepoint <- as.factor(d_clean$Timepoint)
  model <- lmer(NPX ~ Timepoint + (1|ID), data = d_clean[d_clean$Assay == prot,])
  reduced_model <- lmer(NPX ~ (1|ID), data = d_clean[d_clean$Assay == prot,])
  res_p <- anova(reduced_model, model)["model","Pr(>Chisq)"]
  lmm_res[prot,1] <- res_p
}

lmm_res <- lmm_res[order(lmm_res$lmer_pval),]
lmm_res$BH_qval <- p.adjust(lmm_res$lmer_pval)

signif <- lmm_res[lmm_res$BH_qval < 0.05,]
plot_list <- list()

for (cnt in 1:nrow(signif) ){
  prot <- row.names(signif)[cnt]
  plot_list[[cnt]] <- ggplot(d_clean[d_clean$Assay == prot,], aes (y = NPX, x = Timepoint)) + 
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.3) + 
    theme_bw() +
    ggtitle(paste0(prot, ", P = ", formatC(signif[cnt, "lmer_pval"], digits = 3)))
}
pdf(paste0("../plots/", PANEL, "_prot_tp_lmm.pdf"), width = 20, height = 20)
grid.arrange(grobs = plot_list, ncol = 5, nrow = 5)  # 5 columns, 3 rows for 13 plots
dev.off()

# non significant
nonsignif <- lmm_res[lmm_res$lmer_pval > 0.05,]
plot_list <- list()

for (cnt in 1:nrow(nonsignif) ){
  prot <- row.names(nonsignif)[cnt]
  plot_list[[cnt]] <- ggplot(d_clean[d_clean$Assay == prot,], aes (y = NPX, x = Timepoint)) + 
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.3) + 
    theme_bw() +
    ggtitle(paste0(prot, ", P = ", formatC(nonsignif[cnt, "lmer_pval"], digits = 3)))
  
  if (cnt > 15) break
}
pdf(paste0("../plots/", PANEL, "_prot_tp_lmm_NONSIGNIF.pdf"), width = 20, height = 20)
grid.arrange(grobs = plot_list, ncol = 5, nrow = 5)  # 5 columns, 3 rows for 13 plots
dev.off()


#
# Compare overlapping between panels and combine
#

inf <- read.delim("INF_Q-14695_NPX_2024-09-28.txt", as.is = T, check.names = F, sep = "\t")

cvd <- read.delim("CVD_Q-08150_NPX_2024-09-28.txt", as.is = T, check.names = F, sep = "\t")

filter_table <- function(d){
  d <- d[,-1]
  d <- d[d$AssayType != "ext_ctrl",]
  d_clean <- d[d$SampleType == "SAMPLE" & d$AssayType == "assay",]
  d_clean$Timepoint <- gsub(".*_", "", d_clean$SampleID)
  d_clean$ID <- gsub("_.*", "", d_clean$SampleID)
  d_clean <- d_clean[d_clean$Normalization != "EXCLUDED",]
  cat("Total number of samples:", length(unique(d_clean$SampleID)), "\n\tNumber of individuals:", length(unique(d_clean$ID)), "\n\tNumber of proteins: ", length(unique(d_clean$Assay)), "\n")
  return(d_clean)
}

inf2 <- filter_table(inf)
cvd2 <- filter_table(cvd)

#Remove 100_2 - the sample that failed CVD completely and partially INF
inf2 <- inf2[inf2$SampleID != '100_2',]
cvd2 <- cvd2[cvd2$SampleID != '100_2',]

# Merge 2 panels
merged_long <- rbind(inf2[,c("SampleID", "Assay", "NPX", "Panel")], cvd2[,c("SampleID", "Assay", "NPX", "Panel")])
merged_long$Assay2 <- paste0(merged_long$Assay, "_", merged_long$Panel)
merged_wide <-merged_long[,c("SampleID", "Assay2", "NPX")] %>%
  pivot_wider(names_from = Assay2, values_from = NPX)



# Plot a heatmap of NPX values of the combined table
tmp_merged_wide <- as.data.frame(merged_wide[,-1])
row.names(tmp_merged_wide) <- merged_wide$SampleID

prot_panel_annot <- as.data.frame(gsub(".*_", "", colnames(tmp_merged_wide)))
row.names(prot_panel_annot) <- colnames(tmp_merged_wide)
colnames(prot_panel_annot) <- "Panel"

ind_panel_annot <- as.data.frame(gsub(".*_", "", rownames(tmp_merged_wide)))
row.names(ind_panel_annot) <- row.names(tmp_merged_wide)
colnames(ind_panel_annot) <- "Timepoint"


pheatmap(as.matrix(tmp_merged_wide), annotation_col = prot_panel_annot, annotation_row = ind_panel_annot,
         show_rownames = F, show_colnames = F, filename = "../plots/merged_panels_heatmap.pdf")


# Plot overlapping proteins
overlap <- intersect(unique(inf2$Assay), unique(cvd$Assay))


cmp_overlapping <- function(merged_wide, p){
  tmp_d <- merged_wide[,c(paste0(p, "_Inflammation"), paste0(p, "_Cardiometabolic"))]
  colnames(tmp_d) <- c("INF", "CVD")
  
  # counts
  cat("INF:\n")
  print(table(is.na(tmp_d$INF)))
  cat("CVD:\n")
  print(table(is.na(tmp_d$CVD)))
  
  
  cor <- cor(tmp_d$INF, tmp_d$CVD, method = "pearson", use = "complete.obs")
  
  p <- ggplot(tmp_d, aes(x = INF, y = CVD)) + 
    geom_point() +
    geom_smooth(method = "lm") +
    ggtitle(paste0("Protein: ", p, "; Pearson r = ", formatC(cor, digits = 3))) +
    theme_bw()

  return(p)
}

plot_list <- list()

for (p in overlap){
  plot_list[[p]] <- cmp_overlapping(merged_wide, p)
}

pdf("../plots/cmp_overlapping.pdf", width = 15, height = 5)
grid.arrange(grobs = plot_list, ncol = 3, nrow = 1)  
dev.off()


# Combine clean without panel name
# For the overlapping proteins keep the INF version as default and add a suffix to the CVD version:
merged_long2 <- merged_long
merged_long2[merged_long2$Assay2 %in% c(paste0(overlap, "_Cardiometabolic")), "Assay"] <- merged_long2[merged_long2$Assay2 %in% c(paste0(overlap, "_Cardiometabolic")), "Assay2"]
merged_long2$Assay2 <- NULL
merged_long2$Panel <- NULL
merged_long2$Timepoint <- gsub(".*_","",merged_long2$SampleID)
merged_long2$ID <- gsub("_.*","",merged_long2$SampleID)
merged_long2Assay <- gsub("_Cardiometabolic", "_CVD", merged_long2$Assay)
merged_wide2 <-merged_long2[,c("SampleID", "Assay", "NPX")] %>%
  pivot_wider(names_from = Assay, values_from = NPX)

num_cols_with_na <- sum(colSums(is.na(merged_wide2)) > 0)
num_rows_with_na <- sum(rowSums(is.na(merged_wide2)) > 0)

write.table(merged_wide2, file = "olink_clean_CVD+INF.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(merged_long2, file = "olink_clean_CVD+INF_long.txt", quote = F, sep = "\t", row.names = FALSE)


# PCA
tmp_merged_wide2 <- as.data.frame(merged_wide2[,-1])
row.names(tmp_merged_wide2) <- merged_wide2$SampleID
tmp_merged_wide2 <- na.omit(tmp_merged_wide2)

pca <- prcomp(tmp_merged_wide2, scale = T)
pca_to_plot <- as.data.frame(pca$x[,1:4])

pca_to_plot$Timepoint <- gsub(".*_", "", row.names(pca_to_plot))


g1 <- ggplot(pca_to_plot, aes(x = PC1, y = PC2, col = Timepoint)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = my_colors) +
  theme(legend.position="none")

g2 <- ggplot(pca_to_plot, aes(x = PC2, y = PC3, col = Timepoint)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = my_colors) + 
  theme(legend.position="none")

g3 <- ggplot(pca_to_plot, aes(x = PC3, y = PC4, col = Timepoint)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = my_colors)

pdf(paste0("/Users/Dasha/work/Sardinia/W4H/olink//plots/merged_panels.pca.pdf"), width = 10, height = 5)
print(g1 + g2 + g3)
dev.off()

