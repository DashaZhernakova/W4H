library(timeOmics)
library(netOmics)
library(tidyverse)
library(lmms)

setwd("/Users/Dasha/work/Sardinia/W4H/olink/data/")

d_wide <- read.delim("olink_clean_CVD+INF.txt", as.is = T, check.names = F, sep = "\t")
ID <- gsub("_.*", "", d_wide$SampleID)
TP <- gsub(".*_", "", d_wide$SampleID)
d_wide <- cbind(ID, TP, d_wide)

pheno <- read.delim("../../phenotypes/blood_pheno_03102024_log.txt", as.is = T, check.names = F, sep = "\t")
pheno$Record.ID <- gsub("ID_", "", pheno$Record.ID)
colnames(pheno) <- gsub("Record.ID","ID",colnames(pheno))
pheno <- cbind(SampleID = paste0(pheno$ID, "_", pheno$TP), pheno)
pheno$Age <- NULL

lmm_res_prot_tp <- read.delim("../results/prot_vs_tp_poly3_lmm_adj_age_bmi_preg_storage.txt", sep = "\t", as.is = T, check.names = F)
lmm_res_pheno_tp <- read.delim("../results/pheno_prot/pheno_vs_tp_poly3_lmm_adj_age_bmi_preg_storage.txt", sep = "\t", as.is = T, check.names = F)
sign_prots <- lmm_res_prot_tp[lmm_res_prot_tp$BH_pval < 0.05,"prot"]
sign_pheno <- lmm_res_pheno_tp[lmm_res_pheno_tp$BH_pval < 0.05,"pheno"]


prot_d <- d_wide[,sign_prots]
row.names(prot_d) <- d_wide$SampleID
prot_d <- na.omit(prot_d)

pheno_d <- pheno[,sign_pheno]
row.names(pheno_d) <- pheno$SampleID
pheno_d <- na.omit(pheno_d)

shared_ids <- intersect(row.names(prot_d), row.names(pheno_d))
prot_d <- prot_d[shared_ids,]
pheno_d <- pheno_d[shared_ids,]

prot_d <- scale(prot_d)
pheno_d <- scale(pheno_d)

all(row.names(prot_d) == row.names(pheno_d))


#
# Model trajectories using lmmSpline
#

prot_profile <- get_lmm_spline_trajs(prot_d, n_points = 10)
pheno_profile <- get_lmm_spline_trajs(pheno_d, n_points = 10)

# Multi omics clustering
data.lmms <- list("pheno" = pheno_profile,
                  "prot" = prot_profile)
block_pls_clustering(data.lmms)


# Try our trajectories
prot_trajs = read.delim("../results/trajectories/protein_trajectories_109signif.txt", sep = "\t", as.is = T, check.names = F)
prot_trajs <- prot_trajs %>%
  as.data.frame() %>%
  {rownames(.) <- .[[1]]; .[-1]}
prot_profile <- as.data.frame(t(prot_trajs)) 

pheno_trajs = read.delim("../results/trajectories/pheno_traj_adj_covar_signif_pheno.txt", sep = "\t", as.is = T, check.names = F)

pheno_trajs <- pheno_trajs %>%
  as.data.frame() %>%
  {rownames(.) <- .[[1]]; .[-1]}
pheno_profile <- as.data.frame(t(pheno_trajs)) 

# Multi omics clustering
data.lmms <- list("pheno" = pheno_profile,
                  "prot" = prot_profile)
block_pls_clustering(data.lmms, n_pcs = 3)






block_pls_clustering <- function(data.lmms, do_scaling = F, n_pcs = 3){
  cb.lmms <- do.call("cbind", data.lmms)
  colnames(cb.lmms) <- lapply(data.lmms, colnames) %>% unlist
  
  block.pls.res <- block.pls(data.lmms, ncomp = n_pcs, indY = 1, 
                         scale = do_scaling, mode = "canonical")
  
  block.ncomp <- getNcomp(block.pls.res, X = data.lmms,  indY = 1,
                          scale = do_scaling, mode = "canonical")
  block.ncomp$choice.ncomp
  plot(block.ncomp)
  
  block.pls.res <- block.pls(data.lmms, ncomp = n_pcs, indY = 1, 
                             scale = do_scaling, mode = "canonical")
  
  block.pls.cluster <- getCluster(block.pls.res)
  
  plotLong(block.pls.res)
}


get_lmm_spline_trajs <- function(in_data, n_points = 4){
  time <- as.numeric(gsub(".*_", "", row.names(in_data)))
  lmms.output <- lmms::lmmSpline(data = in_data, time = time,
                                 sampleID = rownames(in_data), deri = FALSE,
                                 basis = "p-spline", numCores = 4, timePredict = seq(1,4, length.out = n_points),
                                 keepModels = TRUE)
  modelled.data <- t(slot(lmms.output, 'predSpline'))
  
  data.gathered <- modelled.data %>% as.data.frame() %>% 
    rownames_to_column("time") %>%
    mutate(time = as.numeric(time)) %>%
    pivot_longer(names_to="feature", values_to = 'value', -time)
  
  filter.res <- lmms.filter.lines(data = in_data, 
                                  lmms.obj = lmms.output, time = time)
  
  profile.filtered <- filter.res$filtered
  
  # plot profiles
  ggplot(data.gathered[data.gathered$feature %in% colnames(profile.filtered),], aes(x = time, y = value, color = feature)) + geom_line() +
    theme_bw() + ggtitle("`lmms` profiles") + ylab("Feature expression") +
    xlab("Time") + theme(legend.position = "none")
  
  return(profile.filtered)
}


#
# Cluster trajectories
#

cluster_trajs_pca <- function(profile.filtered){
  pca.res <- pca(X = profile.filtered, ncomp = 7, scale=FALSE, center=FALSE)
  pca.ncomp <- getNcomp(pca.res, max.ncomp = 7, X = profile.filtered, 
                        scale = FALSE, center=FALSE)
  plot(pca.ncomp)
  pca.ncomp$choice.ncomp
  pca.cluster <- getCluster(pca.res)
  head(pca.cluster)
  
  plotLong(pca.res, scale = FALSE, center = FALSE, 
           title = "PCA longitudinal clustering")
  
  # sparse PCA (?)
  tune.spca.res <- tuneCluster.spca(X = profile.filtered, ncomp = 3, 
                                    test.keepX = c(2:10))
  
  #plot(tune.spca.res)
  spca.res <- spca(X = profile.filtered, ncomp = 3, 
                   keepX = tune.spca.res$choice.keepX, scale = FALSE)
  plotLong(spca.res)
  return(list("pca" = pca.res, "spca" = spca.res))
}

