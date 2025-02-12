
mb_raw <- read.delim("../../mb_abundance/W4H_FECAL_merged_abundance_table.txt", sep = "\t", check.names = F, as.is = T, row.names = 1)

mb_raw <- as.data.frame(t(mb_raw))
row.names(mb_raw) <- gsub(".*-","",row.names(mb_raw))
row.names(mb_raw) <- gsub("_metaphlan","",row.names(mb_raw))

mb_raw$UNCLASSIFIED=NULL
mb_raw <- mb_raw[!grepl("Mock", row.names(mb_raw)), ]
colnames(mb_raw)=gsub(".*s__","",colnames(mb_raw)) # simplify species names

# Filter taxa with mean abundance > 0.5
mb_flt_abund <- mb_raw[, colMeans(mb_raw) > 0.5]
dim(mb_flt_abund)

# Calculate the proportion of non-zero values for each taxa
non_zeros <- colSums(mb_flt_abund != 0) / nrow(mb_flt_abund)

# Filter columns where the non-zero proportion is > 0.2
mb_flt_abund_prev <- mb_flt_abund[, non_zeros > 0.2]
dim(mb_flt_abund_prev)

mb_flt_sgb=mb_flt_abund_prev[,grep("t__",colnames(mb_flt_abund_prev))]
colnames(mb_flt_sgb)=gsub(".*s__","",colnames(mb_flt_sgb)) # simplify species names
dim(mb_flt_sgb)
mb_flt_sp=mb_flt_abund_prev[,grepl("s__", colnames(mb_flt_abund_prev)) & !grepl("t__", colnames(mb_flt_abund_prev))]
colnames(mb_flt_sgb)=gsub(".*s__","",colnames(mb_flt_sgb)) # simplify species names
dim(mb_flt_sp)





do_clr_externalWeighting = function(interest_matrix, core_matrix){
  interest_matrix = interest_matrix + min(core_matrix[core_matrix>0])/2
  core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  
  # estimate weighting parameter
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_core = apply(core_matrix, 1, gm_mean)
  
  #do transformation
  data_prepared = cbind(Gmean_core,interest_matrix)
  data_transformed = t(apply(data_prepared,1,function(x){
    log(x / x[1])[-1]
  }))
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}





species <- read.delim("../../mb_abundance/W4H_FECAL_species_clr_transformed_metaphlanAbundance.txt", sep = "\t", check.names = F, as.is = T)
species <- cbind(data.frame(SampleID = gsub(".*\\.","" , species$sample_ID), ID = gsub("_.*", "", species$sample_ID)),TP =  gsub(".*_","", species$sample_ID),  species[,-ncol(species)])
species$ID <- gsub(".*\\.", "",species$ID)

#genus <- read.delim("../../mb_abundance/W4H_FECAL_genus_clr_transformed_metaphlanAbundance.txt", sep = "\t", check.names = F, as.is = T)
#genus <- cbind(data.frame(SampleID = gsub(".*\\.","" , genus$sample_ID), ID = gsub("_.*", "", genus$sample_ID)),TP =  gsub(".*_","", genus$sample_ID),  genus[,-ncol(genus)])
#genus$ID <- gsub(".*\\.", "",genus$ID)
#genus_sel <- genus[,1:23]

colnames(species) <- gsub("s__", "", colnames(species))
species_names <- c("Parabacteroides_merdae", "Brotolimicola_acetigignens", "Sutterella_wadsworthensis", "Parabacteroides_distasonis", "Bacteroides_faecis", "Roseburia_faecis", "Alistipes_shahii", "Anaerobutyricum_hallii", "Anaerostipes_hadrus", "Phocaeicola_dorei", "Bifidobacterium_adolescentis", "Clostridium_sp_AF36_4", "Phocaeicola_massiliensis", "Bacteroides_caccae", "Bacteroides_uniformis", "Faecalibacterium_SGB15346", "GGB9758_SGB15368", "Ruminococcus_bicirculans", "Dorea_longicatena", "Phocaeicola_vulgatus", "Akkermansia_muciniphila", "Coprococcus_eutactus", "Gemmiger_formicilis", "Blautia_massiliensis", "Alistipes_communis", "Bacteroides_stercoris", "Alistipes_putredinis", "Alistipes_onderdonkii", "Eubacterium_rectale", "Faecalibacterium_prausnitzii", "Oscillibacter_sp_ER4", "Bifidobacterium_longum", "Fusicatenibacter_saccharivorans", "Odoribacter_splanchnicus", "Lachnospira_pectinoschiza", "Barnesiella_intestinihominis", "Candidatus_Cibionibacter_quicibialis", "Clostridiaceae_bacterium", "Blautia_wexlerae", "Bacteroides_ovatus", "Roseburia_hominis", "Blautia_faecis", "Ruminococcus_bromii")

species_sel <- species[,c("SampleID", "ID", "TP", species_names)]


run_lmm_pheno_prot_reverse(species_sel, pheno, 'Bacteroides_uniformis', 'Trigl', subset(covariate_data, select = -c(Batch_hormones)), scale = T)


lmm_res_mb <- data.frame(matrix(nrow = (ncol(d_wide) -4) * (ncol(pheno_old) -4), ncol = 8))
colnames(lmm_res_mb) <- c("prot", "pheno", "estimate", "pval", "se", "tval", "N", "N_unique")
cnt <- 1
for (ph in colnames(pheno_old)[4:ncol(pheno_old)]) {
  cat(ph, "\n")
  for (prot in colnames(species_sel)[4:ncol(species_sel)]){
    res <- run_lmm_pheno_prot_reverse(species_sel, pheno, prot, ph, subset(covariate_data, select = -c(Batch_hormones)), scale = T)
    lmm_res_mb[cnt,] <- c(prot, ph, unlist(res))
    
    cnt <- cnt + 1
  }
}

lmm_res_mb <- na.omit(lmm_res_mb) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
lmm_res_mb$BH_pval <- p.adjust(lmm_res_mb$pval)
