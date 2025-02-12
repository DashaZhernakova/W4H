library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(ReactomePA)

hs <- org.Hs.eg.db
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
prot_names <- read.delim("data/olink_protein_names.txt", sep = "\t", as.is = T, check.names = F)

#lmm_res_prot_tp <- read.delim("../results/prot_vs_tp_poly3_lmm_adj_age_bmi_preg_storage.txt", sep = "\t", as.is = T, check.names = F)


#sign_pheno <- unique(lmm_res[lmm_res$BH_pval < 0.05, "pheno"])

plot_list <- list()

#input_data <- as.data.frame(pheno_dist)
#input_data$eucl_coef_similarity <- 1 / (1 + input_data$eucl_coef)
#input_data$eucl_similarity <- 1 / (1 + input_data$eucl)
#fold_change <- 'eucl_similarity'
#fold_change <- 'eucl_coef_similarity'

input_data <- gam_res_lin 
#fold_change <- "est_spline"
fold_change <- "est_noTP"


#input_data <- lmm_res
#fold_change <- 'tval'
enrichment_res <- data.frame()
plot_list <- list()
for (ph in unique(input_data$pheno)){
  cat(ph, "\n")
  subs <- input_data[input_data$pheno == ph, c("prot", fold_change)]

  enrich_res <- run_enrichment_analysis(subs)

  gse_res <- data.frame()
  kegg_res <- data.frame()
  rea_res <- data.frame()
  
  if (nrow(enrich_res[['gse']]@result) > 0){
  if ("ONTOLOGY" %in% colnames(enrich_res[['gse']]@result)) { 
    gse_res <-  as_tibble(enrich_res[['gse']]@result)
  } else {
    gse_res <- cbind(data.frame("ONTOLOGY" = "BP"), as_tibble(enrich_res[['gse']]@result))
  }
  } 
  if (nrow(enrich_res[['kegg']]@result) > 0) kegg_res <- cbind(data.frame("ONTOLOGY" = "KEGG"), as_tibble(enrich_res[['kegg']]@result))
  
  if (nrow(enrich_res[['reactome']]@result) > 0) rea_res <- cbind(data.frame("ONTOLOGY" = "Reactome"), as_tibble(enrich_res[['reactome']]@result))
  
  combined_res <- rbind(gse_res, kegg_res, rea_res)
  enrichment_res <- rbind(enrichment_res, cbind(data.frame("Phenotype" = ph), combined_res))
  tryCatch({
    plot_list[[ph]] <- dotplot(enrich_res[['gse']], showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle(ph) + theme(axis.text.y=element_text(size=8))
    plot_list[[paste0(ph, "_kegg")]] <- dotplot(enrich_res[['kegg']], showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle(paste(ph, " KEGG"))
    plot_list[[paste0(ph, "_reactome")]] <- dotplot(enrich_res[['reactome']], showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle(paste(ph, " Reactome"))
  }, error = function(e) {
    message("No results: ", conditionMessage(e))
  })
}

pdf(paste0(out_basedir,"/plots/enrichment_prot_vs_pheno_linear_noTP_gam.pdf"), height = 15, width = 15)
grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)  
dev.off()

#write.table(enrichment_res, file = "../results/prot_vs_pheno_withTP_lmm_adj_covar_enrichment_results.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(enrichment_res, file = paste0(out_basedir,"prot_vs_pheno_linear_no_TP_gam_GO_KEGG_enrichment_results.txt"), quote = F, sep = "\t", row.names = FALSE)


# msigdb
#fold_change <- "eucl_coef_similarity"
#fold_change <- 'eucl_similarity'
enrichment_res <- data.frame()

for (ph in unique(input_data$pheno)){
  cat(ph, "\n")
  #subs <- input_data[input_data$pheno == ph & input_data$signif_tp == T, c("prot", fold_change)]
  subs <- input_data[input_data$pheno == ph, c("prot", fold_change)]
  
  for (cat in c('H', 'C2', 'C5')){
    enrich_res <- run_enrichment_analysis_msigdb(subs, cat)
    #enrich_res <- run_overrepresentation_analysis_msigdb(subs, cat)
    if (nrow(enrich_res) > 0) {
      enrichment_res <- rbind(enrichment_res, cbind(data.frame("Phenotype" = ph, "Category" = cat), enrich_res))
    }
  }
}

write.table(enrichment_res, file = paste0(out_basedir,"prot_vs_pheno_linear_noTP_gam_msigdb_enrichment_results.txt"), quote = F, sep = "\t", row.names = FALSE)


### Functions

add_entrez_uniprot_ids <- function(subs) {
  df_for_enrich <- left_join(subs, prot_names[,c("Assay", "UniProt")], by = c("prot" = "Assay"))
  df_for_enrich <- unique(df_for_enrich)
  # Convert protein names to Entrez gene ids
  
  # remove duplicated proteins from CVD
  df_for_enrich <- df_for_enrich[! grepl("_CVD", df_for_enrich$prot),]
  
  # for proteins that are transcribed from multiple genes take the first gene
  df_for_enrich$gene <- gsub("_.*", "", df_for_enrich$prot)
  
  # keep only BNP and remove NT-proBNP (as the gene is the same)
  df_for_enrich <- df_for_enrich[df_for_enrich$prot != "NT-proBNP", ] 
  
  entrez_ids <- AnnotationDbi::select(hs, 
                       keys = df_for_enrich$gene,
                       columns = c("ENTREZID", "SYMBOL"),
                       keytype = "SYMBOL")
  
  df_for_enrich$entrez <- entrez_ids$ENTREZID
  
  return(df_for_enrich)
}


run_enrichment_analysis <- function(subs){
  
  df_for_enrich <- add_entrez_uniprot_ids(subs)
  # GO
  gene_list <- df_for_enrich[,2]
  names(gene_list) = as.character(df_for_enrich$entrez)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  
  if (all(gene_list > 0)) {
    scoretype = 'pos'
  } else {
    scoretype = 'std'
  }
  gse <- gseGO(geneList=gene_list, 
               ont ="BP", 
               keyType = "ENTREZID", 
               verbose = F, 
               OrgDb = hs, 
               pAdjustMethod = "BH", 
               pvalueCutoff = 0.1,
               scoreType = scoretype)

  gse_res <-  as_tibble(gse@result)
  
  gse_rea <- gsePathway(gene_list, 
                        pvalueCutoff = 0.1,
                        pAdjustMethod = "BH", 
                        verbose = FALSE,
                        scoreType = scoretype)
  
  
  
  # KEGG
  gene_list <- df_for_enrich[,2]
  names(gene_list) = as.character(df_for_enrich$UniProt)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  gse_kegg <- gseKEGG(geneList=gene_list,
                 keyType = "uniprot", 
                 verbose = F, 
                 pAdjustMethod = "BH", 
                 pvalueCutoff = 0.1,
                 scoreType = scoretype)
  
  
  
  return(list("gse" = gse, "kegg" = gse_kegg, "reactome" = gse_rea))
}

run_enrichment_analysis_msigdb <- function(subs, cat = 'H'){
  df_for_enrich <- add_entrez_uniprot_ids(subs)
  
  
  gene_list <- df_for_enrich[,2]
  names(gene_list) = as.character(df_for_enrich$entrez)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  msigdb_subset <- msigdbr(species = "Homo sapiens", category = cat)
  msigdb_subset_genes <- msigdb_subset %>% 
    dplyr::select(gs_name, entrez_gene)
  msigdb_subset_descr <- msigdb_subset %>% 
    dplyr::select(gs_name, gs_subcat, gs_description)
  msigdb_subset_descr <- unique(msigdb_subset_descr)
    
  em <- GSEA(gene_list, 
              TERM2GENE = msigdb_subset_genes, 
              scoreType = "pos", 
              pvalueCutoff = 0.1,
              pAdjustMethod = "BH",
              verbose = F)
  em <- dplyr::left_join(as_tibble(em@result), msigdb_subset_descr, by = c("ID" = "gs_name"))

  return(em)
}


###############################################################################
# Overrepresentation analysis
###############################################################################

gam_res_prot_tp <- read.delim(paste0(out_basedir, "prot_vs_tp_gam_adj_age_bmi_preg_storage.txt"), as.is = T, sep = "\t")

tmp <- add_entrez_uniprot_ids(gam_res_prot_tp)

for (cl in unique(tmp$trajectory_cluster)){
  signif_genes <- na.omit(tmp[tmp$trajectory_cluster == cl, "entrez"])
  background_genes <- na.omit(tmp[tmp$gam_BH_pval < 0.05, "entrez"])
  if(length(signif_genes) < 10) next
  res <- run_overrepresentation_analysis(signif_genes, background_genes)
  if (nrow(res) > 0) print(cl)
}

genename_to_entrez <- function(gene_names){
  entrez_ids <- select(hs, 
                       keys = gene_names,
                       columns = c("ENTREZID", "SYMBOL"),
                       keytype = "SYMBOL")
  return(entrez_ids$ENTREZID)
}

hs_df_kegg <- hs_msigdb_df %>%
  dplyr::filter(
    gs_subcat == "CP:KEGG" 
  )
hs_df_reactome <- hs_msigdb_df %>%
  dplyr::filter(
    gs_subcat == "CP:REACTOME"
  )

run_overrepresentation_analysis <- function(signif_genes, background_genes, pval_threshold = 0.05){

  #df_for_enrich <- add_entrez_uniprot_ids(subs)
  
  #signif_genes <- df_for_enrich[df_for_enrich[,col_name] < threshold , 'entrez']
  #background_genes <- df_for_enrich$entrez

  if (length(signif_genes) < 10) return(NULL)

  enrich_res_kegg <- enricher(
    gene = signif_genes, # A vector of your genes of interest
    pvalueCutoff = 0.1, # Can choose a FDR cutoff
    pAdjustMethod = "BH", # Method to be used for multiple testing correction
    universe = background_genes,
    TERM2GENE = dplyr::select(
      hs_df_kegg,
      gs_name,
      human_entrez_gene
    )) 
  
  enrich_res_reactome <- enricher(
    gene = signif_genes, # A vector of your genes of interest
    pvalueCutoff = 0.1, # Can choose a FDR cutoff
    pAdjustMethod = "BH", # Method to be used for multiple testing correction
    universe = background_genes,
    TERM2GENE = dplyr::select(
      hs_df_reactome,
      gs_name,
      human_entrez_gene
    )) 
  
  result_df_kegg <- data.frame(enrich_res_kegg@result)
  result_df_reactome <- data.frame(enrich_res_reactome@result)
  
  result_df_kegg <- result_df_kegg[result_df_kegg$pvalue < 0.05,]
  result_df_reactome <- result_df_reactome[result_df_reactome$pvalue < 0.05,]
  
  combined_res <- rbind(result_df_reactome, result_df_kegg)
  return(combined_res)
}



run_overrepresentation_analysis_msigdb <- function(subs, cat){
  
  df_for_enrich <- add_entrez_uniprot_ids(subs)

  signif_genes <- df_for_enrich[df_for_enrich$prot %in% all_prots & df_for_enrich[,2] > quantile(df_for_enrich[,2], 0.95), 'entrez']
  background_genes <- df_for_enrich[df_for_enrich$prot %in% all_prots & df_for_enrich[,2] < quantile(df_for_enrich[,2], 0.95), 'entrez']
  
  if (length(signif_genes) < 10) return(NULL)
  
  msigdb_subset <- msigdbr(species = "Homo sapiens", category = cat)
  msigdb_subset_genes <- msigdb_subset %>% 
    dplyr::select(gs_name, entrez_gene)
  msigdb_subset_descr <- msigdb_subset %>% 
    dplyr::select(gs_name, gs_subcat, gs_description)
  msigdb_subset_descr <- unique(msigdb_subset_descr)
  
  
  enrich_res <- enricher(
    gene = signif_genes, # A vector of your genes of interest
    pvalueCutoff = pval_threshold, # Can choose a FDR cutoff
    pAdjustMethod = "BH", # Method to be used for multiple testing correction
    universe = background_genes,
    TERM2GENE = msigdb_subset_genes
      ) 
  
  enrich_res <- enrich_res[enrich_res$pvalue < 0.05,]
  return(enrich_res)
  
}
