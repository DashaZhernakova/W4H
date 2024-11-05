library(factoextra)

sim_matrix <- prot_sim_wide
dist_matrix <- prot_dist_wide

plot_clusters <- function(cl, method, num_k, colored = F, signif = NULL){
  plot_list = list()
  for (cluster in unique(cl)){
    plot_list[[cluster]] <- plot_traj_many_prots2(prot_trajs, names(cl[cl == cluster]), colored = colored, signif)
  }
  pdf(paste0("../plots/clustering_signif/", method, "_k", num_k, ".pdf"), width = 20, height = 20)
  if (length(unique(cl)) < 10){
    grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)  
  } else {
    grid.arrange(grobs = plot_list, ncol = 4, nrow = 4)  
  }
  dev.off()
  
}

prot_deriv <- as.data.frame(t(apply(prot_coefs_raw, 1, get_derivatives)))

### get pairwise deltas
prot_deltas <- data.frame(matrix(nrow = nrow(prot_trajs), ncol = 3))
colnames(prot_deltas) <- c("d12", "d23", "d34")
row.names(prot_deltas) <- row.names(prot_trajs)
prot_deltas$d12 <- prot_trajs[,2] - prot_trajs[,1]
prot_deltas$d23 <- prot_trajs[,3] - prot_trajs[,2]
prot_deltas$d34 <- prot_trajs[,4] - prot_trajs[,3]


signif = lmm_res_prot_tp[lmm_res_prot_tp$BH_pval < 0.05, "prot"]
signif_dist <- intersect(signif, row.names(dist_matrix))

for (num_k in 3:5){
  #
  # k-means
  #
  
  # trajectories
  method = 'kmeans_traj'
  km <- kmeans(prot_trajs[signif,], num_k, nstart = 50)
  cl <- km$cluster
  plot_clusters(cl, method, num_k)
  
  # coefficients
  method = 'kmeans_coef'
  km <-kmeans(prot_coefs[signif,2:4], num_k, nstart = 50)
  cl <- km$cluster
  plot_clusters(cl, method, num_k)
  
  #method = 'kmeans_coef_raw'
  km <-kmeans(prot_coefs_raw[signif,2:4], num_k, nstart = 50)
  cl <- km$cluster
  plot_clusters(cl, method, num_k)
  
  # 1st derivatives
  method = 'kmeans_deriv1'
  km <-kmeans(prot_deriv[signif,1:10], num_k, nstart = 50)
  cl <- km$cluster
  plot_clusters(cl, method, num_k)
  
  # 2nd derivatives
  method = 'kmeans_deriv2'
  km <-kmeans(prot_deriv[signif,10:20], num_k, nstart = 50)
  cl <- km$cluster
  plot_clusters(cl, method, num_k)
  

  
  #
  # PAM
  #
  #method = 'pam_my_eucl'
  pam_res <- cluster::pam(dist_matrix[signif_dist,signif_dist], num_k, diss = T)
  cl <- pam_res$clustering
  plot_clusters(cl, method, num_k)
  
  method = 'pam_coef2'
  pam_res <- cluster::pam(prot_coefs[signif,2:4], num_k, diss = F)
  cl <- pam_res$clustering
  plot_clusters(cl, method, num_k)  
  
  #method = 'pam_coef_raw'
  #pam_res <- cluster::pam(prot_coefs_raw[,2:4], num_k, diss = F)
  #cl <- pam_res$clustering
  #plot_clusters(cl, method, num_k)  
      
  method = 'pam_traj2'
  pam_res <- cluster::pam(prot_trajs[signif,], num_k, diss = F)
  cl <- pam_res$clustering
  plot_clusters(cl, method, num_k)    
  
  
  #
  # hclust
  #
  
  #method = 'hclust_my_eucl'
  hc <- hclust(as.dist(dist_matrix[signif_dist,signif_dist]), method = 'complete')
  cl <- cutree(hc, k = num_k)
  plot_clusters(cl, method, num_k)    
  
  method = 'hclust_coef'
  hc <- hclust(dist(prot_coefs[signif,2:4], method = 'euclidean'), method = 'complete')
  cl <- cutree(hc, k = num_k)
  plot_clusters(cl, method, num_k)    
  
  method = 'hclust_traj'
  hc <- hclust(dist(prot_trajs[signif,], method = 'euclidean'), method = 'complete')
  cl <- cutree(hc, k = num_k)
  plot_clusters(cl, method, num_k)    
  
  
  #
  # KML
  #
  method = 'kml_traj'
  cld <- clusterLongData(traj = as.matrix(prot_trajs[signif,]), time = as.numeric(colnames(prot_trajs)))
  kml(cld, nbRedrawing = 20, nbClusters = num_k)
  cl <- getClusters(cld, num_k, asInteger = T)
  names(cl) <- row.names(prot_trajs[signif,])
  plot_clusters(cl, method, num_k)    

  
  # Deltas and simplified coefficients
  
  method = 'pam_delta'
  pam_res <- cluster::pam(prot_deltas[signif,], num_k, diss = F)
  cl <- pam_res$clustering
  plot_clusters(cl, method, num_k)  
  
  method = 'pam_simpl'
  pam_res <- cluster::pam(prot_simp[signif,], num_k, diss = F)
  cl <- pam_res$clustering
  plot_clusters(cl, method, num_k,)  
  
  method = 'pam_coef_of'
  pam_res <- cluster::pam(prot_coef_of[signif,], num_k, diss = F)
  cl <- pam_res$clustering
  plot_clusters(cl, method, num_k)  
  
  method = 'hclust_delta'
  hc <- hclust(dist(prot_deltas[signif,], method = 'euclidean'), method = 'complete')
  cl <- cutree(hc, k = num_k)
  plot_clusters(cl, method, num_k)  
  
  method = 'hclust_simpl'
  hc <- hclust(dist(prot_simp[signif,], method = 'euclidean'), method = 'complete')
  cl <- cutree(hc, k = num_k)
  plot_clusters(cl, method, num_k)  
  
  
}




pheatmap(t(prot_coefs_all), cutree_cols = num_k)
pheatmap(sim_matrix, clustering_distance_rows = as.dist(dist_matrix), clustering_distance_cols = as.dist(dist_matrix), cutree_cols = num_k)



### For each prot get x of inflation point, the roots, and x where the 1st derivative is 0

