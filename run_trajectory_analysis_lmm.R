
run_trajectory_analysis_lmm <- function(d_wide, all_prots, out_basedir, n_points = 100){

  n_prots <- length(all_prots)
  prot_coefs <- data.frame(matrix(nrow =  n_prots, ncol = 3))
  row.names(prot_coefs) <- all_prots
  
  prot_trajs <- data.frame(matrix(nrow = n_prots , ncol = n_points))
  row.names(prot_trajs) <- all_prots
  colnames(prot_trajs) <- seq(1,4, length.out = n_points)
  
  prot_dist <- data.frame(matrix(nrow = n_prots * n_prots, ncol = 4))
  colnames(prot_dist) <- c("prot1", "prot2", "eucl",  "eucl_coef")
  
  prot1_passed <- ""
  
  pb = txtProgressBar(min = 0, max = n_prots, initial = 0) 
  stepi = 0
  cnt <- 1
  for (prot1 in all_prots){
    #print(prot1)
    lmm_fit1 <- fit_lmm_poly3_adj_covar(d_wide, prot1, n = n_points, covariates, scale = T)
    prot_coefs[prot1,] <- lmm_fit1$coefficients
    prot_trajs[prot1,] <- lmm_fit1$predicted
    
    prot1_passed <- c(prot1_passed, prot1)
    
    setTxtProgressBar(pb,stepi)
    stepi <- stepi + 1
    for (prot2 in all_prots){
      if (prot2 %in% prot1_passed) next
      
      lmm_fit2 <- fit_lmm_poly3_adj_covar(d_wide, prot2, n = n_points, covariates, scale = T)
      
      # Euclidean distance between trajectories
      eucl_dist <- sum(abs(as.numeric(lmm_fit1$predicted) - as.numeric(lmm_fit2$predicted)))/length(lmm_fit1$predicted)
      
      #Euclidean distance in the lm coefficients space
      eucl_dist_coef <- TSdist::EuclideanDistance(lmm_fit1$coefficients, lmm_fit2$coefficients)
      
      prot_dist[cnt,] <- c(prot1, prot2, eucl_dist, eucl_dist_coef)
      cnt <- cnt + 1
    }
  }
  cat("\n")
  close(pb)
  
  prot_dist <- na.omit(prot_dist) %>%
    mutate(across(-c( prot1, prot2), as.numeric)) 
  dist_matrix <- my_pivot_wider(prot_dist, "prot1", "prot2", "eucl")
  dist_matrix[lower.tri(dist_matrix)] <- t(dist_matrix)[lower.tri(dist_matrix)]
  
  dist_matrix_coef <- my_pivot_wider(prot_dist, "prot1", "prot2", "eucl_coef")
  dist_matrix_coef[lower.tri(dist_matrix_coef)] <- t(dist_matrix_coef)[lower.tri(dist_matrix_coef)]
  
  #convert distance to similarity metric
  prot_dist$eucl_similarity <- 1 / (1 + prot_dist$eucl)
  prot_dist$eucl_similarity_coef <- 1 / (1 + prot_dist$eucl_coef)
  
  sim_matrix <- my_pivot_wider(prot_dist, "prot1", "prot2", "eucl_similarity")
  sim_matrix[lower.tri(sim_matrix)] <- t(sim_matrix)[lower.tri(sim_matrix)]
  
  sim_matrix_coef <- my_pivot_wider(prot_dist, "prot1", "prot2", "eucl_similarity_coef")
  sim_matrix_coef[lower.tri(sim_matrix_coef)] <- t(sim_matrix_coef)[lower.tri(sim_matrix_coef)]
  
  dir.create(file.path(out_basedir, "trajectories_lmm"), showWarnings = FALSE)
  write.table(sim_matrix, file = paste0(out_basedir, "trajectories_lmm/similarity_matrix_eucl_", n_prots, "_prots.txt"), quote = F, sep = "\t", row.names = FALSE)
  write.table(dist_matrix, file = paste0(out_basedir, "trajectories_lmm/distance_matrix_eucl_", n_prots, "_prots.txt"), quote = F, sep = "\t", row.names = FALSE)
  write.table(prot_dist, file = paste0(out_basedir, "trajectories_lmm/protein_distance_eucl_", n_prots, "_prots.txt"), quote = F, sep = "\t", row.names = FALSE)
  write.table(prot_trajs, file = paste0(out_basedir, "trajectories_lmm/protein_trajectories_", n_prots, "_prots.txt"), quote = F, sep = "\t", col.names = NA, row.names = T)
  
  
  # Heatmap and clustering on euclidean distance between curves
  method = 'hclust_my_eucl_dist'
  num_k = 8
  hm <- pheatmap(dist_matrix, cutree_rows = num_k, cutree_cols = num_k, filename = paste0(out_basedir,"trajectories_lmm/", method, "_k", num_k, "_heatmap.pdf"))
  cl  = cutree(hm$tree_row, k = num_k)
  plot_clusters(cl, method, prot_trajs = prot_trajs, num_k, out_path = paste0(out_basedir,"trajectories_lmm/", method, "_k", num_k, "_", n_prots, "_prots.pdf"))
  write.table(as.data.frame(cl), file = paste0(out_basedir, "trajectories_lmm/", method, "_k", num_k, "_", n_prots, "_prots.txt"), quote = F, sep = "\t")
  
  # Get pairwise differences between i and i+1 timepoint and cluster according to them
  prot_deltas <- data.frame(matrix(nrow = nrow(prot_trajs), ncol = ncol(prot_trajs) - 1))
  
  row.names(prot_deltas) <- row.names(prot_trajs)
  
  for (c in 1:(ncol(prot_trajs) - 1)){
    prot_deltas[,c] <- prot_trajs[,(c + 1)] - prot_trajs[,c]
  }
  method = 'pam_delta'
  num_k = 6
  #num_k = 15
  pam_res <- cluster::pam(prot_deltas, num_k, diss = F)
  cl <- pam_res$clustering
  plot_clusters(cl, method, num_k, out_path = paste0(out_basedir,"trajectories_lmm/", method, "_k", num_k, "_", n_prots, "_prots.pdf"))
  write.table(as.data.frame(cl), prot_trajs = prot_trajs, file = paste0(out_basedir, "trajectories_lmm/", method, "_k", num_k, "_", n_prots, "_prots.txt"), quote = F, sep = "\t")
  
}
