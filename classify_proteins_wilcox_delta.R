library(dplyr)

classify_median_wilcox <- function(d_wide, all_prots, wilcox_p_threshold = 0.05){
  subs <- d_wide[,c("ID", "TP", all_prots)]
  subs$TP <- paste0("TP", subs$TP)
  
  subs2 <- subs %>% 
    pivot_longer(cols = 3:ncol(subs), names_to = "prot") %>%
    pivot_wider(names_from = TP, values_from = value, names_sort = TRUE)
  
  classify_change <- function(x1, x2) {
    test <- wilcox.test(x1, x2, paired = TRUE)
    if (test$p.value < wilcox_p_threshold) {
      direction <- ifelse(median(x2, na.rm = T) > median(x1, na.rm = T), "/", "\\")
    } else {
      direction <- "-"
    }
    return(direction)
  }

  classified <- subs2 %>%
    group_by(prot) %>%
    summarize(
      Change1_2 = classify_change(TP1, TP2),
      Change2_3 = classify_change(TP2, TP3),
      Change3_4 = classify_change(TP3, TP4),
      Change1_4 = classify_change(TP1, TP4),
      .groups = "drop"
    ) %>%
    mutate(
      Pattern = if_else(
        Change1_2 == "-" & Change2_3 == "-" & Change3_4 == "-" ,
        paste(Change1_4, sep = " "),
        paste(Change1_2, Change2_3, Change3_4, sep = " ")
      )
    )
  
  classified <- as.data.frame(classified)
  #cl <- classified$Pattern
  #names(cl) <- classified$prot
  # plot_clusters(cl, "", 0,save_pdf = F, add_cluster_name = T)  
  
  return(classified)
}

