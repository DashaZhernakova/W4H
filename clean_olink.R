library(dplyr)
library(ggplot2)

library(tidyr)
library(pheatmap)
library(patchwork)
library(tibble)
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")

setwd("/Users/Dasha/work/Sardinia/W4H/olink/data/")


################################################################################
# 1. Remove technical samples and assays from the data table
################################################################################


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

################################################################################
# 2. Remove 100_2 - the sample that failed CVD completely and partially INF
################################################################################

inf2 <- inf2[inf2$SampleID != '100_2',]
cvd2 <- cvd2[cvd2$SampleID != '100_2',]


################################################################################
# 3. Combine INF and CVD panels
################################################################################

merged_long <- rbind(inf2[,c("SampleID", "Assay", "NPX", "Panel")], cvd2[,c("SampleID", "Assay", "NPX", "Panel")])
merged_long$Assay2 <- paste0(merged_long$Assay, "_", merged_long$Panel)
merged_wide <-merged_long[,c("SampleID", "Assay2", "NPX")] %>%
  pivot_wider(names_from = Assay2, values_from = NPX)


# Combine clean without panel name
# For the overlapping proteins keep the INF version as default and add a suffix to the CVD version:
overlap <- intersect(unique(inf2$Assay), unique(cvd$Assay))
merged_long2 <- merged_long
merged_long2[merged_long2$Assay2 %in% c(paste0(overlap, "_Cardiometabolic")), "Assay"] <- merged_long2[merged_long2$Assay2 %in% c(paste0(overlap, "_Cardiometabolic")), "Assay2"]
merged_long2$Assay2 <- NULL
merged_long2$Panel <- NULL
merged_long2$Timepoint <- gsub(".*_","",merged_long2$SampleID)
merged_long2$ID <- gsub("_.*","",merged_long2$SampleID)
merged_long2$Assay <- gsub("_Cardiometabolic", "_CVD", merged_long2$Assay)
merged_wide2 <-merged_long2[,c("SampleID", "Assay", "NPX")] %>%
  pivot_wider(names_from = Assay, values_from = NPX)

num_cols_with_na <- sum(colSums(is.na(merged_wide2)) > 0)
num_rows_with_na <- sum(rowSums(is.na(merged_wide2)) > 0)

num_rows_with_na
num_cols_with_na

################################################################################
# 4. Write combined cleaned tables
################################################################################

ID <- gsub("_.*", "", merged_wide2$SampleID)
TP <- gsub(".*_", "", merged_wide2$SampleID)

merged_wide2 <- cbind(ID, TP, merged_wide2)

write.table(merged_wide2, file = "olink_clean_CVD+INF.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(merged_long2, file = "olink_clean_CVD+INF_long.txt", quote = F, sep = "\t", row.names = FALSE)

################################################################################
# 5. Set extreme outliers > 4 SDs from the mean to NA
################################################################################

remove_outliers_per_feature <- function(d, sd_cutoff = 4, iqr_cutoff = 3, method = 'zscore') {
  if (method == 'zscore'){
    zscore <- scale(d)  # Compute z-scores
    d[abs(zscore) > sd_cutoff] <- NA  # Replace outliers with NA
  } else if (method == 'IQR'){
    q <- quantile(d, probs = c(0.25, 0.75), na.rm = TRUE)
    iqr <- diff(q)
    d[d < (q[1] - iqr_cutoff * iqr) | d > (q[2] + iqr_cutoff * iqr)] <- NA
  } else {
    stop("Wrong method, should be zscore or IQR")
  }
  return(d)
}

remove_outliers_dataframe <- function(df, sd_cutoff = 5) {
  # Create a copy of the data frame to store the outlier mask
  outlier_summary <- df %>%
     select(-c(SampleID, ID, TP)) %>%  # Exclude specific columns
     summarise(across(everything(), ~any(abs(scale(.)) > sd_cutoff, na.rm = TRUE))) %>%
     pivot_longer(cols = everything(), names_to = "column", values_to = "has_outliers")
  
  # Apply the outlier removal
  df_cleaned <- df %>%
    mutate(across(.cols = -c(SampleID, ID, TP), .fns = ~remove_outliers_per_feature(., sd_cutoff)))
  
  # Return the cleaned data and the outlier mask
  list(cleaned_data = df_cleaned, outlier_mask = outlier_summary)
}

res <- remove_outliers_dataframe(merged_wide2,  sd_cutoff = 5)

cleaned_data <- res$cleaned_data
outliers <- res$outlier_mask[res$outlier_mask$has_outliers == T, ]$column
length(outliers)

write.table(cleaned_data, file = "olink_clean_CVD+INF_rm_outliers_5sd.txt", quote = F, sep = "\t", row.names = FALSE)


plots <- plot_features_with_outliers(merged_wide2, outliers, cutoff = 5)


plot_features_with_outliers <- function(df, features_with_outliers, cutoff = 5, output_dir = NULL, method = 'zscore') {
  
  # Iterate over each feature and create a plot
  plots <- list()
  for (feature in features_with_outliers) {
    feature_values <- df[,feature]
    if (method == 'zscore'){
      zscore <- scale(feature_values)
      threshold_upper <- mean(feature_values, na.rm = TRUE) + cutoff * sd(feature_values, na.rm = TRUE)
      threshold_lower <- mean(feature_values, na.rm = TRUE) - cutoff * sd(feature_values, na.rm = TRUE)
    } else {
      q <- quantile(feature_values, probs = c(0.25, 0.75), na.rm = TRUE)
      iqr <- diff(q)
      
      threshold_lower <- q[1] - cutoff * iqr
      threshold_upper <- q[2] + cutoff * iqr
    }
    # Create the plot
    p <- ggplot(df, aes(x = SampleID, y = .data[[feature]], color = ID)) +
      geom_point() +
      geom_hline(yintercept = c(threshold_lower, threshold_upper), color = "red", linetype = "dashed") +
      labs(
        title = paste("Outliers in", feature),
        x = "SampleID",
        y = feature
      ) +
      theme_minimal() +
      theme(legend.position="none")
    
    # Save the plot if output_dir is specified
    if (!is.null(output_dir)) {
      ggsave(filename = file.path(output_dir, paste0(feature, "_outliers.png")), plot = p)
    }
    
    # Store the plot in a list
    plots[[feature]] <- p
  }
  
  return(plots)
}
