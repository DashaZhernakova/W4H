library(dplyr)

setwd("/Users/Dasha/work/Sardinia/W4H/olink/")

########
## Format diet data
########

#diet <- read.delim("../phenotypes/DietFreeze_DataCleaned_Feb2025/FFQ_DIARI_ms.csv", as.is = T, sep = ",")
diet <- read.delim("../phenotypes/DietFreeze_DataCleaned_Feb2025/FFQ_DIARI_3gg.csv", as.is = T, sep = ",")

diet$ID <- sprintf("%03d", diet$ID)

ffq <- diet[,c(1, grep("FFQ", colnames(diet)))]

write.table(ffq, file = "../phenotypes/DietFreeze_DataCleaned_Feb2025/fmt_dasha/FFQ.txt", quote = F, sep = "\t", row.names = FALSE)

diet_w <- diet[,-grep("FFQ", colnames(diet))]
diet_l <- diet_w %>%
  pivot_longer(
    cols = -ID,
    names_to = c(".value", "TP"),
    names_sep = "g_"
  ) %>%
  mutate(TP = as.numeric(factor(TP, levels = c("first", "second", "third", "fourth"))))

#write.table(diet_l, file = "../phenotypes/DietFreeze_DataCleaned_Feb2025/fmt_dasha/FFQ_DIARI_ms_long.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(diet_l, file = "../phenotypes/DietFreeze_DataCleaned_Feb2025/fmt_dasha/FFQ_DIARI_3gg_long.txt", quote = F, sep = "\t", row.names = FALSE)



#
# OLD
#

diet <- read.delim("FFQ_visit0.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(record_id = "character"))

colnames(diet) <- gsub("pane_consumo_sett", "pane_freq_sett", colnames(diet))
colnames(diet) <- gsub("cerdure_freq_sett", "verdure_freq_sett", colnames(diet))
colnames(diet) <- gsub("caffe\\b", "caffe_consumo", colnames(diet))
colnames(diet) <- gsub("alcol\\b", "alcol_consumo", colnames(diet))


# INTAKE

intake <- cbind(diet$record_id, diet[,grepl("consumo", colnames(diet))])
colnames(intake)[1] <- 'ID'
inake_counts <- intake %>%
  pivot_longer(-ID) %>%
  select(-ID) %>%
  pivot_wider(names_from = value, values_fn = length, names_sort=TRUE)

# remove 3 types of milk consumption because they have different values from the rest
diet$latte_consumo_tipo___1 <- NULL
diet$latte_consumo_tipo___2 <- NULL
diet$latte_consumo_tipo___3 <- NULL

# remove products with less than 10 samples per factor level 
diet$formaggi_consumo <- NULL
diet$uova_consumo <- NULL
diet$pasta_riso_consumo <- NULL
diet$verdure_consumo <- NULL
diet$frutta_consumo <- NULL
diet$olio_oliva_consumo <- NULL
diet$margarina_consumo <- NULL

intake <- cbind(diet$record_id, diet[,grepl("consumo", colnames(diet))])
colnames(intake)[1] <- 'ID'
inake_counts <- intake %>%
  pivot_longer(-ID) %>%
  select(-ID) %>%
  pivot_wider(names_from = value, values_fn = length, names_sort=TRUE)

write.table(intake, "FFQ_visit0_intake.txt", quote = F, sep = "\t", row.names = FALSE)


#
# FREQUENCY WEEKLY
#
freq <- cbind(diet$record_id, diet$caffe_quantita_die, diet[,grepl("freq_sett", colnames(diet))])
colnames(freq)[c(1, 2)] <- c('ID', 'caffe_quantita_die')
freq_counts <- freq %>%
  pivot_longer(-ID) %>%
  select(-ID) %>%
  pivot_wider(names_from = value, values_fn = length, names_sort=TRUE)

diet$margarina_freq_sett <- NULL
diet[diet > 14] <- 14

freq <- cbind(diet$record_id, diet$caffe_quantita_die, diet[,grepl("freq_sett", colnames(diet))])
colnames(freq)[c(1, 2)] <- c('ID', 'caffe_quantita_die')


write.table(freq, "FFQ_visit0_freq_sett.txt", quote = F, sep = "\t", row.names = FALSE)


#
# FREQUENCY OVERALL
#
freq <- cbind(diet$record_id,  diet[,grepl("frequenza", colnames(diet))])
colnames(freq)[1] <- c('ID')
freq_counts <- freq %>%
  pivot_longer(-ID) %>%
  select(-ID) %>%
  pivot_wider(names_from = value, values_fn = length, names_sort=TRUE)

diet$verdure_frequenza <- NULL
diet$latte_frequenza <- NULL
diet$formaggi_frequenza <- NULL
diet$frutta_frequenza <- NULL
diet$olio_oliva_frequenza <- NULL
diet$margarina_frequenza <- NULL

freq <- cbind(diet$record_id,  diet[,grepl("frequenza", colnames(diet))])
colnames(freq)[1] <- c('ID')

write.table(freq, "FFQ_visit0_frequenza.txt", quote = F, sep = "\t", row.names = FALSE)

colnames(diet)[1] <- c('ID')
diet_counts <- diet %>%
  pivot_longer(-ID) %>%
  select(-ID) %>%
  pivot_wider(names_from = value, values_fn = length, names_sort=TRUE)
write.table(diet, "FFQ_visit0_clean.txt", quote = F, sep = "\t", row.names = FALSE)
