d <- read.delim("ukb_pqtls.txt", sep = "\t", as.is = T, check.names = F)

df_summary <- d %>%
  group_by(Assay, cis_trans) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = cis_trans, values_from = count, values_fill = 0)

write.table(df_summary, file = "ukb_pqtls_counts.txt",  quote = F, sep = "\t", row.names = FALSE)
