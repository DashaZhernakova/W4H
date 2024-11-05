setwd("/Users/Dasha/work/Sardinia/W4H/olink/data/")
d <- read.delim("pheno_subset.txt", as.is = T, check.names = F, sep = "\t")
d$data_arruolamento <- as.Date(d$data_arruolamento, format = "%d/%m/%Y")
d$data_nascita <- as.Date(d$data_nascita, format = "%d/%m/%Y")
d$age <- as.numeric((d$data_arruolamento - d$data_nascita)/365.25)
d$BMI <- d$peso_kg / (d$altezza_cm / 100)^2

row.names(d) <- sprintf("%03d", d$record_id)

d <- d %>%
  select(-c(record_id, data_arruolamento, data_nascita, altezza_cm, peso_kg ))

write.table(d, file = "pheno_subset.fmt.txt", quote = F, col.names = NA, sep = "\t")
