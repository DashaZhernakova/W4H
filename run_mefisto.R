library(reticulate)
use_python("/Users/Dasha/miniconda3/bin/python", required = TRUE)
library(MOFA2)



d_long2 <- pivot_longer(d_wide, cols = 4:ncol(d_wide), names_to = 'feature', values_to = 'value')
#d_long2$SampleID <- NULL
d_long2$view <- 'protein'
colnames(d_long2) <- gsub("SampleID", "sample", colnames(d_long2))

pheno_long <- pivot_longer(pheno, cols = 4:ncol(pheno), names_to = 'feature', values_to = 'value')
#pheno_long$SampleID <- NULL
pheno_long$view <- 'lipid'
hormones <- c("PROG","LH","FSH","X17BES","PRL")
pheno_long[pheno_long$feature %in% hormones, "view"] <- 'hormone'
colnames(pheno_long) <- gsub("SampleID", "sample", colnames(pheno_long))

all_combined_long <- rbind(d_long2, pheno_long)
all_combined_long$TP <- as.numeric(all_combined_long$TP)
# Read data
MOFAobject_untrained <- create_mofa(data = all_combined_long)

MOFAobject_untrained <- set_covariates(MOFAobject_untrained, covariates = "TP")

plot_data_overview(MOFAobject_untrained,
                   show_covariate = TRUE,
                   show_dimensions = TRUE) 

# Prepare the MOFA object

data_opts <- get_default_data_options(MOFAobject_untrained)
data_opts$scale_views <- T

model_opts <- get_default_model_options(MOFAobject_untrained)
model_opts$num_factors <- 5

train_opts <- get_default_training_options(MOFAobject_untrained)
train_opts$maxiter <- 100

mefisto_opts <- get_default_mefisto_options(MOFAobject_untrained)

MOFAobject_untrained <- prepare_mofa(MOFAobject_untrained, model_options = model_opts,
                   mefisto_options = mefisto_opts,
                   training_options = train_opts,
                   data_options = data_opts)

# Run MOFA
outfile <- "../results/mefisto/test_model.hdf5"
MOFAobject_untrained <- run_mofa(MOFAobject_untrained, outfile = outfile, use_basilisk = TRUE)
