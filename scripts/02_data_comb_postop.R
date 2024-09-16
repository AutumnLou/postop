###########################################################################
# PhD
# 02_data_comb_postop.R
###########################################################################

# Load the dataset
source("scripts/02_data_msc.R")
source("scripts/02_data_rna.R")

#Combine the Clinical and mRNA data
comb_data <- merge(msc_data, rna_data, by = 0)
rownames(comb_data ) <- comb_data [, 1]
comb_data <- comb_data [, -1]
comb_data[,11:ncol(comb_data)] <- lapply(comb_data[,11:ncol(comb_data)], function(x) {
  if(is.character(x)) as.numeric(x)
})

# Remove intermediary datasets
rm(rna_data, msc_data)
