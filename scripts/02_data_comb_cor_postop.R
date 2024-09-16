###########################################################################
# PhD
# 02_data_comb_cor_postop.R
###########################################################################

# Load the dataset
source("scripts/02_data_comb_postop.R")

# Pre-filter Correlation between mRNA variables
cor1_df <- cor(comb_data[11:8826], comb_data[11:8826])
fc1 <- findCorrelation(cor1_df, cutoff = 0.6, names = TRUE)
cor2_df <- cor(comb_data[11:8826], comb_data[8827:17641])
fc2 <- findCorrelation(cor2_df, cutoff = 0.6, names = TRUE)
cor3_df <- cor(comb_data[11:8826], comb_data[17642:26457])
fc3 <- findCorrelation(cor3_df, cutoff = 0.6, names = TRUE)
cor4_df <- cor(comb_data[8827:17641], comb_data[8827:17641])
fc4 <- findCorrelation(cor4_df, cutoff = 0.6, names = TRUE)
cor5_df <- cor(comb_data[8827:17641], comb_data[17642:26457])
fc5 <- findCorrelation(cor5_df, cutoff = 0.6, names = TRUE)
cor6_df <- cor(comb_data[17642:26457], comb_data[17642:26457])
fc6 <- findCorrelation(cor6_df, cutoff = 0.6, names = TRUE)

fc <- c(fc1, fc2, fc3, fc4, fc5, fc6)
fc <- fc[!duplicated(fc)]

cor_df <- comb_data %>%
  select(!(fc))

# Remove intermediary datasets
rm(comb_data, cor1_df, cor2_df, cor3_df, cor4_df, cor5_df,
   cor6_df, fc1, fc2, fc3, fc4, fc5, fc6, fc)


