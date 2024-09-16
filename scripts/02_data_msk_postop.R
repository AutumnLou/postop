###########################################################################
# MSKCC
# 02_data_msk_postop.R
###########################################################################

# Read in the required packages ------------------------------------------

source("scripts/01_data_clean.R")

# Create a dataset for use with the MSKCC Nomogram
msk_nomo_data <- final_data

# Remove all unused columns
msk_nomo_data <- msk_nomo_data %>%
  
  # Remove Race, PreTxPSA, PathStage
  select(!c(Race, PreTxPSA, PathStage))

# Create new variables
msk_nomo_data <- msk_nomo_data %>%
  mutate(biop_grp = ifelse(BxGGS == 6,
                           1,
                           ifelse(BxGG1 == 3 & BxGG2 == 4,
                                  2,
                                  ifelse(BxGG1 == 4 & BxGG2 == 3,
                                         3,
                                         ifelse(BxGG1 == 4 & BxGG2 == 4,
                                                4, 5)))),
         path_grp = ifelse(PathGGS == 6,
                           1,
                           ifelse(PathGG1 == 3 & PathGG2 == 4,
                                  2,
                                  ifelse(PathGG1 == 4 & PathGG2 == 3,
                                         3,
                                         ifelse(PathGG1 == 4 & PathGG2 == 4,
                                                4, 5)))),
         stage = mapvalues(ClinT_Stage,
                           from = c("T1C", "T2", "T2A",
                                    "T2B", "T2C", "T3",
                                    "T3A", "T3B", "T3C"),
                           to = c("1C", "2A", "2A",
                                  "2B", "2C", "3p",
                                  "3p", "3p", "3p")),
         ece = ifelse(ECE == 'Absent', 0, 1),
         svi = ifelse(SVI == 'Negative', 0, 1),
         sms = ifelse(SMS == 'Negative', 0, 1),
         lni = ifelse(LNI == 'Abnormal_N1', 1,
                      ifelse(LNI == 'Normal_N0', 0, 0.5))
         
  )

# Remove the original columns
msk_nomo_data <- msk_nomo_data %>%
  select(!c(BxGG1, BxGG2, BxGGS,
            PathGG1, PathGG2, PathGGS,
            ClinT_Stage,
            ECE, SVI, SMS, LNI))

# Widen the dataset
msk_wide_data <- msk_nomo_data %>%
  dummy_cols(select_columns = c("biop_grp", "path_grp","stage")) %>%
  
  # Remove original columns
  select(!c(biop_grp, path_grp, stage))

# Remove intermediary datasets
rm(final_data, msk_nomo_data)

# Make a .csv file
write.csv(msk_wide_data, 'data/msk_wide_postop_data.csv')
