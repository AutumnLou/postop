###########################################################################
# MSc
# 02_data_msc.R
###########################################################################

# Read in the required packages ------------------------------------------

source("scripts/01_data_clean.R")

# Create a dataset for use with the MSc model
msc_data <- final_data

# Create new variables
msc_data <- msc_data %>%
  mutate(bgg1 = ifelse(BxGG1 == 3,
                       '3', '4p'),
         pgg1 = ifelse(PathGG1 == 3,
                       '3', '4p'),
         clin_stage = ifelse(ClinT_Stage == 'T1C',
                             'T1',
                             ifelse(ClinT_Stage == 'T2' |
                                      ClinT_Stage == 'T2A' |
                                      ClinT_Stage == 'T2B'|
                                      ClinT_Stage == 'T2C',
                                    'T2', 'T3')),
         path_stage = ifelse(PathStage == 'T2A' |
                               PathStage == 'T2B' |
                               PathStage == 'T2C',
                             'T2',
                             ifelse(PathStage == 'T3A' |
                                      PathStage == 'T3B' |
                                      PathStage == 'T3C',
                                    'T3', 'T4')))
# Correct variable types
msc_data$bgg1 <- as.factor(msc_data$bgg1)
msc_data$pgg1 <- as.factor(msc_data$pgg1)
msc_data$clin_stage <- as.factor(msc_data$clin_stage)
msc_data$path_stage <- as.factor(msc_data$path_stage)

# Remove the original columns
msc_data <- msc_data %>%
  select(!c(BxGG1, PathGG1, ClinT_Stage, PathStage))


# Remove intermediary datasets
rm(final_data)

# Removing BGG1, BGG2, BGGS PGG1 and PGG2
# clin_stage, path_stage and B_PSA from post-op
msc_data <- msc_data %>%
  
  select(!c(bgg1, pgg1, BxGG2, PathGG2,
            BxGGS, clin_stage, path_stage,
            PreDxBxPSA))

# Make a .csv file
write.csv(msc_data, 'data/msc_data.csv') 
