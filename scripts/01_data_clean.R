###########################################################################

# 01_data_clean.R
###########################################################################

# Read in the required packages ------------------------------------------

source("scripts/00_packages.R")

# Load the data
data <- read.csv("data/raw_data.csv", header = TRUE)

# Attach the data
attach(data)

# Make the row names the Sample ID
rownames(data) <- data$Sample.ID
data <- data %>% select(!Sample.ID)

# Reduce number of variables
bcr_data <- data %>%
  # Remove all rows without BCR_Free Time
  filter(!is.na(BCR_FreeTime) & !is.na(BCR_Event)) %>%
  
  # Remove MetSite and Type
  select(!c(MetSite, Type)) %>%
  
  # Remove ERG columns
  select(!c(ERG.fusion.aCGH, ERG.fusion.gex)) %>%
  
  # Removed Num_Nodes_Removed and Num_Nodes_Positive
  select(!c(Num_Nodes_Removed, Num_Nodes_Positive)) %>%
  
  # Removed Nomogram and Copy Number Cluster Columns
  select(!c(Nomogram.NomoPred_ECE, Nomogram.NomoPred_LNI,
            Nomogram.NomoPred_OCD, Nomogram.NomoPred_SVI,
            Nomogram.PFP_PostRP, Copy.number.Cluster)) %>%
  
  # Remove RP type
  select(!c(RP_Type)) %>%
  
  # Remove MetsEvent, SurvTime and Event
  ## as they are post BCR or censored
  select(!c(MetsEvent, SurvTime, Event))


# Remove therapy columns
bcr_data <- bcr_data %>%
  select(!c(NeoAdjRadTx, ChemoTx, HormTx, RadTxType))

# Remove incomplete observations
bcr_comp_data <- na.omit(bcr_data)

# Correct variable types
bcr_comp_data$ClinT_Stage <- as.factor(bcr_comp_data$ClinT_Stage)
bcr_comp_data$SMS <- as.factor(bcr_comp_data$SMS)
bcr_comp_data$SVI <- as.factor(bcr_comp_data$SVI)
bcr_comp_data$LNI <- as.factor(bcr_comp_data$LNI)
bcr_comp_data$PathStage <- as.factor(bcr_comp_data$PathStage)
bcr_comp_data$PathGG1 <- as.factor(bcr_comp_data$PathGG1)
bcr_comp_data$PathGG2 <- as.factor(bcr_comp_data$PathGG2)
bcr_comp_data$PathGGS <- as.factor(bcr_comp_data$PathGGS)

# Change the categories of ECE column
bcr_comp_data$ECE <- ifelse(bcr_comp_data$ECE == "NONE",
                            "Absent","Present")
bcr_comp_data$ECE <- as.factor(bcr_comp_data$ECE)


# Change the categories of the Race column
bcr_comp_data$Race <-
  mapvalues(bcr_comp_data$Race, from = c("Black Non Hispanic", "White Non Hispanic",
                                    "Black Hispanic", "White Hispanic", "Asian",
                                    "Unknown"),
            to = c("2Black", "1White", "2Black", "1White", "3Asian", "4Unknown"))
bcr_comp_data$Race <- as.factor(bcr_comp_data$Race)

# Change levels of Biopsy Gleason score lowest == 3
## and sum == 6 
bcr_comp_data$BxGG1 <- mapvalues(bcr_comp_data$BxGG1,
                                 from = c(2, 3, 4, 5),
                                 to = c(3, 3, 4, 5))
bcr_comp_data$BxGG2 <- mapvalues(bcr_comp_data$BxGG2,
                                 from = c(2, 3, 4, 5),
                                 to = c(3, 3, 4, 5))
bcr_comp_data$BxGGS <- mapvalues(bcr_comp_data$BxGGS,
                                 from = c(5, 6, 7, 8, 9),
                                 to = c(6, 6, 7, 8, 9))
bcr_comp_data$BxGG1 <- as.factor(bcr_comp_data$BxGG1)
bcr_comp_data$BxGG2 <- as.factor(bcr_comp_data$BxGG2)
bcr_comp_data$BxGGS <- as.factor(bcr_comp_data$BxGGS)

# Change the BCR event to a numeric variable
bcr_comp_data$BCR_Event <- ifelse(bcr_comp_data$BCR_Event == "NO",
                                  0, 1)

## Drop Observation with PSA 506 (error) and Drop MET duplicate
final_data <- bcr_comp_data[!(row.names(bcr_comp_data)
                                 %in% c('PCA0045','PCA0187', 'PCA0207')), ]

# Remove intermediary datasets
rm(data, bcr_data, bcr_comp_data)

# Make a .csv file
write.csv(final_data, 'data/clean_bcr_data.csv')
