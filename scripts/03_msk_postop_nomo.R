###########################################################################
# MSKCC
# 03_msk_postop_nomo.R
###########################################################################

# Load the dataset
source("scripts/02_data_msk_postop.R")
attach(msk_wide_data)

# Load all coefficients
intercept <- 5.92175832
age <- 0.00709393
psa <- -0.30482286
psa_s1 <- 0.00267446
psa_s2 <- -0.00737241
bgg2 <- -0.45903283
bgg3 <- -0.8199369
bgg4 <- -0.84800196
bgg5 <- -0.57541563
pgg2 <- -0.67641121
pgg3 <- -1.56956461
pgg4 <- -2.0580583
pgg5 <- -2.19510654
cs_2a <- -0.2453267
cs_2b <- -0.47240073
cs_2c <- -0.31475176
cs_3 <- -0.45849604
ece_coef <- -0.61600111
svi_coef <- -0.47509077
lni_coef <- -1.10568768
sms_coef <- -0.92186578

# Load scaling parameter
scaling_param <- 0.97139714

# c_index <- 0.84057261
# model_no <- 12475

# Load Spline Knots
PSAPreopKnot1 <- 0.2
PSAPreopKnot2 <- 4.71
PSAPreopKnot3 <- 7.22
PSAPreopKnot4 <- 96.53

# Create the spline variables
msk_wide_spline_data <- msk_wide_data %>%
  mutate(psa_k1 = PreDxBxPSA - PSAPreopKnot1,
         psa_k2 = PreDxBxPSA - PSAPreopKnot2,
         psa_k3 = PreDxBxPSA - PSAPreopKnot3,
         psa_k4 = PreDxBxPSA - PSAPreopKnot4) %>%
  mutate(sp1var = pmax(psa_k1, 0)**3 -
           (pmax(psa_k3,0)**3)*(PSAPreopKnot4 - PSAPreopKnot1)/(PSAPreopKnot4 - PSAPreopKnot3) +
           (pmax(psa_k4,0)**3)*(PSAPreopKnot3 - PSAPreopKnot1)/(PSAPreopKnot4 - PSAPreopKnot3),
         sp2var = (pmax((psa_k2),0))**3 -
           ((pmax((psa_k3),0))**3)*((PSAPreopKnot4 - PSAPreopKnot2)/(PSAPreopKnot4 - PSAPreopKnot3)) +
           ((pmax((psa_k4),0))**3)*((PSAPreopKnot3 - PSAPreopKnot2)/(PSAPreopKnot4 - PSAPreopKnot3))
  )

# Reattach data with additional variables
attach(msk_wide_spline_data)

# Build the survival model
surv_model <- intercept + age*DxAge + psa*PreDxBxPSA + psa_s1*sp1var +
  psa_s2*sp2var + bgg2*biop_grp_2 + bgg3*biop_grp_3 + bgg4*biop_grp_4 +
  bgg5*biop_grp_5 + pgg2*biop_grp_2 + pgg3*biop_grp_3 + pgg4*biop_grp_4 +
  pgg5*biop_grp_5 + cs_2a*stage_2A + cs_2b*stage_2B + cs_2c*stage_2C +
  cs_3*stage_3p + ece_coef*ece + svi_coef*svi + sms_coef*sms + lni_coef*lni

# Calculate the prediction
pred_prob_2 <- 1/(1 + (exp(-surv_model)*2)**(1/scaling_param))
pred_prob_5 <- 1/(1 + (exp(-surv_model)*5)**(1/scaling_param))

# Add predictions to dataset
msk_pred_data <- msk_wide_data %>%
  select(BCR_FreeTime, BCR_Event)
msk_pred_data$surv2 <- pred_prob_2
msk_pred_data$surv5 <- pred_prob_5
msk_pred_data$risk <- surv_model


# c-index
c_index_2 <- survcomp::concordance.index(x = 1 - msk_pred_data$surv2,
                             surv.time = msk_pred_data$BCR_FreeTime,
                             surv.event = msk_pred_data$BCR_Event)
c_index_2$c.index

c_index_5 <- survcomp::concordance.index(x = 1 - msk_pred_data$surv5,
                               surv.time = msk_pred_data$BCR_FreeTime,
                               surv.event = msk_pred_data$BCR_Event)
c_index_5$c.index
## 0.7709677

# Make a .csv file
write.csv(msk_pred_data, 'data/msk_pred_data.csv', row.names = FALSE)

