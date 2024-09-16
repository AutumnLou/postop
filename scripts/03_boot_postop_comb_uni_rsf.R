###########################################################################
# PhD - Post Op Journal Article
## RSF Bootstrapped
## mRNA data univariate feature selection
# 03_boot_postop_comb_uni_rsf.R
###########################################################################

# Load the dataset
source("scripts/02_data_comb_postop.R")

# Load Libraries/Functions
source("scripts/mycalibration_function2.R")
source("scripts/mypredictSurvProb_functions.R")
source("scripts/RFE_for_RSF_function.R")
require(RegParallel)

# Prep date for rfe function
colnames(comb_data)[9] <- "time"
colnames(comb_data)[10] <- "status"

#Feature Selection
options(scipen=10)
options(digits=6)
data <- comb_data[,11:ncol(comb_data)]
colnames(data) <- paste0('gene', 1:ncol(data))
rownames(data) <- paste0('sample', 1:nrow(data))
variables = colnames(data)[1:ncol(data)]
data$time <- comb_data[,9]
data$alive <- comb_data[,10]

res4 <- RegParallel(
  data = data,
  formula = 'Surv(time, as.integer(alive)) ~  [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = variables,
  blocksize = 2000,
  p.adjust = "BH"
)


res5 <- res4[!is.na(res4$P),]
res5 <- res5[order(res5$LogRank.adjust, decreasing = FALSE),]
final <- head(res5, 50)
final <- subset(final, LogRank.adjust < 0.05)
probes <- gsub('[^0-9.-]', '', final$Variable)

app_genes_columns <- sort(as.numeric(probes) + 10)

gene_data <- comb_data[, c(1:10, app_genes_columns)]

# Model w/ all 133 obs (app)
set.seed(123)
rsf_app <- rsf_rfe(gene_data)
pred_app <- predict.rfsrc(rsf_app)

# Discrimination
c_app <- 1 - pred_app$err.rate[500]

###########################################################################
# Resampling Framework
###########################################################################
# 100 Bootstrap resamples
set.seed(123)
B <- 100
bootIndex <- createResample(factor(comb_data$status), times = B, list = TRUE)

c_oob <- numeric(B)
o_boot <- numeric(B)
oob_har_res_2 <- data.frame(group = numeric(B*5), mean.p = numeric(B*5),
                            km = numeric(B*5), lower = numeric(B*5),
                            upper = numeric(B*5))
oob_har_res_5 <- data.frame(group = numeric(B*5), mean.p = numeric(B*5),
                            km = numeric(B*5), lower = numeric(B*5),
                            upper = numeric(B*5))
oob_res <- data.frame(pred = numeric(B*5), obs = numeric(B*5),
                      mean.obs = numeric(B*5), low = numeric(B*5),
                      up = numeric(B*5), means.low = numeric(B*5),
                      means.up = numeric(B*5), group = numeric(B*5))


stab.v <- rep.int(0,length(colnames(comb_data[,c(1:8,11:ncol(comb_data))])))
names(stab.v) <- colnames(comb_data[,c(1:8,11:ncol(comb_data))])
stab.uni.v <- rep.int(0,length(colnames(comb_data[,c(1:8,11:ncol(comb_data))])))
names(stab.uni.v) <- colnames(comb_data[,c(1:8,11:ncol(comb_data))])
selection.v <- data.frame(matrix(ncol = length(colnames(comb_data[,c(1:8,11:ncol(comb_data))])), nrow = 0))
colnames(selection.v) <- colnames(comb_data[,c(1:8,11:ncol(comb_data))])
haz <- data.frame(matrix(ncol = 133, nrow = 0))
colnames(haz) <- rownames(comb_data)
mod_fits <- list()

for(i in 1:length(bootIndex)) {
  print(i)
  # define Training / Validation data
  bag <- comb_data[bootIndex[[i]],]
  oob <- comb_data[-bootIndex[[i]],]
  
  #Feature Selection
  options(scipen=10)
  options(digits=6)
  data <- bag[,11:ncol(bag)]
  colnames(data) <- paste0('gene', 1:ncol(data))
  rownames(data) <- paste0('sample', 1:nrow(data))
  variables = colnames(data)[1:ncol(data)]
  data$time <- bag[,9]
  data$alive <- bag[,10]
  
  res4 <- RegParallel(
    data = data,
    formula = 'Surv(time, as.integer(alive)) ~  [*]',
    FUN = function(formula, data)
      coxph(formula = formula,
            data = data,
            ties = 'breslow',
            singular.ok = TRUE),
    FUNtype = 'coxph',
    variables = variables,
    blocksize = 2000,
    p.adjust = "BH"
  )
  
  
  res5 <- res4[!is.na(res4$P),]
  res5 <- res5[order(res5$LogRank.adjust, decreasing = FALSE),]
  final <- head(res5, 50)
  probes <- gsub('[^0-9.-]', '', final$Variable)
  
  genes_columns <- sort(as.numeric(probes) + 10)
  
  bag_gene_data <- bag[, c(1:10, genes_columns)]
  
  # Univariate Stability
  genes <- colnames(comb_data[,c(genes_columns)])
  for(k in 1:length(genes)) {
    stab.uni.v[which(names(stab.uni.v) == genes[k])] <- stab.uni.v[which(names(stab.uni.v) == genes[k])] + 1
  }
  
  # RSF Selection
  set.seed(123)
  rsf_boot <- rsf_rfe(bag_gene_data)
  mod_fits[[i]] <- rsf_boot
  pred_boot <- predict.rfsrc(rsf_boot)
  
  # Stability
  terms <- names(rsf_boot$importance)
  for(j in 1:length(terms)) {
    stab.v[which(names(stab.v) == terms[j])] <- stab.v[which(names(stab.v) == terms[j])] + 1
  }
  for(j in 1:length(terms)) {
    selection.v[i,which(colnames(selection.v) == terms[j])] <- 1
  }
  
  # Validation
  oob_val <- predict.rfsrc(rsf_boot, newdata = oob)
  orig_val <- predict.rfsrc(rsf_boot, newdata = comb_data)
  
  # For ROC and DCA
  preds <- oob_val$predicted
  names(preds) <- rownames(oob)
  preds = (preds - min(preds))/(max(preds)-min(preds))
  for(k in 1:length(preds)) {
    haz[i,which(colnames(haz) == names(preds)[k])] <- preds[k]
  }
  
  ## OOB 2
  oob_surv2yr <-  pec::predictSurvProb(rsf_boot, newdata = oob, times = 24)
  oob_cal_dat <- data.frame(BCR_FreeTime = oob$BCR_FreeTime,
                            BCR_Event = oob$BCR_Event,
                            surv2yr = oob_surv2yr,
                            row.names = rownames(oob))
  oob_cal_dat <- oob_cal_dat[order(oob_cal_dat$surv2yr),]
  oob_cal_dat$group <- cut(1:nrow(oob_cal_dat), 5, labels = FALSE)
  oob_cal_dat <- oob_cal_dat %>% 
    group_by(group) %>% 
    mutate(mean.p = mean(surv2yr),
           km = summary(survfit(Surv(BCR_FreeTime,BCR_Event)~1), times = 24, extend = TRUE)$surv,
           lower = summary(survfit(Surv(BCR_FreeTime,BCR_Event)~1), times = 24, extend = TRUE)$lower,
           upper = summary(survfit(Surv(BCR_FreeTime,BCR_Event)~1), times = 24, extend = TRUE)$upper) %>%
    summarise(mean.p = mean(mean.p),
              km = mean(km),
              lower = mean(lower),
              upper = mean(upper))
  
  oob_har_res_2[(((i-1)*5)+1):(i*5),] <- oob_cal_dat
  
  ## OOB 5
  oob_surv5yr <- pec::predictSurvProb(rsf_boot, newdata = oob, times = 60)
  oob_cal_dat <- data.frame(BCR_FreeTime = oob$BCR_FreeTime,
                            BCR_Event = oob$BCR_Event,
                            surv5yr = oob_surv5yr,
                            row.names = rownames(oob))
  oob_cal_dat <- oob_cal_dat[order(oob_cal_dat$surv5yr),]
  oob_cal_dat$group <- cut(1:nrow(oob_cal_dat), 5, labels = FALSE)
  oob_cal_dat <- oob_cal_dat %>% 
    group_by(group) %>% 
    mutate(mean.p = mean(surv5yr),
           km = summary(survfit(Surv(BCR_FreeTime,BCR_Event)~1), times = 60, extend = TRUE)$surv,
           lower = summary(survfit(Surv(BCR_FreeTime,BCR_Event)~1), times = 60, extend = TRUE)$lower,
           upper = summary(survfit(Surv(BCR_FreeTime,BCR_Event)~1), times = 60, extend = TRUE)$upper) %>%
    summarise(mean.p = mean(mean.p),
              km = mean(km),
              lower = mean(lower),
              upper = mean(upper))
  
  oob_har_res_5[(((i-1)*5)+1):(i*5),] <- oob_cal_dat
  
  ## OOB
  oob_km <- survfit(Surv(oob$BCR_FreeTime,oob$BCR_Event)~1)
  oob_Surv.output = Surv(oob$BCR_FreeTime,oob$BCR_Event)
  oob_sc = aggregate.predictSurvProb(pec::predictSurvProb(rsf_boot, newdata = oob, times = oob_km$time), oob_km$time)
  oob_cal <- as.data.frame(my.calibration2(oob_sc, oob_km, oob_Surv.output))
  oob_cal <- oob_cal %>% mutate(group = rownames(oob_cal))
  
  oob_res[(((i-1)*5)+1):(i*5),] <- oob_cal
  
  # store results
  c_oob[i] <- 1 - oob_val$err.rate[500]
  o_boot[i] <- (1 - pred_boot$err.rate[500]) - (1 - orig_val$err.rate[500])
}

###########################################################################
# Results
###########################################################################
# Univariate Stability
uni_var <- colnames(comb_data[,c(app_genes_columns)])
stab.uni.v <- sort(stab.uni.v, decreasing = TRUE)
barchart(head(stab.uni.v, 10),
         xlab = "Number of selections accross resamples",
         xlim = 0:100)

# Stability
app_var <- names(rsf_app$importance)
stab.v <- sort(stab.v, decreasing = TRUE)
barchart(head(stab.v, 10),
         xlab = "Number of uses accross resamples",
         xlim = 0:100)

# Selection
selection.v <- ifelse(is.na(selection.v), 0, 1)
selection.v <- as.data.frame(selection.v)
selection.v <- selection.v %>%
  dplyr::select(which(colSums(selection.v)!=0))

# OOB Discrimination
mean(c_oob)
quantile(c_oob, probs = c(0.025, 0.975))

# Harrels Discrimination
o <- mean(o_boot)
c_app - o
quantile(c_app - o_boot, probs = c(0.025, 0.975))

# Calibration
oob_har_res_2 <- oob_har_res_2 %>%
  group_by(group) %>%
  summarise(mean.p = mean(mean.p),
            km = mean(km),
            lower = mean(lower),
            upper = mean(upper))

my.calibration.plot(oob_har_res_2$mean.p,
                    oob_har_res_2$km,
                    oob_har_res_2$lower,
                    oob_har_res_2$upper,
                    xlab="Predicted 2-year Survival",
                    ylab="Observed 2-year Survival")

oob_har_res_5 <- oob_har_res_5 %>%
  group_by(group) %>%
  summarise(mean.p = mean(mean.p),
            km = mean(km),
            lower = mean(lower),
            upper = mean(upper))

my.calibration.plot(oob_har_res_5$mean.p,
                    oob_har_res_5$km,
                    oob_har_res_5$lower,
                    oob_har_res_5$upper,
                    xlab="Predicted 5-year Survival",
                    ylab="Observed 5-year Survival")

oob_res <- oob_res %>%
  group_by(group) %>%
  summarise(pred = mean(pred, na.rm = T),
            obs = mean(obs, na.rm = T),
            mean.obs = mean(mean.obs, na.rm = T),
            low = mean(low, na.rm = T),
            up = mean(up, na.rm = T),
            means.low = mean(means.low, na.rm = T),
            means.up = mean(means.up, na.rm = T))

my.calibration.plot2(oob_res,
                     xlab="Predicted",
                     ylab="Observed")

# save discrimination data
dis_dat <- data.frame(oob = c_oob,
                      har = c_app - o_boot)
roc_dat <- data.frame(time = comb_data$time,
                      status = comb_data$status,
                      marker = colMeans(haz, na.rm = TRUE))
write.csv(selection.v, 'data/var_select_comb_uni_rsf.csv')
stab.v.df <- as.data.frame(stab.v[which(stab.v != 0)])
