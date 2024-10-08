###########################################################################
# PhD - Post Op Journal Article
## Traditional Cox Bootstrapped
## mRNA data correlation pre-filter
## and univariate feature selection
# 03_boot_postop_comb_cor_uni_cox.R
###########################################################################

# Load the dataset
source("scripts/02_data_comb_cor_postop.R")
source("scripts/mycalibration_function2.R")
#source("scripts/mycalibration_plot_function.R")
require(RegParallel)


###########################################################################
# Prep and calculation of Apparent C-Index
###########################################################################

#Feature Selection
options(scipen=10)
options(digits=6)
data <- cor_df[,11:2594]
colnames(data) <- paste0('gene', 1:ncol(data))
rownames(data) <- paste0('sample', 1:nrow(data))
variables = colnames(data)[1:ncol(data)]
data$time <- cor_df[,9]
data$alive <- cor_df[,10]

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

gene_data <- cor_df[, c(1:10, app_genes_columns)]

# Model w/ all 133 obs (app)
app_start_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ 1, data = gene_data, x = TRUE)
app_full_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ ., data = gene_data, x = TRUE)
app_fit_step <- MASS::stepAIC(app_start_cox, direction = "both",
                     scope = app_full_cox$formula,
                     steps = 10)

# Discrimination
c_app <- concordance(app_fit_step)$concordance

###########################################################################
# Resampling Framework
###########################################################################
# 100 Bootstrap resamples
set.seed(123)
B <- 100
bootIndex <- createResample(factor(cor_df$BCR_Event), times = B, list = TRUE)

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
stab.v <- rep.int(0,length(colnames(cor_df[,c(1:8,11:ncol(cor_df))])))
names(stab.v) <- colnames(cor_df[,c(1:8,11:ncol(cor_df))])
stab.uni.v <- rep.int(0,length(colnames(cor_df[,c(1:8,11:ncol(cor_df))])))
names(stab.uni.v) <- colnames(cor_df[,c(1:8,11:ncol(cor_df))])
selection.v <- data.frame(matrix(ncol = length(colnames(cor_df[,c(1:8,11:ncol(cor_df))])), nrow = 0))
colnames(selection.v) <- colnames(cor_df[,c(1:8,11:ncol(cor_df))])
haz <- data.frame(matrix(ncol = 133, nrow = 0))
colnames(haz) <- rownames(cor_df)
mod_fits <- list()

for(i in 1:length(bootIndex)) {
  print(i)
  # define Training / Validation data
  bag <- cor_df[bootIndex[[i]],]
  oob <- cor_df[-bootIndex[[i]],]
  
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
  final <- subset(final, LogRank.adjust < 0.05)
  probes <- gsub('[^0-9.-]', '', final$Variable)
  
  genes_columns <- sort(as.numeric(probes) + 10)
  
  bag_gene_data <- bag[, c(1:10, genes_columns)]
  
  # Univariate Stability
  genes <- colnames(cor_df[,c(genes_columns)])
  for(k in 1:length(genes)) {
    stab.uni.v[which(names(stab.uni.v) == genes[k])] <- stab.uni.v[which(names(stab.uni.v) == genes[k])] + 1
  }
  
  # Stepwise Selection
  start_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ 1, data = bag_gene_data, x = TRUE)
  
  full_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ ., data = bag_gene_data, x = TRUE)
  
  fit_step <- MASS::stepAIC(start_cox, direction = "both", scope = full_cox$formula,
                   trace = 0, steps = 10)
  mod_fits[[i]] <- fit_step
  
  # Stability
  terms <- attr(terms(fit_step), "term.labels")
  for(j in 1:length(terms)) {
    stab.v[which(names(stab.v) == terms[j])] <- stab.v[which(names(stab.v) == terms[j])] + 1
  }
  for(j in 1:length(terms)) {
    selection.v[i,which(colnames(selection.v) == terms[j])] <- 1
  }
  
  # Validation
  oob_val <- concordance(fit_step, newdata = oob)
  orig_val <- concordance(fit_step, newdata = cor_df)
  
  # For ROC and DCA
  preds <- predict(fit_step, newdata = oob)
  preds <- (preds - min(preds))/(max(preds)-min(preds))
  for(k in 1:length(preds)) {
    haz[i,which(colnames(haz) == names(preds)[k])] <- preds[k]
  }
  
  ## OOB 2
  oob_surv <- survfit(fit_step, newdata = oob)
  oob_surv_time <- which(abs(oob_surv$time-24)==min(abs(oob_surv$time-24)))
  oob_surv2yr <- oob_surv$surv[oob_surv_time,]
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
  oob_surv <- survfit(fit_step, newdata = oob)
  oob_surv_time <- which(abs(oob_surv$time-60)==min(abs(oob_surv$time-60)))
  oob_surv5yr <- oob_surv$surv[oob_surv_time,]
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
  oob_sc = aggregate(survfit(fit_step, newdata = oob))
  oob_cal <- as.data.frame(my.calibration2(oob_sc, oob_km, oob_Surv.output))
  oob_cal <- oob_cal %>% mutate(group = rownames(oob_cal))
  
  oob_res[(((i-1)*5)+1):(i*5),] <- oob_cal
  
  # store results
  c_oob[i] <- oob_val$concordance
  o_boot[i] <- concordance(fit_step)$concordance - orig_val$concordance
 }

###########################################################################
# Results
###########################################################################

# Univariate Stability
uni_var <- colnames(cor_df[,c(app_genes_columns)])
stab.uni.v <- sort(stab.uni.v, decreasing = TRUE)
barchart(head(stab.uni.v, 10),
         xlab = "Number of selections accross resamples",
         xlim = 0:100)

# Stability
app_var <- attr(terms(app_fit_step), "term.labels")
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
roc_dat <- data.frame(time = cor_df$BCR_FreeTime,
                      status = cor_df$BCR_Event,
                      marker = colMeans(haz, na.rm = TRUE))
write.csv(selection.v, 'data/var_select_comb_cor_uni_cox.csv')
stab.v.df <- as.data.frame(stab.v[which(stab.v != 0)])
