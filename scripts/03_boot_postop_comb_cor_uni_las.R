###########################################################################
# PhD - Post Op Journal Article
## LASSO Cox Bootstrapped
## mRNA data correlation pre-filter
## and univariate Feature Selection
# 03_boot_postop_comb_cor_las.R
###########################################################################

# Load the dataset
source("scripts/02_data_comb_cor_postop.R")
source("scripts/mycalibration_function2.R")
#source("scripts/mycalibration_plot_function.R")
require(RegParallel)

###########################################################################
# Prep and calculation of Apparent C-Index
###########################################################################
# Lasso
## Data preparation
x_full <- model.matrix( ~ .-1, cor_df[,c(1:8,11:ncol(cor_df))])
y_full <- Surv(time = cor_df$BCR_FreeTime,
               event = cor_df$BCR_Event)

#Feature Selection
options(scipen=10)
options(digits=6)
data <- cor_df[,11:ncol(cor_df)]
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

x_app <- model.matrix( ~ .-1, gene_data[,c(1:8,11:ncol(gene_data))])
y_app <- Surv(time = gene_data$BCR_FreeTime,
              event = gene_data$BCR_Event)

# Model w/ all 133 obs (app)
## Lasso Selection
set.seed(123)
cv_lasso_cox_app <- cv.glmnet(x = x_app,
                              y = y_app,
                              family="cox")
plot(cv_lasso_cox_app)

lasso_cox_app <- glmnet(x = x_app,
                        y = y_app,
                        family="cox",
                        lambda = cv_lasso_cox_app$lambda.1se)
pred_app <- predict(lasso_cox_app, newx = x_app)

# Discrimination
c_app <- glmnet::Cindex(pred_app, y = y_app)

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

stab.v <- rep.int(0,length(colnames(x_full)))
names(stab.v) <- colnames(x_full)
stab.uni.v <- rep.int(0,length(colnames(x_full)))
names(stab.uni.v) <- colnames(x_full)
selection.v <- data.frame(matrix(ncol = length(colnames(x_full)), nrow = 0))
colnames(selection.v) <- colnames(x_full)
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
  oob_gene_data <- oob[, c(1:10, genes_columns)]
  orig_gene_data <- cor_df[, c(1:10, genes_columns)]
  
  # Univariate Stability
  genes <- colnames(cor_df[,c(genes_columns)])
  for(k in 1:length(genes)) {
    stab.uni.v[which(names(stab.uni.v) == genes[k])] <- stab.uni.v[which(names(stab.uni.v) == genes[k])] + 1
  }
  
  # Lasso data prep
  x <- model.matrix( ~ .-1, bag[,c(1:8,11:ncol(bag_gene_data))])
  y <- Surv(time = bag_gene_data$BCR_FreeTime,
            event = bag_gene_data$BCR_Event)
  test_x <- model.matrix( ~ .-1, oob[,c(1:8,11:ncol(oob_gene_data))])
  test_y <- Surv(time = oob_gene_data$BCR_FreeTime,
                 event = oob_gene_data$BCR_Event)
  orig_x <- model.matrix( ~ .-1, orig_gene_data[,c(1:8,11:ncol(orig_gene_data))])
  orig_y <- Surv(time = orig_gene_data$BCR_FreeTime,
                 event = orig_gene_data$BCR_Event)
  
  # Lasso Modelling
  set.seed(123)
  cv_lasso_cox <- cv.glmnet(x = x,
                            y = y,
                            family="cox")
  lasso_cox <- glmnet(x = x,
                      y = y,
                      family="cox",
                      lambda = cv_lasso_cox$lambda.1se)
  mod_fits[[i]] <- lasso_cox
  train_pred <- predict(lasso_cox, newx = x)
  
  # Stability
  terms <- rownames(lasso_cox$beta)[lasso_cox$beta@i +1]
  for(j in 1:length(terms)) {
    stab.v[which(names(stab.v) == terms[j])] <- stab.v[which(names(stab.v) == terms[j])] + 1
  }
  for(j in 1:length(terms)) {
    selection.v[i,which(colnames(selection.v) == terms[j])] <- 1
  }
  
  # Validation
  oob_val <- predict(lasso_cox, newx = test_x)
  orig_val <- predict(lasso_cox, newx = orig_x)
  
  # For ROC and DCA
  preds <- c(oob_val)
  names(preds) <- rownames(oob)
  preds <- (preds - min(preds))/(max(preds)-min(preds))
  for(k in 1:length(preds)) {
    haz[i,which(colnames(haz) == names(preds)[k])] <- preds[k]
  }
  
  ## OOB 2
  oob_surv <- survfit(lasso_cox, x = x, y = y, newx = test_x)
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
  oob_surv <- survfit(lasso_cox, x = x, y = y, newx = test_x)
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
  oob_sc = aggregate(survfit(lasso_cox, x = x, y = y, newx = test_x))
  oob_cal <- as.data.frame(my.calibration2(oob_sc, oob_km, oob_Surv.output))
  oob_cal <- oob_cal %>% mutate(group = rownames(oob_cal))
  
  oob_res[(((i-1)*5)+1):(i*5),] <- oob_cal
  
  # store results
  c_oob[i] <- glmnet::Cindex(oob_val, y = test_y)
  o_boot[i] <- glmnet::Cindex(train_pred, y = y) - glmnet::Cindex(orig_val, y = y_full)
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
app_var <- rownames(lasso_cox_app$beta)[lasso_cox_app$beta@i +1]
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
write.csv(selection.v, 'data/var_select_comb_cor_uni_las.csv')
stab.v.df <- as.data.frame(stab.v[which(stab.v != 0)])
