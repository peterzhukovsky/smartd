library(haven)
library(glmnet)
library(pROC)
library(readxl)
library(dplyr)

# marker building
dataset <- read_sav("EMBARC_BUP_paper_data_6August2019_forJASP_N87_Analysis_YA.sav")
head(dataset)

# Buproprion predictions
data.bup = dataset[dataset$Stage2TX=="BUP",] %>% 
  select(Phase2_responder, resp_bias_tot_sess1, reward_sensitivity, drugXresp_biNACCseed_rACCeffect_p001FDR)
data.bup = data.bup[complete.cases(data.bup),]

# loo with three different methods
rf_loo = sapply(1:nrow(data.bup), function(i){
  rf = randomForest::randomForest(as.factor(Phase2_responder) ~ resp_bias_tot_sess1 + reward_sensitivity + drugXresp_biNACCseed_rACCeffect_p001FDR, data = data.bup[-i,])
  predict = predict(rf, data.bup[i,,drop = F], type = "prob")[2]
  c(data.bup$Phase2_responder[i], predict)
})

glm_loo = sapply(1:nrow(data.bup), function(i){
  fit = glm(Phase2_responder ~ resp_bias_tot_sess1 + reward_sensitivity + drugXresp_biNACCseed_rACCeffect_p001FDR, family = "binomial", data = data.bup[-i,])
  predict = predict(fit, data.bup[i,,drop = F], type = "response")
  c(data.bup$Phase2_responder[i], predict)
})

lasso_loo = sapply(1:nrow(data.bup), function(i){
  X = as.matrix(data.bup[-i, 2:4])
  y = data.bup$Phase2_responder[-i]
  fit = cv.glmnet(x = X, y = y, family = "binomial")
  predict = predict(fit, as.matrix(data.bup[i, 2:4, drop = F]), type = "response", s = "lambda.min")
  c(data.bup$Phase2_responder[i], predict)
})

rf_roc = roc(rf_loo[1,], rf_loo[2,])
glm_roc = roc(glm_loo[1,], glm_loo[2,])
lasso_roc = roc(lasso_loo[1,], lasso_loo[2,])

# 
bup_roc_plt = data.frame( sen = c(rf_roc$sensitivities,
                                  glm_roc$sensitivities,
                                  lasso_roc$sensitivities),
                          fpr = 1 - c(rf_roc$specificities,
                                  glm_roc$specificities,
                                  lasso_roc$specificities),
                          method = rep(c("RF", "GLM", "LASSO"),
                                       c(length(rf_roc$sensitivities), 
                                         length(glm_roc$sensitivities),
                                         length(lasso_roc$sensitivities))),
                          marker = "BUP") %>%
  arrange(fpr, sen)

# need to use glm based on the loo res
bup_thresh = glm_roc$thresholds[order( (1-glm_roc$sensitivities)^2 + (1-glm_roc$specificities)^2 )[1]]
bup_model = glm(Phase2_responder ~ resp_bias_tot_sess1 + reward_sensitivity + drugXresp_biNACCseed_rACCeffect_p001FDR, 
                family = "binomial", data = data.bup)

# Sertraline predictions, using Phase1 data
data.phase1 = read.csv("EMBARC_imputed_data.csv")
names(data.phase1)[1] = "ID"

data.phase1.fmri = read_excel("embDB_MDD_and_HC_20191012.xlsx")
data.phase1.fmri = data.phase1.fmri%>%select(Subject, biCC_w_rACC)
names(data.phase1.fmri) = c("ID", "drugXresp_biNACCseed_rACCeffect_p001FDR")
data.phase1.merged = left_join(data.phase1, data.phase1.fmri)
data.phase1.merged$responder =
  1*((data.phase1.merged$w8_score_17 - data.phase1.merged$w0_score_17)/
       data.phase1.merged$w0_score_17 < -0.5)
data.phase1.m.ser = data.phase1.merged %>% filter(tx ==1) %>% select(ID, flanker_acc, demo_employ_status, drugXresp_biNACCseed_rACCeffect_p001FDR,
                                                                     responder, w0_score_17, neo2_score_ne)
data.phase1.m.ser = data.phase1.m.ser[complete.cases(data.phase1.m.ser),]

# 10 fold CV
set.seed(3)
all_res = lapply(1:1e3, function(x){
  folds = caret::createFolds(data.phase1.m.ser$responder)
  cv_ser_rf = do.call(rbind, lapply(folds, function(fold){
    fit = randomForest::randomForest(as.factor(responder) ~ flanker_acc+demo_employ_status+drugXresp_biNACCseed_rACCeffect_p001FDR+
                                       neo2_score_ne+w0_score_17, data = data.phase1.m.ser[-fold,])
    predict = predict(fit, data.phase1.m.ser[fold,], type = "prob")[,2]
    cbind(data.phase1.m.ser$responder[fold], predict)
  }))
  
  cv_ser_glm = do.call(rbind, lapply(folds, function(fold){
    fit = glm(responder ~ flanker_acc+demo_employ_status+drugXresp_biNACCseed_rACCeffect_p001FDR+
                neo2_score_ne+w0_score_17,
              family = "binomial", data = data.phase1.m.ser[-fold,])
    predict = predict(fit, data.phase1.m.ser[fold,], type = "response")
    cbind(data.phase1.m.ser$responder[fold], predict)
  }))
  
  cv_ser_lasso = do.call(rbind, lapply(folds, function(fold){
    X = as.matrix(data.phase1.m.ser %>% select(-responder))[-fold,]
    y = data.phase1.m.ser$responder[-fold]
    fit = cv.glmnet(x = X, y = y, family = "binomial")
    predict = predict(fit, as.matrix(data.phase1.m.ser %>% select(-responder))[fold,], type = "response", s = "lambda.min")
    cbind(data.phase1.m.ser$responder[fold], predict)
  }))
  
  ser_rf_roc = roc(cv_ser_rf[,1], cv_ser_rf[,2], direction = "<")
  ser_glm_roc = roc(cv_ser_glm[,1], cv_ser_glm[,2], direction = "<")
  ser_lasso_roc = roc(cv_ser_lasso[,1], cv_ser_lasso[,2], direction = "<")
  list((ser_rf_roc), (ser_glm_roc), (ser_lasso_roc))
})

# average ROC
average_pROC_roc <- function(roc_list,
                             grid_length = 500,
                             conf_level = 0.95) {
  # browser()
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required but not installed.")
  }
  
  # Basic checks
  if (!is.list(roc_list) || length(roc_list) < 2) {
    stop("roc_list must be a list of at least 2 pROC::roc objects.")
  }
  
  # Common FPR grid
  fpr_grid <- seq(0, 1, length.out = grid_length)
  
  # Interpolate TPR for each ROC curve on the common FPR grid
  tpr_mat <- sapply(roc_list, function(roc_obj) {
    pROC::coords(
      roc_obj,
      x        = fpr_grid,
      input    = "fpr",
      ret      = "tpr",
      # make sure we get a simple numeric vector
      transpose = FALSE
    )$tpr
  })
  
  # Ensure matrix structure
  tpr_mat <- as.matrix(tpr_mat)
  
  # Mean TPR and SE across ROC curves
  mean_tpr <- rowMeans(tpr_mat, na.rm = TRUE)
  se_tpr   <- apply(tpr_mat, 1, sd, na.rm = TRUE) / sqrt(ncol(tpr_mat))
  
  # Pointwise confidence intervals
  alpha <- 1 - conf_level
  zval  <- qnorm(1 - alpha / 2)
  lower_tpr <- mean_tpr - zval * se_tpr
  upper_tpr <- mean_tpr + zval * se_tpr
  
  # Clamp to [0, 1]
  lower_tpr <- pmax(lower_tpr, 0)
  upper_tpr <- pmin(upper_tpr, 1)
  
  # AUC summary
  aucs      <- sapply(roc_list, function(r) as.numeric(r$auc))
  auc_mean  <- mean(aucs, na.rm = TRUE)
  auc_sd    <- sd(aucs,   na.rm = TRUE)
  
  # Return everything useful
  list(
    fpr_grid  = fpr_grid,
    mean_tpr  = mean_tpr,
    lower_tpr = lower_tpr,
    upper_tpr = upper_tpr,
    tpr_mat   = tpr_mat,
    aucs      = aucs,
    auc_mean  = auc_mean,
    auc_sd    = auc_sd,
    conf_level = conf_level
  )
}

rf_roc_all = lapply(all_res, function(x) x[[1]])
glm_roc_all = lapply(all_res, function(x) x[[2]])
lasso_roc_all = lapply(all_res, function(x) x[[3]])

rf_roc_mean = average_pROC_roc(rf_roc_all)
glm_roc_mean = average_pROC_roc(glm_roc_all)
lasso_roc_mean = average_pROC_roc(lasso_roc_all)

ser_auc_curve = data.frame(sen = c(rf_roc_mean$mean_tpr,
                                   glm_roc_mean$mean_tpr,
                                   lasso_roc_mean$mean_tpr),
                           fpr = c(rf_roc_mean$fpr_grid,
                                   glm_roc_mean$fpr_grid,
                                   lasso_roc_mean$fpr_grid),
                           method = rep(c("RF", "GLM", "LASSO"),
                                        c(length(rf_roc_mean$mean_tpr),
                                          length(glm_roc_mean$mean_tpr),
                                          length(lasso_roc_mean$mean_tpr))),
                           marker = "SER") %>%
  arrange(fpr, sen)

# plot all ROC curves
auc_curves_all = rbind(ser_auc_curve, bup_roc_plt)
p1 = ggplot(data = auc_curves_all,
       aes(x = fpr, y = sen, color = method)) + geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) + theme_bw() +
  facet_grid(.~marker) +
  xlab('1 - Specificity') + ylab("Sensitivity") +
  theme(axis.text = element_text(size = 12, color = "black"), 
        axis.title = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

# based on AUC, we use glm
threshold_all = sapply(all_res, function(x){
  roc_res = x[[2]]
  roc_res$thresholds[order((1-roc_res$sensitivities)^2 + (1-roc_res$specificities)^2)[1]]
  })
ser_threshold = mean(threshold_all)
ser_model = glm(responder ~ flanker_acc+demo_employ_status+drugXresp_biNACCseed_rACCeffect_p001FDR+
                  neo2_score_ne+w0_score_17, family = "binomial", data = data.phase1.m.ser)

# all subjects with both phase 1 and phase 2 data
# Marker status composition
data_merge = left_join( dataset, data.phase1 %>% select(ProjectSpecificID, flanker_acc, demo_employ_status, neo2_score_ne, w0_score_17),
                        by = c("ID" = "ProjectSpecificID"))
data_merge_use = data_merge %>% select(ID, flanker_acc, demo_employ_status, neo2_score_ne, w0_score_17, 
                                       drugXresp_biNACCseed_rACCeffect_p001FDR, resp_bias_tot_sess1, reward_sensitivity)

# BUP status
p_ser_all = predict(ser_model, data_merge_use, type = "response")
p_bup_all = predict(bup_model, data_merge_use, type = "response")

ser_status = p_ser_all > ser_threshold
bup_status = p_bup_all > bup_thresh

table(SER = ser_status, BUP = bup_status)

# SMART-D response analysis
data_raw = read.csv("SMART_MADRS_SHAPS_Aug25.csv")
data_raw = data_raw %>% filter(X != "")
data_trt = read.csv("trt_assignment_090425.csv")
data_assign = read_xlsx("SMARTDataTx_9.4.25.xlsx")
colnames(data_assign) = c("ID", "Trt", "Dose")
data_assign$Trt = ifelse(data_assign$Trt == "BUPROPION", 1, 0)

data_combined = left_join(data_raw, data_trt, by = c("screening_visit_arm_1.1.record_id" = "ID"))
data_combined = left_join(data_combined, data_assign, by = c("screening_visit_arm_1.1.record_id" = "ID") )
# consistency
data_combined = data_combined %>% mutate(con = (BUP == 1 & Trt == 1) | (SER == 1 & Trt == 0))

library(lme4)
library(lmerTest)

data_c_long = data_combined %>% pivot_longer( cols = madrs_w0:madrs_w8, names_to = "time", values_to = "val") %>%
  mutate(time = as.numeric(gsub("madrs_w", "", time))) %>%
  mutate(group = as.factor(SER + BUP))

# comparing two trt
lmm1 = lmer(val ~ time*Trt + sex + age + (time|X), data = data_c_long, REML = F)
anova(lmm1)

# comparing different -- vs. +- vs ++
lmm2 = lmer(val ~ sex + age + time*group + (time|X), data = data_c_long, REML = F)
anova(lmm2)

# comparing consistent vs. inconsistent
lmm3 = lmer(val ~ sex + age + time*con + (time|X), data = data_c_long, REML = F)
anova(lmm3)

# all possible subgroups 
data_c_long = data_c_long %>% mutate(group_full = case_when( BUP == 0 & SER == 0 ~ "BUP-/SER-",
                                                             BUP == 0 & SER == 1 ~ "BUP-/SER+",
                                                             BUP == 1 & SER == 0 ~ "BUP+/SER-",
                                                             BUP == 1 & SER == 1 ~ "BUP+/SER+"))

lmm4 = lmer(val ~ sex + age + time*group_full + (time|X), data = data_c_long, REML = F)
anova(lmm4)
sjPlot::tab_model(lmm4, show.intercept = F)