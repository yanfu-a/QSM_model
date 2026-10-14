pacman::p_load(openxlsx,ppcor,data.table,tidyfst,tidyverse,plyr,stringr,vroom,writexl,vroom,lubridate,Hmisc,ppcor,plyr,survival,mgcv,purrr)
# 261014 Subject design ---------------------------------------------------
# Sesstion 0 data prepare -------------------------------------------------
## Part 1 QSM data import  -----------------------------------------------------------
dt_QSM <- fread("combined_REF416_HUASHAN418_QSM82_covariates.csv")
colnames(dt_QSM) <- gsub("-", "_", colnames(dt_QSM), fixed = TRUE)
# Sesstion 1 Surreal-GAN -------------------------------------------------
## Part 1 残差化分析 --------------------------------------------
### Part 1.1 REF组残差化分析 --------------------------------------------
dt_REF_GROUP <- dt_QSM %>% filter(model_label == "REF")

roi_cols <- dput(colnames(dt_REF_GROUP)[c(2:83)])
dt_REF_GROUP$sex <- as.factor(dt_REF_GROUP$sex)
dt_REF_final <- dt_REF_GROUP
ref_baseline_params <- list()
for (roi in roi_cols) {
  # A. 仅在 REF 组拟合线性回归模型: QSM ~ 性别 + 站点
  formula_str <- paste(roi, "~ sex + age + site")
  fit <- lm(as.formula(formula_str), data = dt_REF_GROUP)
  # B. 计算“预期偏差”并提取纯净残差
  # 预期偏差 = 模型的预测值 - 模型的截距 (截距代表纯净的物理基准线)
  expected_effect <- predict(fit, newdata = dt_REF_GROUP) - coef(fit)["(Intercept)"]
  # 残差 = 原始 QSM - 预期偏差 (抹平了性别和机器差异)
  resid_qsm <- dt_REF_GROUP[[roi]] - expected_effect
  # C. 计算该 ROI 纯净残差的 均值 和 标准差 (Z-score 基准)
  ref_mean <- mean(resid_qsm, na.rm = TRUE)
  ref_sd <- sd(resid_qsm, na.rm = TRUE)
  # D. 对 REF 组进行 Z-score 标准化
  dt_REF_final[[roi]] <- (resid_qsm - ref_mean) / ref_sd
  # E. 【关键】把这个脑区的回归模型、均值、标准差保存下来
  ref_baseline_params[[roi]] <- list(
    model = fit,
    mean = ref_mean,
    sd = ref_sd
  )
}
dt_REF_final$diagnosis <- -1

### Part 1.2 TRA组残差化分析 --------------------------------------------
dt_TRA_GROUP <- dt_QSM %>% filter(!model_label == "REF")

dt_TRA_GROUP$sex <- as.factor(dt_TRA_GROUP$sex)
# dt_MNI_TAR$p54_i2 <- as.factor(dt_MNI_TAR$p54_i2)
dt_TRA_final <- dt_TRA_GROUP
for (roi in roi_cols) {
  # A. 调出 REF 组在这个脑区的专属回归模型、均值和标准差
  ref_fit <- ref_baseline_params[[roi]]$model
  ref_mean <- ref_baseline_params[[roi]]$mean
  ref_sd <- ref_baseline_params[[roi]]$sd
  # B. 用 REF 的模型去预测 TAR 组受试者因“性别”和“机器”带来的预期偏倚量
  # 注意：这里减去截距是为了保留真实世界的物理 QSM 均值差异
  expected_effect_TRA <- predict(ref_fit, newdata = dt_TRA_GROUP) - coef(ref_fit)["(Intercept)"]
  # C. 提取 TAR 组纯净残差 (剥离了机器误差和性别差异)
  resid_qsm_TRA <- dt_TRA_GROUP[[roi]] - expected_effect_TRA
  # D. 用 REF 组的均值和标准差，对 TAR 组进行 Z-score 标准化
  dt_TRA_final[[roi]] <- (resid_qsm_TRA - ref_mean) / ref_sd
}
dt_TRA_final$diagnosis <- 1

## Part 2 Surreal-GAN learning data --------------------------------------------
dt_QSM_GROUP_S2 <- rbind(dt_REF_final, dt_TRA_final)
dt_QSM_GROUP_S2 <- dt_QSM_GROUP_S2 %>% filter(!is.na(L_bankssts))
colnames(dt_QSM_GROUP_S2)[1] <- "participant_id"
dt_QSM_GROUP_S2 <- dt_QSM_GROUP_S2[,c(1,94,2:83)]
fwrite(dt_QSM_GROUP_S2,"dt_Model_learning_1014.csv",row.names = FALSE)