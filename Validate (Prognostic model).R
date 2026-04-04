# =============================================================================
# Section 1: Load Libraries
# =============================================================================
{
  library(survival)
  library(survivalROC)
  library(dplyr)
  library(ggplot2)
  library(caret)
  library(tidyverse)
  library(readxl)
  library(glmnet)
  library(patchwork)
}
# =============================================================================
# Section 2: Load and Filter Data
# =============================================================================
TCGA_data <- read_excel("TCGA-data.xlsx")

# Keep squamous cell carcinoma only
TCGA_data <- TCGA_data[
  grep(
    "squamous cell carcinoma",
    TCGA_data$primary_diagnosis,
    ignore.case = TRUE
  ),
]

# Exclude Stage IVB (distant metastasis)
TCGA_data <- TCGA_data[
  !grepl("IVB", TCGA_data$figo_stage, ignore.case = TRUE),
]

TCGA_data <- TCGA_data %>%
  dplyr::select(
    id = submitter_id,
    stage = figo_stage,
    days_to_last_follow_up,
    days_to_death,
    status = vital_status,
    age = age_at_index,
    Pharmacy_treat = treatments_pharmaceutical_treatment_or_therapy,
    Radiation_treat = treatments_radiation_treatment_or_therapy,
    LAMA4,
    TNFAIP6,
    SULF1,
    ADAM12,
    VCAM1,
    CTHRC1,
    GREM1,
    COL5A2,
    SPARCL1
  )

# =============================================================================
# Section 3: Data Preparation
# =============================================================================
TCGA_data <- TCGA_data %>%
  transmute(
    id,
    stage,
    stage_gr = stage %>%
      str_extract('[IV]+') %>%
      factor(levels = c("I", "II", "III", "IV")),
    across(
      c(
        LAMA4,
        SULF1,
        TNFAIP6,
        ADAM12,
        VCAM1,
        CTHRC1,
        GREM1,
        COL5A2,
        SPARCL1,
        age
      ),
      as.numeric
    ),
    across(c(Pharmacy_treat, Radiation_treat), as.factor),
    duration = coalesce(days_to_death, days_to_last_follow_up) / 365.25,
    status
  ) %>%
  mutate(across(
    c(LAMA4, SULF1, TNFAIP6, ADAM12, VCAM1, CTHRC1, GREM1, COL5A2, SPARCL1),
    ~ log2(. + 1)
  )) %>%
  filter(duration > 0, !is.na(stage_gr))

TCGA_data$status <- ifelse(TCGA_data$status == "Dead", 1, 0)

# =============================================================================
# Section 4: Train-Test Split (70/30)
# =============================================================================
set.seed(123)
split_index <- createDataPartition(
  as.factor(TCGA_data$status),
  times = 1,
  p = 0.7,
  list = FALSE
)

train_df <- TCGA_data[split_index, ]
test_df <- TCGA_data[-split_index, ]

{
  cat("Training set:", nrow(train_df), "patients\n")
  cat("Test set:    ", nrow(test_df), "patients\n")
  cat("Event rate - Train:", mean(train_df$status), "\n")
  cat("Event rate - Test: ", mean(test_df$status), "\n")
}
# =============================================================================
# Section 5: Fit Cox Models
# =============================================================================
genes <- c(
  "LAMA4",
  "SULF1",
  "TNFAIP6",
  "ADAM12",
  "VCAM1",
  "CTHRC1",
  "GREM1",
  "COL5A2",
  "SPARCL1"
)

# Survival outcome
y_train <- Surv(train_df$duration, train_df$status)

# Model matrices — built from training data; test matrix aligned to same columns
f_CAFs <- ~ LAMA4 +
  SULF1 +
  TNFAIP6 +
  ADAM12 +
  VCAM1 +
  CTHRC1 +
  GREM1 +
  COL5A2 +
  SPARCL1 +
  0
f_full <- ~ stage_gr +
  age +
  Pharmacy_treat +
  Radiation_treat +
  LAMA4 +
  SULF1 +
  TNFAIP6 +
  ADAM12 +
  VCAM1 +
  CTHRC1 +
  GREM1 +
  COL5A2 +
  SPARCL1 +
  0

X_train_CAFs <- model.matrix(f_CAFs, data = train_df)
X_train_full <- model.matrix(f_full, data = train_df)

# Align test matrices to training columns (handles missing factor levels in test)
X_test_CAFs <- model.matrix(f_CAFs, data = test_df)[, colnames(X_train_CAFs)]
X_test_full <- model.matrix(f_full, data = test_df)[, colnames(X_train_full)]

# penalty.factor: 1 = penalized (genes), 0 = unpenalized (clinical)
# derived from column names — robust to factor expansion
pf_full <- ifelse(colnames(X_train_full) %in% genes, 1, 0)

# Model matrices for clinical-only ridge
f_clin <- ~ stage_gr + age + Pharmacy_treat + Radiation_treat + 0
X_train_clin <- model.matrix(f_clin, data = train_df)
X_test_clin <- model.matrix(f_clin, data = test_df)[, colnames(X_train_clin)]

# Model 1: Clinical only — ridge Cox
set.seed(123)
cv_clinical <- cv.glmnet(
  X_train_clin,
  y_train,
  family = "cox",
  alpha = 0,
  nfolds = 5
)
ridge_clinical <- glmnet(
  X_train_clin,
  y_train,
  family = "cox",
  alpha = 0,
  lambda = cv_clinical$lambda.min
)

# Model 2: CAFs only — ridge Cox (alpha = 0), lambda by 5-fold CV
set.seed(123)
cv_CAFs <- cv.glmnet(
  X_train_CAFs,
  y_train,
  family = "cox",
  alpha = 0,
  nfolds = 5
)
ridge_CAFs <- glmnet(
  X_train_CAFs,
  y_train,
  family = "cox",
  alpha = 0,
  lambda = cv_CAFs$lambda.min
)

# Model 3: Clinical + CAFs — clinical unpenalized, genes penalized
set.seed(123)
cv_full <- cv.glmnet(
  X_train_full,
  y_train,
  family = "cox",
  alpha = 0,
  nfolds = 5,
  penalty.factor = pf_full
)
ridge_full <- glmnet(
  X_train_full,
  y_train,
  family = "cox",
  alpha = 0,
  lambda = cv_full$lambda.min,
  penalty.factor = pf_full
)

# =============================================================================
# Section 6: Cox PH Assumption Diagnostics
# Note: cox.zph requires a coxph object; fit a standard Cox for diagnostics only
# =============================================================================
cox_clinical_diag <- coxph(
  Surv(duration, status) ~ stage_gr + age + Pharmacy_treat + Radiation_treat,
  data = train_df
)
zph_clinical <- cox.zph(cox_clinical_diag)
print(zph_clinical)
par(mfrow = c(2, 2))
plot(zph_clinical)
par(mfrow = c(1, 1))

# =============================================================================
# Section 7: Predict Risk Scores
# =============================================================================

# Clinical ridge — linear predictor
train_df$risk_clinical <- as.numeric(predict(
  ridge_clinical,
  newx = X_train_clin,
  type = "link"
))
test_df$risk_clinical <- as.numeric(predict(
  ridge_clinical,
  newx = X_test_clin,
  type = "link"
))

# CAFs ridge — linear predictor (log scale, monotone with risk)
train_df$risk_CAFs <- as.numeric(predict(
  ridge_CAFs,
  newx = X_train_CAFs,
  type = "link"
))
test_df$risk_CAFs <- as.numeric(predict(
  ridge_CAFs,
  newx = X_test_CAFs,
  type = "link"
))

# Full ridge
train_df$risk_full <- as.numeric(predict(
  ridge_full,
  newx = X_train_full,
  type = "link"
))
test_df$risk_full <- as.numeric(predict(
  ridge_full,
  newx = X_test_full,
  type = "link"
))

# =============================================================================
# Section 8: Time-dependent ROC per Model (Train vs Test, at 3 years)
# =============================================================================
predict_time <- 3

# --- 8a. Full model (Clinical + CAFs) ---
roc_train_full <- survivalROC(
  Stime = train_df$duration,
  status = train_df$status,
  marker = train_df$risk_full,
  predict.time = predict_time,
  method = "KM"
)
roc_test_full <- survivalROC(
  Stime = test_df$duration,
  status = test_df$status,
  marker = test_df$risk_full,
  predict.time = predict_time,
  method = "KM"
)
par(mfrow = c(1, 2))
plot(
  roc_train_full$FP,
  roc_train_full$TP,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = paste0(
    "Train | Clinical + CAFs (Ridge)\n",
    predict_time,
    "-yr OS  AUC = ",
    round(roc_train_full$AUC, 3)
  )
)
abline(a = 0, b = 1, lty = 2, col = "gray")
plot(
  roc_test_full$FP,
  roc_test_full$TP,
  type = "l",
  col = "red",
  lwd = 2,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = paste0(
    "Test | Clinical + CAFs (Ridge)\n",
    predict_time,
    "-yr OS  AUC = ",
    round(roc_test_full$AUC, 3)
  )
)
abline(a = 0, b = 1, lty = 2, col = "gray")
par(mfrow = c(1, 1))

# --- 8b. CAFs only ---
roc_train_CAFs <- survivalROC(
  Stime = train_df$duration,
  status = train_df$status,
  marker = train_df$risk_CAFs,
  predict.time = predict_time,
  method = "KM"
)
roc_test_CAFs <- survivalROC(
  Stime = test_df$duration,
  status = test_df$status,
  marker = test_df$risk_CAFs,
  predict.time = predict_time,
  method = "KM"
)
par(mfrow = c(1, 2))
plot(
  roc_train_CAFs$FP,
  roc_train_CAFs$TP,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = paste0(
    "Train | CAFs only (Ridge)\n",
    predict_time,
    "-yr OS  AUC = ",
    round(roc_train_CAFs$AUC, 3)
  )
)
abline(a = 0, b = 1, lty = 2, col = "gray")
plot(
  roc_test_CAFs$FP,
  roc_test_CAFs$TP,
  type = "l",
  col = "red",
  lwd = 2,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = paste0(
    "Test | CAFs only (Ridge)\n",
    predict_time,
    "-yr OS  AUC = ",
    round(roc_test_CAFs$AUC, 3)
  )
)
abline(a = 0, b = 1, lty = 2, col = "gray")
par(mfrow = c(1, 1))

# --- 8c. Clinical only ---
roc_train_clin <- survivalROC(
  Stime = train_df$duration,
  status = train_df$status,
  marker = train_df$risk_clinical,
  predict.time = predict_time,
  method = "KM"
)
roc_test_clin <- survivalROC(
  Stime = test_df$duration,
  status = test_df$status,
  marker = test_df$risk_clinical,
  predict.time = predict_time,
  method = "KM"
)
par(mfrow = c(1, 2))
plot(
  roc_train_clin$FP,
  roc_train_clin$TP,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = paste0(
    "Train | Clinical only\n",
    predict_time,
    "-yr OS  AUC = ",
    round(roc_train_clin$AUC, 3)
  )
)
abline(a = 0, b = 1, lty = 2, col = "gray")
plot(
  roc_test_clin$FP,
  roc_test_clin$TP,
  type = "l",
  col = "red",
  lwd = 2,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = paste0(
    "Test | Clinical only\n",
    predict_time,
    "-yr OS  AUC = ",
    round(roc_test_clin$AUC, 3)
  )
)
abline(a = 0, b = 1, lty = 2, col = "gray")
par(mfrow = c(1, 1))

# =============================================================================
# Section 9: Compare 3 Models on Test Set (at 5 years)
# =============================================================================
predict_time <- 5

roc_clin <- survivalROC(
  Stime = test_df$duration,
  status = test_df$status,
  marker = test_df$risk_clinical,
  predict.time = predict_time,
  method = "KM"
)
roc_CAFS <- survivalROC(
  Stime = test_df$duration,
  status = test_df$status,
  marker = test_df$risk_CAFs,
  predict.time = predict_time,
  method = "KM"
)
roc_full <- survivalROC(
  Stime = test_df$duration,
  status = test_df$status,
  marker = test_df$risk_full,
  predict.time = predict_time,
  method = "KM"
)

plot(
  roc_clin$FP,
  roc_clin$TP,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = paste0("Test Set: ", predict_time, "-Year OS | Comparison of 3 Models")
)
lines(roc_CAFS$FP, roc_CAFS$TP, col = "green", lwd = 2)
lines(roc_full$FP, roc_full$TP, col = "red", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")
legend(
  "bottomright",
  legend = c(
    paste0("Clinical only   AUC = ", round(roc_clin$AUC, 3)),
    paste0("CAFs only       AUC = ", round(roc_CAFS$AUC, 3)),
    paste0("Clinical + CAFs AUC = ", round(roc_full$AUC, 3))
  ),
  col = c("blue", "green", "red"),
  lwd = 2
)

# =============================================================================
# Section 10: Time-dependent AUC (1-5 years, Train vs Test)
# =============================================================================
time_points <- 1:5

compute_auc_df <- function(df, split_label) {
  do.call(
    rbind,
    lapply(time_points, function(tp) {
      roc_clin <- survivalROC(
        Stime = df$duration,
        status = df$status,
        marker = df$risk_clinical,
        predict.time = tp,
        method = "KM"
      )
      roc_CAFS <- survivalROC(
        Stime = df$duration,
        status = df$status,
        marker = df$risk_CAFs,
        predict.time = tp,
        method = "KM"
      )
      roc_full <- survivalROC(
        Stime = df$duration,
        status = df$status,
        marker = df$risk_full,
        predict.time = tp,
        method = "KM"
      )
      data.frame(
        time = rep(tp, 3),
        model = c("Clinical only", "CAFs only", "Clinical + CAFs"),
        AUC = c(roc_clin$AUC, roc_CAFS$AUC, roc_full$AUC),
        split = split_label
      )
    })
  )
}

auc_df <- rbind(
  compute_auc_df(train_df, "Train"),
  compute_auc_df(test_df, "Test")
)

auc_plot_theme <- list(
  geom_line(linewidth = 1.2),
  geom_point(size = 3),
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60"),
  scale_x_continuous(breaks = time_points),
  scale_y_continuous(limits = c(0.4, 1)),
  scale_color_manual(
    values = c(
      "Clinical only" = "#1F77B4",
      "CAFs only" = "#2CA02C",
      "Clinical + CAFs" = "#D62728"
    )
  ),
  labs(x = "Time (Years)", y = "AUC", color = "Model"),
  theme_bw(),
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )
)

p_train <- ggplot(
  subset(auc_df, split == "Train"),
  aes(x = time, y = AUC, color = model, group = model)
) +
  auc_plot_theme +
  labs(title = "Time-dependent AUC — Train (Ridge Cox)")

p_test <- ggplot(
  subset(auc_df, split == "Test"),
  aes(x = time, y = AUC, color = model, group = model)
) +
  auc_plot_theme +
  labs(title = "Time-dependent AUC — Test (Ridge Cox)")

p_train +
  p_test +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")



