###############################################################################
# GBP1 vs GBP5 Single-Gene Diagnostic Performance Comparison
# DeLong test + NRI + IDI
# Uses merged discovery matrix (post-ComBat) + external validation
###############################################################################

library(pROC)
library(PredictABEL)  # for NRI/IDI
library(ggplot2)

# ---- Load merged expression matrix ----
# Assumes: expr_mat = merged ComBat-corrected matrix (genes x samples)
#          meta = sample metadata with TB/Control labels
# Adjust paths to your 2hao machine setup

setwd("F:/GBP1_pipeline_2hao")

# Load the merged data (post step 6-7)
load("6/merged_expr.RData")  # adjust filename
# Expected: expr (matrix), meta (data.frame with 'condition' column: Treat/Control)

# ---- 1. Extract GBP1 and GBP5 expression ----
gbp1_expr <- as.numeric(expr["GBP1", ])
gbp5_expr <- as.numeric(expr["GBP5", ])
tb_status <- ifelse(meta$condition == "Treat", 1, 0)

# ---- 2. Single-gene logistic regression models ----
gbp1_glm <- glm(tb_status ~ gbp1_expr, family = binomial)
gbp5_glm <- glm(tb_status ~ gbp5_expr, family = binomial)

gbp1_pred <- predict(gbp1_glm, type = "response")
gbp5_pred <- predict(gbp5_glm, type = "response")

# ---- 3. ROC curves ----
roc_gbp1 <- roc(tb_status, gbp1_pred, ci = TRUE)
roc_gbp5 <- roc(tb_status, gbp5_pred, ci = TRUE)

# ---- 4. DeLong test: H0: AUC_GBP1 = AUC_GBP5 ----
delong_test <- roc.test(roc_gbp1, roc_gbp5, method = "delong")

# ---- 5. NRI and IDI ----
# NRI: Net Reclassification Improvement
# IDI: Integrated Discrimination Improvement
nri_result <- reclassification(data = data.frame(tb_status, gbp1_pred, gbp5_pred),
                                cOutcome = 1,
                                predrisk1 = gbp5_pred,  # old model
                                predrisk2 = gbp1_pred,  # new model (GBP1 added to or vs GBP5)
                                cutoff = c(0, 0.3, 0.7, 1))

# Alternative: use PredictABEL for IDI
idi_result <- improveProb(x1 = gbp5_pred,
                           x2 = gbp1_pred,
                           y = tb_status)

# ---- 6. Sensitivity/Specificity at Youden threshold ----
youden_gbp1 <- coords(roc_gbp1, "best", ret = c("threshold", "sensitivity", "specificity"))
youden_gbp5 <- coords(roc_gbp5, "best", ret = c("threshold", "sensitivity", "specificity"))

# ---- 7. Also run on external validation sets ----
# For GSE19444 and GSE39940 separately (if data loaded)

# ---- 8. Compile results ----
results <- list(
    auc_gbp1 = c(AUC = as.numeric(roc_gbp1$auc), CI_low = as.numeric(ci(roc_gbp1)[1]), CI_high = as.numeric(ci(roc_gbp1)[3])),
    auc_gbp5 = c(AUC = as.numeric(roc_gbp5$auc), CI_low = as.numeric(ci(roc_gbp5)[1]), CI_high = as.numeric(ci(roc_gbp5)[3])),
    delta_auc = as.numeric(roc_gbp1$auc - roc_gbp5$auc),
    delong_p = delong_test$p.value,
    delong_z = delong_test$statistic,
    nri = nri_result,
    idi = idi_result,
    sens_gbp1 = youden_gbp1["sensitivity"],
    spec_gbp1 = youden_gbp1["specificity"],
    sens_gbp5 = youden_gbp5["sensitivity"],
    spec_gbp5 = youden_gbp5["specificity"]
)

# ---- 9. Save ----
saveRDS(results, "17/gbp1_vs_gbp5_compare.rds")

# Print summary
cat(sprintf("\n========== GBP1 vs GBP5 Single-Gene Comparison ==========\n"))
cat(sprintf("GBP1 AUC: %.3f (95%% CI: %.3f-%.3f)\n",
            results$auc_gbp1[1], results$auc_gbp1[2], results$auc_gbp1[3]))
cat(sprintf("GBP5 AUC: %.3f (95%% CI: %.3f-%.3f)\n",
            results$auc_gbp5[1], results$auc_gbp5[2], results$auc_gbp5[3]))
cat(sprintf("Delta AUC (GBP1 - GBP5): %.4f\n", results$delta_auc))
cat(sprintf("DeLong Z = %.3f, P = %.4f\n", results$delong_z, results$delong_p))
cat(sprintf("GBP1 Sens=%.3f, Spec=%.3f (Youden)\n", results$sens_gbp1, results$spec_gbp1))
cat(sprintf("GBP5 Sens=%.3f, Spec=%.3f (Youden)\n", results$sens_gbp5, results$spec_gbp5))
cat(sprintf("NRI: %.4f, IDI: %.4f\n", results$nri, results$idi))
cat(sprintf("=========================================================\n"))
