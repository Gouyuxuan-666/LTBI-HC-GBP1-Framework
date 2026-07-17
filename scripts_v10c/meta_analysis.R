###############################################################################
# GBP1 Meta-Analysis: A) Cross-cohort expression + B) Diagnostic accuracy
###############################################################################
setwd("F:/GBP1_pipeline_2hao")

library(limma)
library(metafor)
library(pROC)

# ---- Load merged data ----
load("6/merged_expr.RData")  # expr=merged_combat, type_vec, batch_vec, common_genes

# ================================================================
# A: GBP1 Cross-Cohort Expression Meta-Analysis (Forest Plot)
# ================================================================
cat("\n========== A: GBP1 Expression Meta-Analysis ==========\n")

# Run per-cohort limma to get GBP1 logFC + SE for each dataset
cohorts <- unique(batch_vec)
meta_results <- data.frame(
  cohort = character(),
  n_total = integer(),
  n_tb = integer(),
  n_control = integer(),
  gbp1_logFC = numeric(),
  gbp1_SE = numeric(),
  gbp1_P = numeric(),
  stringsAsFactors = FALSE
)

for (coh in cohorts) {
  idx <- batch_vec == coh
  sub_expr <- merged_combat[, idx]
  sub_type <- type_vec[idx]

  if (sum(sub_type == "Treat") < 3 || sum(sub_type == "Control") < 3) next

  design <- model.matrix(~ sub_type)
  fit <- lmFit(sub_expr, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf)

  if ("GBP1" %in% rownames(tt)) {
    gbp1 <- tt["GBP1", ]
    meta_results <- rbind(meta_results, data.frame(
      cohort = coh,
      n_total = ncol(sub_expr),
      n_tb = sum(sub_type == "Treat"),
      n_control = sum(sub_type == "Control"),
      gbp1_logFC = gbp1$logFC,
      gbp1_SE = abs(gbp1$logFC) / abs(gbp1$t),
      gbp1_P = gbp1$P.Value,
      stringsAsFactors = FALSE
    ))
  }
}

# Also add external validation cohorts (if available)
# GSE19444 and GSE39940 would be loaded separately

cat(sprintf("Cohorts with GBP1 data: %d\n", nrow(meta_results)))
print(meta_results[, c("cohort", "n_total", "gbp1_logFC", "gbp1_P")])

# Random-effects meta-analysis
if (nrow(meta_results) >= 4) {
  ma <- rma(yi = gbp1_logFC, sei = gbp1_SE, data = meta_results, method = "REML")

  # Forest plot
  pdf("meta_analysis/Fig_GBP1_Forest.pdf", width = 10, height = 6)
  forest(ma, slab = meta_results$cohort,
         xlab = "GBP1 log2 Fold Change (TB vs Control)",
         header = "Cohort",
         mlab = "RE Model",
         ilab = cbind(meta_results$n_tb, meta_results$n_control),
         ilab.xpos = c(min(ma$yi) - 1.5, min(ma$yi) - 0.8),
         cex = 0.9)
  text(min(ma$yi) - 1.15, nrow(meta_results) + 2, "TB", cex = 0.7, font = 2)
  text(min(ma$yi) - 0.55, nrow(meta_results) + 2, "Control", cex = 0.7, font = 2)
  dev.off()

  cat(sprintf("\nMeta-analysis results:\n"))
  cat(sprintf("  Pooled logFC: %.3f (95%% CI: %.3f-%.3f)\n", ma$b, ma$ci.lb, ma$ci.ub))
  cat(sprintf("  I²: %.1f%%, Q P = %.4f\n", ma$I2, ma$QEp))
  cat(sprintf("  Tau²: %.4f\n", ma$tau2))

  # Save results
  sink("meta_analysis/meta_expression_summary.txt")
  cat(sprintf("GBP1 Cross-Cohort Expression Meta-Analysis\n"))
  cat(sprintf("==========================================\n"))
  cat(sprintf("Cohorts: %d\n", nrow(meta_results)))
  cat(sprintf("Pooled logFC: %.3f (95%% CI: %.3f-%.3f)\n", ma$b, ma$ci.lb, ma$ci.ub))
  cat(sprintf("I²: %.1f%%\n", ma$I2))
  cat(sprintf("Q statistic: %.2f, df=%d, P=%.4f\n", ma$QE, ma$k-1, ma$QEp))
  cat(sprintf("Tau²: %.4f\n", ma$tau2))
  cat(sprintf("P for pooled effect: %.2e\n", ma$pval))
  cat(sprintf("\nPer-cohort:\n"))
  print(meta_results)
  sink()
}

# ================================================================
# B: Diagnostic Accuracy Meta-Analysis
# ================================================================
cat("\n========== B: Diagnostic Accuracy Meta-Analysis ==========\n")

# For each cohort, train LDA+RF and get AUC
dir.create("meta_analysis", showWarnings = FALSE)

library(caret)
library(randomForest)

# Load the DEG-turquoise intersection genes for ML
load("8/wgcna_results.RData")
deg_tab <- read.table("7/DEG_full.txt", header = TRUE, sep = "\t")
degs <- read.table("7/DEG_filtered.txt", header = TRUE, sep = "\t")

ml_genes_all <- intersect(rownames(degs), turquoise_genes)
fc_all <- abs(deg_tab[ml_genes_all, "logFC"]); names(fc_all) <- ml_genes_all
ml_genes_top200 <- names(sort(fc_all, decreasing = TRUE))[1:min(200, length(fc_all))]
ml_expr <- t(merged_combat[ml_genes_top200, ])

acc_results <- data.frame(
  cohort = character(),
  auc = numeric(),
  sens = numeric(),
  spec = numeric(),
  n_samples = integer(),
  stringsAsFactors = FALSE
)

set.seed(20240715)
for (coh in cohorts) {
  idx <- batch_vec == coh
  if (sum(idx) < 10) next

  coh_expr <- as.data.frame(ml_expr[idx, ])
  coh_labels <- factor(ifelse(type_vec[idx] == "Treat", "TB", "Control"))

  if (length(unique(coh_labels)) < 2) next
  if (min(table(coh_labels)) < 3) next

  # Leave-one-cohort-out training: train on all OTHER cohorts
  train_idx <- batch_vec != coh
  train_x <- as.data.frame(ml_expr[train_idx, ])
  train_y <- factor(ifelse(type_vec[train_idx] == "Treat", "TB", "Control"))

  # LDA feature selection on training
  lda_subset <- train_x[, 1:min(30, ncol(train_x)), drop = FALSE]
  lda_fit <- tryCatch(MASS::lda(x = lda_subset, grouping = train_y), error = function(e) NULL)
  if (is.null(lda_fit)) next
  lda_scores <- abs(lda_fit$scaling[, 1])
  sel_genes <- names(sort(lda_scores, decreasing = TRUE))[1:min(20, length(lda_scores))]

  # RF on selected genes
  ctrl <- trainControl(method = "cv", number = 3, classProbs = TRUE, summaryFunction = twoClassSummary)
  rf_model <- tryCatch(
    caret::train(x = train_x[, sel_genes, drop = FALSE], y = train_y,
                 method = "rf", trControl = ctrl, metric = "ROC", tuneLength = 2),
    error = function(e) NULL
  )
  if (is.null(rf_model)) next

  # Predict on held-out cohort
  test_pred <- predict(rf_model, newdata = coh_expr[, sel_genes, drop = FALSE], type = "prob")[, "TB"]
  coh_roc <- roc(coh_labels, test_pred, quiet = TRUE)

  youden <- coords(coh_roc, "best", ret = c("sensitivity", "specificity"))

  acc_results <- rbind(acc_results, data.frame(
    cohort = coh,
    auc = as.numeric(coh_roc$auc),
    sens = as.numeric(youden[1]),
    spec = as.numeric(youden[2]),
    n_samples = sum(idx),
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  %s: AUC=%.3f, Sens=%.3f, Spec=%.3f (n=%d)\n",
              coh, as.numeric(coh_roc$auc), youden[1], youden[2], sum(idx)))
}

# AUC forest plot
if (nrow(acc_results) >= 3) {
  acc_results$auc_SE <- sqrt(acc_results$auc * (1 - acc_results$auc) / acc_results$n_samples)
  acc_results$cohort_label <- paste0(acc_results$cohort, " (n=", acc_results$n_samples, ")")

  ma_auc <- rma(yi = auc, sei = auc_SE, data = acc_results, method = "REML")

  pdf("meta_analysis/Fig_AUC_Forest.pdf", width = 9, height = 5)
  forest(ma_auc, slab = acc_results$cohort_label,
         xlab = "AUC (leave-one-cohort-out validation)",
         header = "Cohort", mlab = "RE Model",
         refline = 0.5, cex = 0.85)
  dev.off()

  cat(sprintf("\nAUC Meta-analysis:\n"))
  cat(sprintf("  Pooled AUC: %.3f (95%% CI: %.3f-%.3f)\n", ma_auc$b, ma_auc$ci.lb, ma_auc$ci.ub))
  cat(sprintf("  I²: %.1f%%\n", ma_auc$I2))

  # SROC-like scatter plot
  pdf("meta_analysis/Fig_SROC.pdf", width = 7, height = 6)
  plot(acc_results$spec, acc_results$sens,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Specificity", ylab = "Sensitivity",
       main = "GBP1 Diagnostic Performance Across Cohorts",
       pch = 16, cex = acc_results$n_samples / 30,
       col = rgb(0.2, 0.4, 0.8, 0.7))
  text(acc_results$spec, acc_results$sens + 0.05, acc_results$cohort, cex = 0.7)
  abline(h = 0.8, lty = 2, col = "gray")
  abline(v = 0.8, lty = 2, col = "gray")
  legend("bottomleft", c("WHO TPP triage: Sens>=90%, Spec>=70%"),
         lty = 2, col = "gray", cex = 0.7, bty = "n")
  dev.off()

  write.csv(acc_results, "meta_analysis/auc_per_cohort.csv", row.names = FALSE)
}

cat("\nDone! Outputs in meta_analysis/\n")
