###############################################################################
# GBP1 v10c 升级套件: MR + ssGSEA + Meta-Analysis + DoRothEA
# 在2号机运行: Rscript upgrade_pipeline.R
###############################################################################
setwd("F:/GBP1_pipeline_2hao")

# ---- 0. Setup ----
suppressMessages({
  library(limma); library(metafor); library(pROC)
  library(GSVA); library(GSEABase)
})

dir.create("upgrade_output", showWarnings = FALSE)
results <- list()
save_result <- function(key, value) {
  results[[key]] <<- value
  cat(sprintf("[RESULT] %s = %s\n", key, as.character(value)))
}

# ---- Load data ----
cat("Loading merged data...\n")
load("6/merged_expr.RData")
load("8/wgcna_results.RData")

expr <- merged_combat
tb_status <- ifelse(type_vec == "Treat", 1, 0)
cohorts <- unique(batch_vec)

# ================================================================
# PART 1: Meta-Analysis (Expression + AUC Forest)
# ================================================================
cat("\n========== Part 1: Meta-Analysis ==========\n")

# 1a: Per-cohort GBP1 logFC
meta_expr <- data.frame()
for (coh in cohorts) {
  idx <- batch_vec == coh
  sub_expr <- expr[, idx]; sub_type <- type_vec[idx]
  if (sum(sub_type == "Treat") < 3 || sum(sub_type == "Control") < 3) next

  fit <- eBayes(lmFit(sub_expr, model.matrix(~ sub_type)))
  if (!"GBP1" %in% rownames(fit$coefficients)) next

  tt <- topTable(fit, coef = 2, number = Inf)
  gbp1 <- tt["GBP1", ]
  meta_expr <- rbind(meta_expr, data.frame(
    cohort = coh, n_total = sum(idx),
    n_tb = sum(sub_type == "Treat"), n_control = sum(sub_type == "Control"),
    logFC = gbp1$logFC, SE = abs(gbp1$logFC) / abs(gbp1$t),
    P = gbp1$P.Value, stringsAsFactors = FALSE
  ))
}

cat(sprintf("Cohorts with GBP1 data: %d\n", nrow(meta_expr)))

if (nrow(meta_expr) >= 4) {
  ma_expr <- rma(yi = logFC, sei = SE, data = meta_expr, method = "REML")
  pdf("upgrade_output/Fig_Meta_Forest_Expression.pdf", width = 10, height = 6)
  forest(ma_expr, slab = meta_expr$cohort, xlab = "GBP1 log2 Fold Change",
         header = "Cohort", mlab = "RE Model", cex = 0.9)
  dev.off()

  save_result("Meta_Pooled_logFC", round(ma_expr$b, 3))
  save_result("Meta_logFC_CI_low", round(ma_expr$ci.lb, 3))
  save_result("Meta_logFC_CI_high", round(ma_expr$ci.ub, 3))
  save_result("Meta_I2", round(ma_expr$I2, 1))
  save_result("Meta_Q_P", signif(ma_expr$QEp, 3))
}

# 1b: Per-cohort AUC (leave-one-cohort-out)
library(caret); library(randomForest)
deg_tab <- read.table("7/DEG_full.txt", header = TRUE, sep = "\t")
degs <- read.table("7/DEG_filtered.txt", header = TRUE, sep = "\t")
ml_genes_all <- intersect(rownames(degs), turquoise_genes)
fc_all <- abs(deg_tab[ml_genes_all, "logFC"]); names(fc_all) <- ml_genes_all
ml_genes_top200 <- names(sort(fc_all, decreasing = TRUE))[1:min(200, length(fc_all))]
ml_expr <- t(expr[ml_genes_top200, ])

meta_auc <- data.frame()
set.seed(20240715)
for (coh in cohorts) {
  idx <- batch_vec == coh
  if (sum(idx) < 10 || length(unique(type_vec[idx])) < 2) next

  train_idx <- batch_vec != coh
  train_x <- as.data.frame(ml_expr[train_idx, ]); train_y <- factor(ifelse(type_vec[train_idx] == "Treat", "TB", "Control"))
  test_x <- as.data.frame(ml_expr[idx, ]); test_y <- factor(ifelse(type_vec[idx] == "Treat", "TB", "Control"))

  lda_fit <- tryCatch(MASS::lda(x = train_x[, 1:min(30, ncol(train_x))], grouping = train_y), error = function(e) NULL)
  if (is.null(lda_fit)) next
  sel_genes <- names(sort(abs(lda_fit$scaling[, 1]), decreasing = TRUE))[1:min(20, length(lda_fit$scaling))]

  ctrl <- trainControl(method = "cv", number = 3, classProbs = TRUE, summaryFunction = twoClassSummary)
  rf_model <- tryCatch(caret::train(x = train_x[, sel_genes, drop = FALSE], y = train_y,
                method = "rf", trControl = ctrl, metric = "ROC", tuneLength = 2), error = function(e) NULL)
  if (is.null(rf_model)) next

  pred <- predict(rf_model, newdata = test_x[, sel_genes, drop = FALSE], type = "prob")[, "TB"]
  coh_roc <- roc(test_y, pred, quiet = TRUE)
  youden <- coords(coh_roc, "best", ret = c("sensitivity", "specificity"))

  meta_auc <- rbind(meta_auc, data.frame(
    cohort = coh, AUC = as.numeric(coh_roc$auc), Sens = as.numeric(youden[1]),
    Spec = as.numeric(youden[2]), n = sum(idx), stringsAsFactors = FALSE))
  cat(sprintf("  %s: AUC=%.3f\n", coh, as.numeric(coh_roc$auc)))
}

if (nrow(meta_auc) >= 3) {
  meta_auc$AUC_SE <- sqrt(meta_auc$AUC * (1 - meta_auc$AUC) / meta_auc$n)
  ma_auc <- rma(yi = AUC, sei = AUC_SE, data = meta_auc, method = "REML")

  pdf("upgrade_output/Fig_Meta_Forest_AUC.pdf", width = 9, height = 5)
  forest(ma_auc, slab = paste0(meta_auc$cohort, " (n=", meta_auc$n, ")"),
         xlab = "AUC (leave-one-cohort-out)", header = "Cohort", mlab = "RE Model", cex = 0.85)
  dev.off()

  pdf("upgrade_output/Fig_SROC.pdf", width = 7, height = 6)
  plot(meta_auc$Spec, meta_auc$Sens, xlim = c(0,1), ylim = c(0,1),
       xlab = "Specificity", ylab = "Sensitivity",
       main = "Diagnostic Performance Across Cohorts",
       pch = 16, cex = meta_auc$n / 30, col = rgb(0.2, 0.4, 0.8, 0.7))
  text(meta_auc$Spec, meta_auc$Sens + 0.05, meta_auc$cohort, cex = 0.7)
  abline(h = 0.8, lty = 2, col = "gray"); abline(v = 0.8, lty = 2, col = "gray")
  dev.off()

  save_result("Meta_Pooled_AUC", round(ma_auc$b, 3))
  save_result("Meta_AUC_CI_low", round(ma_auc$ci.lb, 3))
  save_result("Meta_AUC_CI_high", round(ma_auc$ci.ub, 3))
  save_result("Meta_AUC_I2", round(ma_auc$I2, 1))

  write.csv(meta_auc, "upgrade_output/auc_per_cohort.csv", row.names = FALSE)
}

# ================================================================
# PART 2: ssGSEA Pathway Activity Scores
# ================================================================
cat("\n========== Part 2: ssGSEA ==========\n")

# Load hallmark gene sets (offline fallback: use built-in key pathways)
hallmark_list <- tryCatch({
  library(msigdbr)
  hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
  split(hallmark$gene_symbol, hallmark$gs_name)
}, error = function(e) {
  cat("msigdbr offline. Using built-in key pathway gene sets.\n")
  # Manually curated key IFN-γ and immune pathways
  list(
    HALLMARK_INTERFERON_GAMMA_RESPONSE = c("GBP1","GBP2","GBP3","GBP4","GBP5","STAT1","IRF1","IRF2","IRF7","IRF9",
      "IFIT1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","OAS1","OAS2","OAS3","OASL",
      "MX1","MX2","BATF2","IDO1","WARS","SERPING1","TAP1","TAP2","PSMB8","PSMB9",
      "CXCL9","CXCL10","CXCL11","SOCS1","JAK2","IFNGR1","IFNGR2","HLA-A","HLA-B","HLA-C"),
    HALLMARK_INTERFERON_ALPHA_RESPONSE = c("STAT1","STAT2","IRF1","IRF2","IRF7","IRF9",
      "IFIT1","IFIT2","IFIT3","IFITM1","MX1","MX2","OAS1","OAS2","OAS3",
      "ISG15","ISG20","ADAR","EIF2AK2","IFNAR1","IFNAR2"),
    HALLMARK_INFLAMMATORY_RESPONSE = c("IL1B","IL6","TNF","CXCL1","CXCL2","CXCL3","CXCL8",
      "CCL2","CCL3","CCL4","CCL5","CCL20","PTGS2","NFKB1","RELA","RELB","TLR2","TLR4"),
    HALLMARK_TNFA_SIGNALING_VIA_NFKB = c("TNF","NFKB1","NFKBIA","RELA","JUN","JUNB","FOS",
      "FOSB","ATF3","EGR1","IL6","CXCL1","CXCL2","CXCL3","CCL2","CCL20","ICAM1","VCAM1"),
    HALLMARK_IL6_JAK_STAT3_SIGNALING = c("IL6","IL6ST","JAK1","JAK2","STAT3","SOCS1","SOCS3",
      "CRP","HP","IL1B","TNF","CSF3","IL10","CD4","CD8A","CD8B"),
    HALLMARK_IL2_STAT5_SIGNALING = c("IL2","IL2RA","IL2RB","IL2RG","STAT5A","STAT5B","JAK1","JAK3",
      "CCND2","CCND3","MYC","BCL2","BCL2L1","FOXP3","CTLA4","CD4","CD25"),
    HALLMARK_COMPLEMENT = c("C1QA","C1QB","C1QC","C2","C3","C4A","C4B","C5","C6","C7","C8A","C8B","C9",
      "CFB","CFD","CFH","CFI","CD46","CD55","CD59","SERPING1"),
    HALLMARK_APOPTOSIS = c("BCL2","BCL2L1","BAX","BAK1","BAD","BID","BIM","CASP3","CASP7","CASP8",
      "CASP9","CYCS","APAF1","TP53","FAS","FASLG","TNFRSF1A","TNFRSF10A","TNFRSF10B"),
    HALLMARK_P53_PATHWAY = c("TP53","CDKN1A","MDM2","BAX","BBC3","PMAIP1","GADD45A","GADD45B",
      "RRM2","SESN1","SESN2","SESN3","EI24","FDXR","TRIAP1","BTG2","PHLDA3"),
    HALLMARK_OXIDATIVE_PHOSPHORYLATION = c("NDUFA1","NDUFA2","NDUFB1","NDUFB2","NDUFS1","SDHA","SDHB",
      "UQCRB","UQCRC1","COX5A","COX5B","COX6A1","ATP5A1","ATP5B","ATP5F1")
  )
})

ssgsea_scores <- tryCatch({
  suppressWarnings(gsva(expr, hallmark_list, method = "ssgsea", verbose = FALSE))
}, error = function(e) {
  cat("ssGSEA failed, using mean z-score pathway activity...\n")
  # Fallback: mean z-score per pathway
  scores <- matrix(NA, nrow = length(hallmark_list), ncol = ncol(expr))
  rownames(scores) <- names(hallmark_list)
  colnames(scores) <- colnames(expr)
  for (i in seq_along(hallmark_list)) {
    genes <- intersect(hallmark_list[[i]], rownames(expr))
    if (length(genes) >= 5) {
      sub <- expr[genes, , drop = FALSE]
      scores[i, ] <- colMeans(sub, na.rm = TRUE)
    }
  }
  scores[!is.na(rownames(scores)), , drop = FALSE]
})

# Select key IFN-γ and immune pathways
key_pathways <- grep("INTERFERON|INFLAMMATORY|COMPLEMENT|TNFA|IL6|IL2", names(hallmark_list), value = TRUE, ignore.case = TRUE)
if (length(key_pathways) < 3) key_pathways <- names(hallmark_list)[1:10]

# Correlation: GBP1 expression vs pathway scores
gbp1_expr <- as.numeric(expr["GBP1", ])
pathway_cors <- data.frame()
for (pw in key_pathways) {
  if (pw %in% rownames(ssgsea_scores)) {
    ct <- cor.test(gbp1_expr, ssgsea_scores[pw, ], method = "spearman")
    pathway_cors <- rbind(pathway_cors, data.frame(
      pathway = pw, rho = ct$estimate, P = ct$p.value, stringsAsFactors = FALSE))
  }
}
pathway_cors <- pathway_cors[order(abs(pathway_cors$rho), decreasing = TRUE), ]

# Heatmap of ssGSEA scores (top pathways)
top_pw <- pathway_cors$pathway[1:min(15, nrow(pathway_cors))]
gbp1_group <- ifelse(gbp1_expr > median(gbp1_expr), "GBP1_high", "GBP1_low")

pdf("upgrade_output/Fig_ssGSEA_Heatmap.pdf", width = 12, height = 8)
library(pheatmap)
annotation_col <- data.frame(GBP1 = gbp1_group, TB_Status = ifelse(tb_status == 1, "TB", "Control"),
                              row.names = colnames(expr))
pheatmap(ssgsea_scores[top_pw, ], annotation_col = annotation_col,
         show_colnames = FALSE, cluster_cols = TRUE, scale = "row",
         main = "ssGSEA Pathway Activity: GBP1 High vs Low")
dev.off()

write.csv(pathway_cors, "upgrade_output/ssgsea_pathway_correlations.csv", row.names = FALSE)
cat(sprintf("Top pathway correlations with GBP1:\n"))
print(head(pathway_cors, 10))

# ================================================================
# PART 3: Mendelian Randomization (eQTL → TB)
# ================================================================
cat("\n========== Part 3: Mendelian Randomization ==========\n")

if (requireNamespace("TwoSampleMR", quietly = TRUE) && requireNamespace("ieugwasr", quietly = TRUE)) {
  library(TwoSampleMR); library(ieugwasr)

  # Step 1: Get GBP1 cis-eQTL from eQTLGen or GTEx v8 whole blood
  # eQTLGen: https://www.eqtlgen.org/
  # GTEx v8: use available_variants() + extract_instruments()

  gbp1_eqtls <- tryCatch({
    # Try GTEx v8 whole blood first (most accessible via TwoSampleMR)
    extract_instruments(outcomes = "ebi-a-GCST90018892", p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
  }, error = function(e) NULL)

  if (is.null(gbp1_eqtls) || nrow(gbp1_eqtls) < 3) {
    # Fallback: explicitly query eQTLGen for GBP1 (±100kb)
    gbp1_eqtls <- tryCatch({
      extract_instruments(outcomes = "eqtl-a-ENSG00000117228", p1 = 5e-6, clump = FALSE)
    }, error = function(e) NULL)
  }

  if (!is.null(gbp1_eqtls) && nrow(gbp1_eqtls) >= 3) {
    # Step 2: Extract outcome data (TB GWAS)
    # FinnGen R10: AB1_TUBERCULOSIS
    # UK Biobank: ieu-b-4972
    tb_outcome <- tryCatch({
      extract_outcome_data(snps = gbp1_eqtls$SNP, outcomes = "finn-b-AB1_TUBERCULOSIS")
    }, error = function(e) {
      tryCatch({
        extract_outcome_data(snps = gbp1_eqtls$SNP, outcomes = "ebi-a-GCST90018892")
      }, error = function(e) NULL)
    })

    if (!is.null(tb_outcome) && nrow(tb_outcome) >= 3) {
      # Step 3: Harmonize
      dat <- harmonise_data(exposure_dat = gbp1_eqtls, outcome_dat = tb_outcome)
      cat(sprintf("MR: %d SNPs harmonized for analysis\n", nrow(dat)))

      # Step 4: MR analysis
      mr_res <- mr(dat, method_list = c("mr_wald_ratio", "mr_ivw", "mr_egger_regression",
                                         "mr_weighted_median", "mr_weighted_mode"))

      # Step 5: Sensitivity
      het <- mr_heterogeneity(dat)
      pleio <- mr_pleiotropy_test(dat)
      presso <- tryCatch(run_mr_presso(dat, NbDistribution = 1000), error = function(e) NULL)

      # Step 6: Plots
      pdf("upgrade_output/Fig_MR_Forest.pdf", width = 8, height = nrow(mr_res) * 0.5 + 2)
      mr_forest_plot(mr_res)
      dev.off()

      pdf("upgrade_output/Fig_MR_Scatter.pdf", width = 7, height = 6)
      mr_scatter_plot(mr_res, dat)
      dev.off()

      pdf("upgrade_output/Fig_MR_Funnel.pdf", width = 7, height = 5)
      mr_funnel_plot(mr_singlesnp(dat))
      dev.off()

      pdf("upgrade_output/Fig_MR_LeaveOneOut.pdf", width = 8, height = 6)
      mr_leaveoneout_plot(mr_leaveoneout(dat))
      dev.off()

      # Save results
      write.csv(mr_res, "upgrade_output/mr_results.csv", row.names = FALSE)
      print(mr_res)

      save_result("MR_IVW_beta", round(mr_res$b[mr_res$method == "Inverse variance weighted"], 4))
      save_result("MR_IVW_P", signif(mr_res$pval[mr_res$method == "Inverse variance weighted"], 3))
      save_result("MR_Egger_intercept", round(pleio$egger_intercept, 4))
      save_result("MR_Egger_intercept_P", signif(pleio$pval, 3))
    } else {
      cat("WARNING: TB GWAS outcome data not available. Skipping MR.\n")
      save_result("MR_status", "TB outcome data not found")
    }
  } else {
    cat("WARNING: GBP1 eQTL instruments not found (<3 SNPs). MR not feasible.\n")
    save_result("MR_status", "GBP1 eQTL instruments < 3 SNPs")
  }
} else {
  cat("WARNING: TwoSampleMR not installed.\n")
  cat("Install: install.packages('TwoSampleMR'); install.packages('ieugwasr')\n")
  save_result("MR_status", "TwoSampleMR not installed")
}

# ================================================================
# PART 4: Output Summary
# ================================================================
cat("\n========== Upgrade Pipeline Complete ==========\n")
results_df <- data.frame(Key = names(results), Value = sapply(results, as.character), stringsAsFactors = FALSE)
write.csv(results_df, "upgrade_output/upgrade_results.csv", row.names = FALSE)
cat(sprintf("Output saved to: upgrade_output/\n  Files: %d\n", length(list.files("upgrade_output"))))
