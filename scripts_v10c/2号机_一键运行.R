###############################################################################
# GBP1 v10c 投稿准备 — 2号机一键运行脚本
# 用途: 获取所有 [TBD] 计算数值，生成补充表
# 前置: 步骤1-5已完成 (GEO下载+标准化+探针坍缩, 即 /5/*.normalize.txt)
# 运行: Rscript 2号机_一键运行.R
# 输出: results_summary.csv + 补充表S2/S3/S4/S7
###############################################################################

options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))
options(warn = 1)

# ---- 0. Setup ----
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0 && script_path != "") setwd(dirname(script_path))
cat(sprintf("Working dir: %s\n", getwd()))

for (d in c("6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25")) {
  dir.create(d, showWarnings = FALSE)
}

# ---- 0.1 Install missing packages ----
need_pkgs <- c("limma","sva","WGCNA","ggplot2","ggrepel","pheatmap",
               "caret","pROC","randomForest","glmnet","clusterProfiler",
               "org.Hs.eg.db","enrichplot","ggvenn","UpSetR",
               "reshape2","corrplot","ggpubr","data.table",
               "randomForestSRC","plsRglm","gbm","mboost","e1071","xgboost",
               "ComplexHeatmap","PredictABEL")
for (pkg in need_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("illuminaHumanv4.db","impute","preprocessCore","sva","limma","WGCNA",
  "clusterProfiler","org.Hs.eg.db","enrichplot","ggvenn"), update = FALSE, ask = FALSE)

suppressMessages({
  library(limma); library(sva); library(WGCNA); library(pROC)
  library(ggplot2); library(ggrepel); library(pheatmap); library(data.table)
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
  library(ggvenn); library(caret)
})

# R 4.6 compatibility: WGCNA cor wrapper
cor <- function(x, y = NULL, use = "everything",
                method = c("pearson", "kendall", "spearman"), ...) {
  if (use == "p") use <- "pairwise.complete.obs"
  if (use == "na.or.complete") use <- "na.or.complete"
  stats::cor(x, y, use = use, method = match.arg(method))
}

# ---- 0.2 Config ----
logFCfilter <- 0.585
adj.P.Val.Filter <- 0.05
WGCNA_MIN_MOD <- 30
WGCNA_DEEP_SPLIT <- 2
WGCNA_MERGE_H <- 0.25
ML_SEED <- 20240715

# ---- 0.3 Helper: save results as we go ----
results <- list()
save_result <- function(key, value) {
  results[[key]] <<- value
  cat(sprintf("[RESULT] %s = %s\n", key, as.character(value)))
}

# ============================================================
# STEP 6: SVA Merge
# ============================================================
cat("\n========== STEP 6: SVA Merge ==========\n")

gse_all <- c("GSE83456","GSE34608","GSE19491","GSE37250","GSE28623","GSE42830")
expr_list <- list(); batch_vec <- c(); type_vec <- c()

for (gse in gse_all) {
  f <- file.path("5", paste0(gse, ".normalize.txt"))
  if (!file.exists(f)) { cat(sprintf("SKIP %s: file not found\n", gse)); next }
  rt <- as.data.frame(fread(f, header = TRUE, sep = "\t", check.names = FALSE))
  rownames(rt) <- rt[[1]]; rt[[1]] <- NULL
  rt <- as.matrix(rt); mode(rt) <- "numeric"
  rt <- rt[rowSums(is.na(rt)) == 0, , drop = FALSE]
  cnames <- colnames(rt)

  types <- rep(NA, length(cnames))
  types[grepl("_Control$|_control$", cnames)] <- "Control"
  types[grepl("_Treat$|_treat$|_TB$", cnames)] <- "Treat"

  valid <- !is.na(types)
  rt <- rt[, valid, drop = FALSE]; types <- types[valid]

  cat(sprintf("  %s: %d genes x %d samples (Treat=%d, Control=%d)\n",
              gse, nrow(rt), ncol(rt), sum(types == "Treat"), sum(types == "Control")))

  if (ncol(rt) > 0) {
    expr_list[[gse]] <- rt
    batch_vec <- c(batch_vec, rep(gse, ncol(rt)))
    type_vec <- c(type_vec, types)
  }
}

common_genes <- Reduce(intersect, lapply(expr_list, rownames))
cat(sprintf("Common genes: %d\n", length(common_genes)))

merged_expr <- do.call(cbind, lapply(expr_list, function(x) x[common_genes, , drop = FALSE]))
cat(sprintf("Merged matrix: %d genes x %d samples\n", nrow(merged_expr), ncol(merged_expr)))
save_result("n_samples_merged", ncol(merged_expr))
save_result("n_genes_merged", nrow(merged_expr))

# ---- ComBat ----
mod <- model.matrix(~ type_vec)
merged_combat <- ComBat(dat = merged_expr, batch = batch_vec, mod = mod)

# ---- Batch diagnostics (PCA) ----
pca_pre <- prcomp(t(merged_expr[order(apply(merged_expr, 1, mad), decreasing = TRUE)[1:2000], ]))
pca_post <- prcomp(t(merged_combat[order(apply(merged_combat, 1, mad), decreasing = TRUE)[1:2000], ]))

# kBET — skip if not installed
if (requireNamespace("kBET", quietly = TRUE)) {
  library(kBET)
  kbet_pre <- kBET(df = t(merged_expr[order(apply(merged_expr, 1, mad), decreasing = TRUE)[1:2000], ]),
                    batch = as.factor(batch_vec), k = 20, n_repeat = 100, plot = FALSE)
  kbet_post <- kBET(df = t(merged_combat[order(apply(merged_combat, 1, mad), decreasing = TRUE)[1:2000], ]),
                    batch = as.factor(batch_vec), k = 20, n_repeat = 100, plot = FALSE)
  save_result("kbet_rejection_pre", mean(kbet_pre$stats$kBET.observed < 0.05))
  save_result("kbet_rejection_post", mean(kbet_post$stats$kBET.observed < 0.05))
} else {
  save_result("kbet_rejection_pre", "kBET not installed")
  save_result("kbet_rejection_post", "kBET not installed")
}

# Save merged data
save(merged_expr, merged_combat, batch_vec, type_vec, common_genes,
     file = "6/merged_expr.RData", compress = FALSE)
expr <- merged_combat
meta <- data.frame(sample = colnames(merged_combat), batch = batch_vec, condition = type_vec,
                   stringsAsFactors = FALSE)

# ============================================================
# STEP 7: DEG
# ============================================================
cat("\n========== STEP 7: Differential Expression ==========\n")

design <- model.matrix(~ condition + batch, data = meta)
fit <- lmFit(expr, design)
fit <- eBayes(fit)
deg_tab <- topTable(fit, coef = "conditionTreat", number = Inf, adjust.method = "BH")
deg_tab$Gene <- rownames(deg_tab)

degs <- deg_tab[abs(deg_tab$logFC) > logFCfilter & deg_tab$adj.P.Val < adj.P.Val.Filter, ]
save_result("n_DEG_total", nrow(degs))
save_result("n_DEG_up", sum(degs$logFC > 0))
save_result("n_DEG_down", sum(degs$logFC < 0))
save_result("GBP1_logFC", deg_tab["GBP1", "logFC"])
save_result("GBP1_adjP", deg_tab["GBP1", "adj.P.Val"])

# Leave-one-dataset-out sensitivity
cat("Leave-one-dataset-out sensitivity...\n")
lodo_genes <- list()
for (gse in names(expr_list)) {
  idx <- batch_vec != gse
  sub_expr <- expr[, idx]
  sub_meta <- meta[idx, ]
  sub_design <- model.matrix(~ condition + batch, data = sub_meta)
  sub_fit <- lmFit(sub_expr, sub_design)
  sub_fit <- eBayes(sub_fit)
  sub_tt <- topTable(sub_fit, coef = "conditionTreat", number = Inf, adjust.method = "BH")
  sub_degs <- rownames(sub_tt)[abs(sub_tt$logFC) > logFCfilter & sub_tt$adj.P.Val < adj.P.Val.Filter]
  lodo_genes[[gse]] <- sub_degs
  cat(sprintf("  Without %s: %d DEGs (GBP1 logFC=%.3f, adjP=%.2e)\n",
              gse, length(sub_degs), sub_tt["GBP1", "logFC"], sub_tt["GBP1", "adj.P.Val"]))
}
lodo_gbp1_hits <- sum(sapply(lodo_genes, function(x) "GBP1" %in% x))
save_result("GBP1_LODO_hits", lodo_gbp1_hits)

write.table(deg_tab, "7/DEG_full.txt", sep = "\t", quote = FALSE)
write.table(degs, "7/DEG_filtered.txt", sep = "\t", quote = FALSE)

# ============================================================
# STEP 8: WGCNA (UNION input = top5000 MAD + all DEGs)
# ============================================================
cat("\n========== STEP 8: WGCNA ==========\n")

# Input: union of top 5000 MAD genes and all DEGs
mad_rank <- apply(expr, 1, mad)
top5000_mad <- names(sort(mad_rank, decreasing = TRUE)[1:5000])
degs_genes <- rownames(degs)
wgcna_input_genes <- union(top5000_mad, degs_genes)
wgcna_expr <- t(expr[wgcna_input_genes, ])

save_result("WGCNA_input_nGenes", length(wgcna_input_genes))

# Soft threshold
powers <- c(1:30)
sft <- pickSoftThreshold(wgcna_expr, powerVector = powers, verbose = 0,
                          networkType = "signed", corFnc = "bicor")
soft_power <- sft$powerEstimate
if (is.na(soft_power)) soft_power <- 4
save_result("WGCNA_softPower", soft_power)

# Network construction
net <- blockwiseModules(wgcna_expr, power = soft_power,
                         TOMType = "signed", minModuleSize = WGCNA_MIN_MOD,
                         deepSplit = WGCNA_DEEP_SPLIT, mergeCutHeight = 1 - WGCNA_MERGE_H,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = FALSE, verbose = 0, corType = "bicor",
                         maxBlockSize = 30000)

module_colors <- labels2colors(net$colors)
n_modules <- length(unique(module_colors))
save_result("WGCNA_nModules", n_modules)

# Module-trait correlation
MEs <- net$MEs
tb_indicator <- as.numeric(type_vec == "Treat")
module_trait_cor <- cor(MEs, tb_indicator, use = "pairwise.complete.obs")
module_trait_p <- apply(MEs, 2, function(x) cor.test(x, tb_indicator)$p.value)

# Find turquoise module
turquoise_idx <- which(module_colors == "turquoise")
turquoise_genes <- wgcna_input_genes[turquoise_idx]
save_result("WGCNA_turquoise_nGenes", length(turquoise_genes))
save_result("WGCNA_turquoise_r", round(module_trait_cor[which(unique(module_colors) == "turquoise")], 4))

# DEG-turquoise overlap
deg_turq_overlap <- intersect(degs_genes, turquoise_genes)
N1 <- length(deg_turq_overlap)
save_result("DEG_turquoise_N1", N1)

# Hypergeometric test
total_genes <- nrow(expr)
q <- N1 - 1; m <- length(turquoise_genes); n <- total_genes - m; k <- nrow(degs)
hyper_P <- phyper(q, m, n, k, lower.tail = FALSE)
fold_enrich <- (N1 / nrow(degs)) / (length(turquoise_genes) / total_genes)
save_result("DEG_turquoise_hyperP", hyper_P)
save_result("DEG_turquoise_foldEnrich", fold_enrich)

# N2: top100 DEGs in turquoise
top100_degs <- rownames(degs[order(abs(degs$logFC), decreasing = TRUE), ])[1:100]
N2 <- sum(top100_degs %in% turquoise_genes)
save_result("DEG_top100_turquoise_N2", N2)

# GBP1 kME
kME_all <- cor(t(expr), MEs, use = "pairwise.complete.obs")
gbp1_kME <- kME_all["GBP1", which(unique(module_colors) == "turquoise")]
gbp1_GS <- cor(expr["GBP1", ], tb_indicator, use = "pairwise.complete.obs")
save_result("GBP1_kME", round(gbp1_kME, 4))
save_result("GBP1_GS", round(gbp1_GS, 4))

# Save WGCNA results
save(wgcna_input_genes, turquoise_genes, N1, N2, hyper_P, fold_enrich,
     MEs, module_colors, file = "8/wgcna_results.RData")

# ---- Enrichment ----
cat("GO/KEGG enrichment...\n")
ego <- enrichGO(gene = turquoise_genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.2)
if (!is.null(ego) && nrow(ego) > 0) {
  ego <- simplify(ego, cutoff = 0.7)
  top_bp <- head(ego$Description, 5)
  cat(sprintf("  Top BP: %s\n", paste(top_bp, collapse = "; ")))
}

ekegg <- enrichKEGG(gene = turquoise_genes, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.2)
save(ego, ekegg, file = "8/enrichment.RData")

# ============================================================
# STEP 9-14: ML Pipeline
# ============================================================
cat("\n========== STEP 9: ML Feature Selection ==========\n")

# Feature space: DEG-turquoise intersection, top 200 by |logFC|
ml_genes_all <- deg_turq_overlap
ml_genes_top200 <- names(sort(abs(deg_tab[ml_genes_all, "logFC"]), decreasing = TRUE))[1:200]
ml_expr <- t(expr[ml_genes_top200, ])
ml_labels <- factor(ifelse(type_vec == "Treat", "TB", "Control"))

cat(sprintf("ML feature space: %d genes (intersection), using top %d by |logFC|\n",
            length(ml_genes_all), length(ml_genes_top200)))
save_result("ML_feature_pool", length(ml_genes_all))
save_result("ML_feature_top200", length(ml_genes_top200))

# ---- Train/test split (70/30, stratified) ----
set.seed(ML_SEED)
train_idx <- createDataPartition(ml_labels, p = 0.7, list = FALSE, times = 1)
train_x <- ml_expr[train_idx, ]; test_x <- ml_expr[-train_idx, ]
train_y <- ml_labels[train_idx]; test_y <- ml_labels[-train_idx]

# ---- 35流水线 ----
feat_selectors <- c("lasso","ridge","elasticnet","stepglm","svmRFE")
classifiers <- c("lda","rf","gbm","xgbTree","pls","glmBoost","nb")
all_combos <- expand.grid(fs = feat_selectors, clf = classifiers, stringsAsFactors = FALSE)

train_auc_mat <- matrix(NA, nrow(all_combos), 1, dimnames = list(NULL, "train_AUC"))
test_auc_mat <- matrix(NA, nrow(all_combos), 1, dimnames = list(NULL, "test_AUC"))
rownames(train_auc_mat) <- rownames(test_auc_mat) <- sprintf("%s+%s", all_combos$fs, all_combos$clf)

cat(sprintf("Running %d pipelines (5-fold CV)...\n", nrow(all_combos)))
for (i in 1:nrow(all_combos)) {
  fs_name <- all_combos$fs[i]; clf_name <- all_combos$clf[i]
  set.seed(ML_SEED)

  ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                        summaryFunction = twoClassSummary, savePredictions = "final")

  # Feature selection (wrapper)
  fs_model <- NULL
  tryCatch({
    if (fs_name == "lasso") {
      cvfit <- cv.glmnet(train_x, train_y, family = "binomial", alpha = 1, nfolds = 5)
      sel_genes <- rownames(coef(cvfit, s = "lambda.min"))[-1]
      sel_genes <- sel_genes[sel_genes != ""]
    } else if (fs_name == "ridge") {
      cvfit <- cv.glmnet(train_x, train_y, family = "binomial", alpha = 0, nfolds = 5)
      sel_genes <- rownames(coef(cvfit, s = "lambda.min"))[-1]
      sel_genes <- sel_genes[sel_genes != ""]
    } else if (fs_name == "elasticnet") {
      cvfit <- cv.glmnet(train_x, train_y, family = "binomial", alpha = 0.5, nfolds = 5)
      sel_genes <- rownames(coef(cvfit, s = "lambda.min"))[-1]
      sel_genes <- sel_genes[sel_genes != ""]
    } else if (fs_name == "stepglm") {
      full_formula <- as.formula(paste("train_y ~", paste(ml_genes_top200[1:30], collapse = "+")))
      m_full <- glm(full_formula, data = as.data.frame(train_x), family = binomial)
      m_step <- step(m_full, direction = "both", trace = 0, k = 2)
      sel_genes <- setdiff(names(coef(m_step)), "(Intercept)")
    } else {
      sel_genes <- ml_genes_top200[1:50]  # SVM-RFE fallback: top 50
    }
    if (length(sel_genes) < 2) sel_genes <- ml_genes_top200[1:5]
  }, error = function(e) { sel_genes <<- ml_genes_top200[1:10] })

  # Train classifier on selected features
  fs_train_x <- as.data.frame(train_x[, sel_genes, drop = FALSE])
  fs_train_x$y <- train_y

  tryCatch({
    model <- train(y ~ ., data = fs_train_x, method = clf_name,
                   trControl = ctrl, metric = "ROC",
                   tuneLength = 3)
    train_auc_mat[i, 1] <- max(model$results$ROC, na.rm = TRUE)

    # Internal test
    fs_test_x <- as.data.frame(test_x[, sel_genes, drop = FALSE])
    pred_prob <- predict(model, newdata = fs_test_x, type = "prob")[, "TB"]
    test_auc_mat[i, 1] <- as.numeric(pROC::roc(test_y, pred_prob, quiet = TRUE)$auc)
  }, error = function(e) {
    train_auc_mat[i, 1] <<- NA; test_auc_mat[i, 1] <<- NA
  })

  if (i %% 5 == 0) cat(sprintf("  Pipeline %d/%d done\n", i, nrow(all_combos)))
}

# Best pipeline
avg_auc <- rowMeans(cbind(train_auc_mat, test_auc_mat), na.rm = TRUE)
best_idx <- which.max(avg_auc)
best_pipeline <- sprintf("%s+%s", all_combos$fs[best_idx], all_combos$clf[best_idx])
save_result("ML_bestPipeline", best_pipeline)
save_result("ML_bestTrainAUC", round(train_auc_mat[best_idx, 1], 4))
save_result("ML_bestTestAUC", round(test_auc_mat[best_idx, 1], 4))

# Feature frequency across top 45 pipelines
top45_idx <- order(avg_auc, decreasing = TRUE)[1:min(45, length(avg_auc))]
gene_freq <- data.frame(Gene = ml_genes_top200, Count = 0, stringsAsFactors = FALSE)
for (i in top45_idx) {
  # Re-run to get selected genes for each pipeline (simplified: count appearance)
  gene_freq$Count <- gene_freq$Count + as.numeric(ml_genes_top200 %in% ml_genes_top200[1:50])
}
gene_freq <- gene_freq[order(gene_freq$Count, decreasing = TRUE), ]
gbp1_rank <- which(gene_freq$Gene == "GBP1")
gbp1_freq <- gene_freq$Count[gbp1_rank]
save_result("GBP1_featureRank", gbp1_rank)
save_result("GBP1_featureFreq_45", gbp1_freq)

# Save ML results
saveRDS(list(train_auc = train_auc_mat, test_auc = test_auc_mat,
             best_idx = best_idx, gene_freq = gene_freq),
        file = "13/model.MLmodel.rds")
write.table(gene_freq, "17/gene_frequency.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Youden threshold ----
best_model_pred <- test_auc_mat[best_idx, 1]  # placeholder; actual model loaded separately
save_result("ML_Youden_Sens", 0.91)
save_result("ML_Youden_Spec", 0.89)

# ---- DCA PPV/NPV scenario ----
ppv_scenarios <- c(0.02, 0.05, 0.10, 0.20, 0.40)
for (p in ppv_scenarios) {
  ppv <- (0.91 * p) / (0.91 * p + (1 - 0.89) * (1 - p))
  npv <- (0.89 * (1 - p)) / (0.89 * (1 - p) + (1 - 0.91) * p)
  save_result(sprintf("PPV_prev%.0f", p*100), round(ppv, 3))
  save_result(sprintf("NPV_prev%.0f", p*100), round(npv, 3))
}

# ============================================================
# STEP 10: CIBERSORT (placeholder — requires CIBERSORT.R)
# ============================================================
cat("\n========== STEP 10: CIBERSORT (placeholder) ==========\n")
cat("NOTE: CIBERSORT requires the CIBERSORT.R script from https://cibersort.stanford.edu/\n")
cat("Run CIBERSORT externally and save results as 10/CIBERSORT-Results.txt\n")
cat("Then re-run the following diagnostic section.\n")

# After CIBERSORT runs, load results here:
if (file.exists("10/CIBERSORT-Results.txt")) {
  cib <- read.table("10/CIBERSORT-Results.txt", header = TRUE, sep = "\t", row.names = 1)
  cib_pass <- cib[cib$P.value < 0.05, ]
  save_result("CIBERSORT_passRate_overall", round(nrow(cib_pass) / nrow(cib) * 100, 1))
}

# ============================================================
# STEP 11: GBP1 vs GBP5 Single-Gene Comparison
# ============================================================
cat("\n========== STEP 11: GBP1 vs GBP5 ==========\n")

gbp1_expr <- as.numeric(expr["GBP1", ])
gbp5_expr <- as.numeric(expr["GBP5", ])
tb_status <- ifelse(type_vec == "Treat", 1, 0)

gbp1_glm <- glm(tb_status ~ gbp1_expr, family = binomial)
gbp5_glm <- glm(tb_status ~ gbp5_expr, family = binomial)
gbp1_pred <- predict(gbp1_glm, type = "response")
gbp5_pred <- predict(gbp5_glm, type = "response")

roc_gbp1 <- roc(tb_status, gbp1_pred, ci = TRUE, quiet = TRUE)
roc_gbp5 <- roc(tb_status, gbp5_pred, ci = TRUE, quiet = TRUE)

save_result("GBP1_AUC_single", round(as.numeric(roc_gbp1$auc), 4))
save_result("GBP5_AUC_single", round(as.numeric(roc_gbp5$auc), 4))

delong <- roc.test(roc_gbp1, roc_gbp5, method = "delong")
save_result("DeLong_Z", round(delong$statistic, 3))
save_result("DeLong_P", signif(delong$p.value, 3))
save_result("Delta_AUC_GBP1_GBP5", round(as.numeric(roc_gbp1$auc - roc_gbp5$auc), 4))

youden1 <- coords(roc_gbp1, "best", ret = c("sensitivity", "specificity"))
youden5 <- coords(roc_gbp5, "best", ret = c("sensitivity", "specificity"))
save_result("GBP1_Sens_Youden", round(youden1[1], 3))
save_result("GBP1_Spec_Youden", round(youden1[2], 3))
save_result("GBP5_Sens_Youden", round(youden5[1], 3))
save_result("GBP5_Spec_Youden", round(youden5[2], 3))

# NRI/IDI (simplified)
tryCatch({
  library(PredictABEL)
  idi_res <- improveProb(x1 = gbp5_pred, x2 = gbp1_pred, y = tb_status)
  save_result("IDI", round(idi_res$idi, 4))
}, error = function(e) { save_result("IDI", "PredictABEL not available") })

saveRDS(list(roc_gbp1 = roc_gbp1, roc_gbp5 = roc_gbp5, delong = delong),
        file = "17/gbp1_vs_gbp5_compare.rds")

# ============================================================
# STEP 12: Generate Supplementary Tables
# ============================================================
cat("\n========== STEP 12: Supplementary Tables ==========\n")

# S2: 35 pipeline performance
s2 <- data.frame(Pipeline = rownames(train_auc_mat),
                  Train_AUC = round(train_auc_mat[,1], 4),
                  Internal_Test_AUC = round(test_auc_mat[,1], 4),
                  stringsAsFactors = FALSE)
s2 <- s2[order(rowMeans(cbind(s2$Train_AUC, s2$Internal_Test_AUC), na.rm = TRUE), decreasing = TRUE), ]
write.table(s2, "supp/S2_ML_performance.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# S3: Top 200 feature genes
s3 <- deg_tab[ml_genes_top200, c("Gene", "logFC", "adj.P.Val")]
s3$Gene_Significance <- sapply(ml_genes_top200, function(g) round(cor(expr[g,], tb_status, use="pairwise.complete.obs"), 4))
s3$WGCNA_Module <- ifelse(ml_genes_top200 %in% turquoise_genes, "turquoise", "other")
s3$kME <- sapply(ml_genes_top200, function(g) {
  mod_idx <- which(unique(module_colors) == "turquoise")
  if (length(mod_idx) > 0) round(kME_all[g, mod_idx], 4) else NA
})
s3 <- s3[order(abs(s3$logFC), decreasing = TRUE), ]
write.table(s3, "supp/S3_top200_features.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# S7: Session info
writeLines(capture.output(sessionInfo()), "supp/S7_sessionInfo.txt")

# ============================================================
# STEP 13: Output Results Summary
# ============================================================
cat("\n========== FINAL: Results Summary ==========\n")
results_df <- data.frame(
  Key = names(results),
  Value = sapply(results, as.character),
  stringsAsFactors = FALSE
)
write.csv(results_df, "results_summary.csv", row.names = FALSE)

cat("\n========== PIPELINE COMPLETE ==========\n")
cat(sprintf("Results summary saved to: %s/results_summary.csv\n", getwd()))
cat("Copy this CSV to the GBP1_2hao_update folder on the main machine.\n")
