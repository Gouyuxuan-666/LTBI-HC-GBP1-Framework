###############################################################################
# GBP1-TB External Validation Script
# Datasets: GSE39940 (334 samples, African children, GPL10558)
#           GSE19444 (54 samples, UK adults, GPL6947, Berry Nature 2010)
# Run AFTER run_all.R completes: Rscript external_validation.R
# Uses GEOquery for fully reproducible download (no local data dependency)
###############################################################################

setwd("F:/GBP1_pipeline_2hao")
# ---- Setup ----
options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))

if (!requireNamespace("GEOquery", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install("GEOquery", update=FALSE, ask=FALSE)
}
suppressMessages({
  library(GEOquery); library(limma); library(pROC)
  library(illuminaHumanv4.db)  # GPL10558 probe mapping
})

# Load trained model from individual files (no dependency on pipeline_results.rds)
source("refer.ML.R")

model <- readRDS("13/model.MLmodel.rds")
auc_tab <- read.table("14/model.AUCmatrix.txt", header=TRUE, sep="\t", row.names=1)
gene_freq <- read.table("17/gene_frequency.txt", header=TRUE)
disease_genes <- head(gene_freq$Gene, 200)

avg_auc <- rowMeans(auc_tab, na.rm=TRUE)
best_method <- names(which.max(avg_auc))
best_model <- model[[best_method]]
cat(sprintf("Best model: %s | Disease genes: %d\n", best_method, length(disease_genes)))

# ---- Function: parse GEO series matrix WITHOUT internet ----
parse_series_matrix <- function(gzfile, platform_pkg) {
  cat(sprintf("  Reading %s...\n", basename(gzfile)))
  lines <- readLines(gzfile)

  # GSM IDs: one line with all IDs quoted
  gsm_line <- lines[grep("^!Sample_geo_accession", lines)][1]
  gsm_parts <- strsplit(gsm_line, '\t')[[1]]
  gsm_ids <- gsub('"', '', gsm_parts[-1])
  gsm_ids <- trimws(gsm_ids)
  n_samples <- length(gsm_ids)
  cat(sprintf("  Found %d samples\n", n_samples))

  # Parse disease group from characteristics (tab-separated on single line!)
  char_lines <- lines[grep("^!Sample_characteristics_ch1", lines)]
  all_char_vals <- c()
  for (cl in char_lines) {
    parts <- strsplit(cl, '\t')[[1]]
    vals <- gsub('"', '', parts[-1])
    vals <- trimws(vals)
    all_char_vals <- rbind(all_char_vals, vals)
  }
  if (length(char_lines) == 1) all_char_vals <- matrix(all_char_vals, nrow=1)

  # Parse from Sample_title too (also tab-separated)
  title_line <- lines[grep("^!Sample_title", lines)][1]
  title_parts <- strsplit(title_line, '\t')[[1]]
  titles <- trimws(gsub('"', '', title_parts[-1]))

  is_tb <- rep(FALSE, n_samples)
  is_control <- rep(FALSE, n_samples)

  for (i in 1:n_samples) {
    # Build full description from all available fields
    fields <- tolower(c(titles[i], all_char_vals[, i]))
    full <- paste(fields, collapse=" ")
    cat(sprintf("  Sample %d: %s\n", i, substr(full, 1, 80)))

    is_tb[i] <- grepl("active tb|active tuberculosis|pulmonary tuberculosis|tb patient|tuberculosis", full) &
      !grepl("latent", full)
    is_control[i] <- grepl("healthy|normal|control", full) &
      !grepl("latent|ltbi|sarcoidosis|pneumonia|cancer|other disease", full)
  }

  # Parse matrix data
  tbl_start <- which(lines == "!series_matrix_table_begin") + 1
  tbl_end <- which(lines == "!series_matrix_table_end") - 1
  cat(sprintf("  Parsing matrix rows %d-%d...\n", tbl_start, tbl_end))

  mat_lines <- lines[tbl_start:tbl_end]
  header <- strsplit(mat_lines[1], "\t")[[1]]
  header <- gsub('"', '', header)

  probe_ids <- c()
  exp_list <- list()
  for (i in 2:length(mat_lines)) {
    parts <- strsplit(mat_lines[i], "\t")[[1]]
    pid <- gsub('"', '', parts[1])
    vals <- as.numeric(gsub('"', '', parts[-1]))
    probe_ids <- c(probe_ids, pid)
    exp_list[[pid]] <- vals
  }

  # Build expression matrix
  exp_mat <- do.call(rbind, exp_list)
  colnames(exp_mat) <- header[-1]
  rownames(exp_mat) <- probe_ids
  cat(sprintf("  Matrix: %d probes x %d samples\n", nrow(exp_mat), ncol(exp_mat)))

  # Probe-to-gene mapping via Bioconductor
  if (platform_pkg == "illuminaHumanv4.db") {
    mapped <- select(illuminaHumanv4.db, keys=rownames(exp_mat), columns="SYMBOL", keytype="PROBEID")
  } else {
    library(illuminaHumanv3.db)
    mapped <- select(illuminaHumanv3.db, keys=rownames(exp_mat), columns="SYMBOL", keytype="PROBEID")
  }
  mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
  exp_mat <- exp_mat[mapped$PROBEID, , drop=FALSE]
  rownames(exp_mat) <- mapped$SYMBOL
  exp_mat <- avereps(exp_mat)

  cat(sprintf("  Active TB: %d, Healthy: %d\n", sum(is_tb), sum(is_control)))

  keep <- is_tb | is_control
  list(expr=exp_mat[, keep, drop=FALSE], label=ifelse(is_tb[keep], 1, 0),
       n_tb=sum(is_tb), n_hc=sum(is_control))
}

# ---- Process validation datasets ----
process_geo <- function(gse_id, platform_pkg, disease_genes) {
  cat(sprintf("\n--- Processing %s ---\n", gse_id))
  gzfile <- file.path(getwd(), paste0(gse_id, "_series_matrix.txt.gz"))
  if (!file.exists(gzfile)) stop("File not found: ", gzfile)

  res <- parse_series_matrix(gzfile, platform_pkg)
  common_genes <- intersect(disease_genes, rownames(res$expr))
  cat(sprintf("  Common disease genes: %d / %d\n", length(common_genes), length(disease_genes)))
  res$expr <- res$expr[common_genes, , drop=FALSE]
  res$gse_id <- gse_id
  res
}

# ---- Download & Process Validation Datasets ----
dir.create("extval", showWarnings=FALSE)

# Dataset 1: GSE39940 (GPL10558, African children)
if (!requireNamespace("illuminaHumanv3.db", quietly=TRUE)) {
  BiocManager::install("illuminaHumanv3.db", update=FALSE, ask=FALSE)
}

val1 <- tryCatch(
  process_geo("GSE39940", "illuminaHumanv4.db", disease_genes),
  error=function(e) { cat(sprintf("GSE39940 failed: %s\n", e$message)); NULL }
)

# Dataset 2: GSE19444 (GPL6947, Berry Nature 2010 validation)
val2 <- tryCatch(
  process_geo("GSE19444", "illuminaHumanv3.db", disease_genes),
  error=function(e) { cat(sprintf("GSE19444 failed: %s\n", e$message)); NULL }
)

# ---- Apply ML Model & Calculate AUC ----
all_results <- list()
auc_matrix <- data.frame(row.names=names(model))

for (val_name in c("val1", "val2")) {
  val <- get(val_name)
  if (is.null(val)) next

  cat(sprintf("\n========== Validating on %s ==========\n", val$gse_id))
  cat(sprintf("Samples: %d (TB:%d, HC:%d)\n", length(val$label), val$n_tb, val$n_hc))

  # Scale validation data
  val_mat <- scaleData(t(val$expr), centerFlags=TRUE, scaleFlags=TRUE)
  val_lab <- data.frame(Type=val$label, Cohort=val$gse_id, row.names=colnames(val$expr))

  # Apply best model and calculate AUC
  for (method in names(model)) {
    fit <- model[[method]]
    feat <- ExtractVar(fit)
    feat <- intersect(feat, colnames(val_mat))
    if (length(feat) < 2) { auc_matrix[method, val$gse_id] <- NA; next }

    auc_val <- tryCatch({
      rs <- CalPredictScore(fit, val_mat[, feat, drop=FALSE])
      as.numeric(auc(roc(val_lab$Type, rs[rownames(val_mat)])))
    }, error=function(e) NA)
    auc_matrix[method, val$gse_id] <- auc_val
  }

  # Save results
  all_results[[val$gse_id]] <- list(
    gse_id=val$gse_id,
    n_samples=length(val$label),
    n_tb=val$n_tb,
    n_hc=val$n_hc,
    auc=auc_matrix[best_method, val$gse_id]
  )
  cat(sprintf("Best model (%s) AUC: %.3f\n", best_method,
    auc_matrix[best_method, val$gse_id]))

  # ROC curve for best model
  fit <- model[[best_method]]
  feat <- ExtractVar(fit)
  feat <- intersect(feat, colnames(val_mat))
  rs <- CalPredictScore(fit, val_mat[, feat, drop=FALSE])

  pdf(sprintf("extval/ROC_%s.pdf", val$gse_id), width=7, height=6)
  roc_obj <- roc(val_lab$Type, rs[rownames(val_mat)])
  plot.roc(roc_obj, col="#DE2D26", lwd=2.5, legacy.axes=TRUE,
    print.auc=TRUE, print.auc.cex=1.4, auc.polygon=TRUE,
    auc.polygon.col=adjustcolor("#DE2D26", 0.2),
    main=paste0("External Validation: ", val$gse_id))
  dev.off()
  cat(sprintf("  ROC_%s.pdf saved\n", val$gse_id))
}

# ---- Save Summary ----
valid_sets <- names(all_results)
cat(sprintf("\n========== EXTERNAL VALIDATION SUMMARY ==========\n"))
for (s in valid_sets) {
  r <- all_results[[s]]
  cat(sprintf("%s: %d samples (TB:%d, HC:%d), Best AUC=%.3f\n",
    r$gse_id, r$n_samples, r$n_tb, r$n_hc, r$auc))
}

write.table(cbind(Method=rownames(auc_matrix), auc_matrix),
  "extval/external_validation_AUC.txt", sep="\t", quote=FALSE, row.names=FALSE)

saveRDS(list(results=all_results, auc_matrix=auc_matrix, best_method=best_method),
  "extval/external_validation.rds")
cat("\nSaved extval/external_validation.rds\n")
