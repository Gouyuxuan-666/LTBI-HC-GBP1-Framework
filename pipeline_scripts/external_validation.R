###############################################################################
# GBP1-TB External Validation Script
# Datasets: GSE39940 (334 samples, African children, GPL10558)
#           GSE19444 (54 samples, UK adults, GPL6947, Berry Nature 2010)
# Run AFTER run_all.R completes: Rscript external_validation.R
# Uses GEOquery for fully reproducible download (no local data dependency)
###############################################################################

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

# ---- Function: process a GEO dataset into gene expression matrix ---
process_geo <- function(gse_id, platform_pkg, disease_genes) {
  cat(sprintf("\n--- Processing %s ---\n", gse_id))
  gse <- getGEO(gse_id, GSEMatrix=TRUE, AnnotGPL=TRUE, getGPL=TRUE)
  eset <- gse[[1]]

  # Probe-to-gene mapping using Bioconductor annotation package
  probes <- rownames(exprs(eset))
  if (platform_pkg == "illuminaHumanv4.db") {
    mapped <- select(illuminaHumanv4.db, keys=probes, columns="SYMBOL", keytype="PROBEID")
  } else if (platform_pkg == "illuminaHumanv3.db") {
    library(illuminaHumanv3.db)
    mapped <- select(illuminaHumanv3.db, keys=probes, columns="SYMBOL", keytype="PROBEID")
  }

  # Merge probes to gene symbols (max probe per gene)
  mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
  exp_dat <- exprs(eset)[mapped$PROBEID, , drop=FALSE]
  rownames(exp_dat) <- mapped$SYMBOL
  exp_dat <- avereps(exp_dat)

  # Extract disease genes
  common_genes <- intersect(disease_genes, rownames(exp_dat))
  cat(sprintf("  Common disease genes: %d / %d\n", length(common_genes), length(disease_genes)))
  exp_dat <- exp_dat[common_genes, , drop=FALSE]

  # Parse sample groups from phenotype
  pd <- pData(eset)
  cat(sprintf("  Phenotype columns: %s\n", paste(colnames(pd), collapse=", ")))

  # Try to auto-detect group column
  grp <- NULL
  for (col in c("characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2",
                "source_name_ch1", "title")) {
    if (!is.null(pd[[col]])) {
      vals <- tolower(pd[[col]])
      if (any(grepl("tb|tuberculosis|active|control|healthy|normal|latent|ltbi", vals))) {
        grp <- pd[[col]]
        cat(sprintf("  Using column '%s' for groups\n", col))
        break
      }
    }
  }
  if (is.null(grp)) {
    cat("  WARNING: Could not auto-detect groups. Using first characteristics column.\n")
    grp <- pd[, grep("characteristics", colnames(pd))[1]]
  }

  # Classify: Active TB vs Healthy Control (exclude LTBI, other diseases)
  grp_lower <- tolower(grp)
  is_tb <- grepl("active tb|active tuberculosis|tuberculosis[^a-z]|pulmonary tb|pulmonary tuberculosis|tb patient", grp_lower)
  is_control <- grepl("healthy control|healthy|normal|control[^a-z]|hc[^a-z]", grp_lower) &
    !grepl("latent|ltbi|sarcoidosis|pneumonia|lung cancer|cancer|other disease|od[^a-z]", grp_lower)

  cat(sprintf("  Active TB: %d, Healthy Control: %d\n", sum(is_tb), sum(is_control)))

  keep <- is_tb | is_control
  if (sum(keep) == 0) stop("No Active TB or Control samples found")
  if (sum(is_tb) == 0) warning("No Active TB samples found in ", gse_id)
  if (sum(is_control) == 0) warning("No control samples found in ", gse_id)

  exp_dat <- exp_dat[, keep, drop=FALSE]
  label <- ifelse(is_tb[keep], 1, 0)

  return(list(expr=exp_dat, label=label, gse_id=gse_id, n_tb=sum(is_tb), n_hc=sum(is_control)))
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
