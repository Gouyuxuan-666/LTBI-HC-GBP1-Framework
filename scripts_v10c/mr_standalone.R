###############################################################################
# GBP1 Mendelian Randomization: cis-eQTL -> TB risk
# 本机运行 (需要外网访问 OpenGWAS + eQTLGen)
###############################################################################
setwd("C:/Users/1/Desktop/GBP1_2hao_update")

# ---- Proxy setup for R ----
options(internet.proxy = "http://127.0.0.1:7993")
Sys.setenv(http_proxy = "http://127.0.0.1:7993")
Sys.setenv(https_proxy = "http://127.0.0.1:7993")
Sys.setenv(HTTP_PROXY = "http://127.0.0.1:7993")
Sys.setenv(HTTPS_PROXY = "http://127.0.0.1:7993")
httr::set_config(httr::use_proxy("http://127.0.0.1:7993"))

# OpenGWAS JWT authentication (required since May 2024)
# Generate a token from https://api.opengwas.io/
token <- tryCatch(ieugwasr::get_opengwas_jwt(), error = function(e) {
  cat("JWT failed. Trying without authentication...\n")
  NULL
})
if (!is.null(token)) ieugwasr::user(token)

library(httr)
library(TwoSampleMR)
library(ieugwasr)
library(ggplot2)

dir.create("mr_output", showWarnings = FALSE)
results <- list()

# ================================================================
# Step 1: Get GBP1 cis-eQTL instruments
# ================================================================
cat("\n=== Step 1: Finding GBP1 cis-eQTL instruments ===\n")

# GBP1 Ensembl ID: ENSG00000117228
# Try GTEx v8 whole blood eQTL first

# Method A: Query IEU OpenGWAS for eQTL data
gbp1_eqtl_exposure <- tryCatch({
  # GTEx v8 Whole Blood eQTL
  extract_instruments(outcomes = "eqtl-a-ENSG00000117228", p1 = 5e-6, clump = TRUE, r2 = 0.01, kb = 1000)
}, error = function(e) {
  cat("  OpenGWAS eQTL failed:", e$message, "\n")
  NULL
})

# Method B: If above fails, try broader search
if (is.null(gbp1_eqtl_exposure) || nrow(gbp1_eqtl_exposure) < 1) {
  cat("  Trying alternative eQTL sources...\n")
  # Query available eQTL datasets for GBP1
  gbp1_eqtl_exposure <- tryCatch({
    extract_instruments(outcomes = "eqtl-a-ENSG00000117228", p1 = 1e-4, clump = TRUE, r2 = 0.01, kb = 1000)
  }, error = function(e) {
    cat("  All eQTL queries failed.\n")
    NULL
  })
}

if (!is.null(gbp1_eqtl_exposure) && nrow(gbp1_eqtl_exposure) >= 1) {
  cat(sprintf("  Found %d GBP1 eQTL instruments\n", nrow(gbp1_eqtl_exposure)))
  # Ensure required columns
  gbp1_eqtl_exposure$exposure <- "GBP1 expression"
  gbp1_eqtl_exposure$id.exposure <- "eqtl-a-ENSG00000117228"
} else {
  cat("  WARNING: No GBP1 eQTL instruments found.\n")
  cat("  MR not feasible with current data sources.\n")

  # Save status
  results$MR_status <- "No GBP1 eQTL instruments found"
  write.csv(data.frame(Key = names(results), Value = unlist(results)),
            "mr_output/mr_status.csv", row.names = FALSE)
  quit(save = "no", status = 0)
}

# ================================================================
# Step 2: Extract TB GWAS outcome data
# ================================================================
cat("\n=== Step 2: Extracting TB outcome data ===\n")

# Try multiple TB GWAS datasets
tb_outcomes <- c(
  "finn-b-AB1_TUBERCULOSIS",      # FinnGen AB1 Tuberculosis
  "ebi-a-GCST90018892",            # TB GWAS
  "ieu-b-4972"                     # UK Biobank TB
)

tb_dat <- NULL
for (outcome_id in tb_outcomes) {
  tb_dat <- tryCatch({
    extract_outcome_data(snps = gbp1_eqtl_exposure$SNP, outcomes = outcome_id)
  }, error = function(e) {
    cat(sprintf("  %s: %s\n", outcome_id, e$message))
    NULL
  })
  if (!is.null(tb_dat) && nrow(tb_dat) >= 1) {
    cat(sprintf("  Success: %s (%d SNPs matched)\n", outcome_id, nrow(tb_dat)))
    break
  }
}

if (is.null(tb_dat) || nrow(tb_dat) < 1) {
  cat("  ERROR: No TB outcome data available for these SNPs.\n")
  results$MR_status <- "No TB outcome data"
  write.csv(data.frame(Key = names(results), Value = unlist(results)),
            "mr_output/mr_status.csv", row.names = FALSE)
  quit(save = "no", status = 0)
}

# ================================================================
# Step 3: Harmonise and run MR
# ================================================================
cat("\n=== Step 3: Harmonising and running MR ===\n")

dat <- harmonise_data(exposure_dat = gbp1_eqtl_exposure, outcome_dat = tb_dat)
cat(sprintf("  %d SNPs after harmonisation\n", nrow(dat)))

if (nrow(dat) < 1) {
  cat("  ERROR: No SNPs survived harmonisation.\n")
  quit(save = "no", status = 0)
}

# MR analysis
mr_res <- mr(dat, method_list = c("mr_wald_ratio", "mr_ivw", "mr_egger_regression",
                                   "mr_weighted_median", "mr_weighted_mode"))
cat("\n=== MR Results ===\n")
print(mr_res)

# Sensitivity
het <- mr_heterogeneity(dat)
pleio <- mr_pleiotropy_test(dat)
cat("\nHeterogeneity:\n"); print(het)
cat("\nPleiotropy:\n"); print(pleio)

# ================================================================
# Step 4: Generate plots
# ================================================================
cat("\n=== Step 4: Generating plots ===\n")

if (nrow(mr_res) > 1) {
  # Forest plot
  tryCatch({
    p1 <- mr_forest_plot(mr_res)
    ggsave("mr_output/Fig_MR_Forest.pdf", p1[[1]], width = 8, height = 5)
    cat("  Forest plot saved\n")
  }, error = function(e) cat("  Forest plot failed:", e$message, "\n"))

  # Scatter plot
  tryCatch({
    p2 <- mr_scatter_plot(mr_res, dat)
    ggsave("mr_output/Fig_MR_Scatter.pdf", p2[[1]], width = 7, height = 6)
    cat("  Scatter plot saved\n")
  }, error = function(e) cat("  Scatter plot failed:", e$message, "\n"))
}

# Funnel plot
if (nrow(dat) >= 3) {
  tryCatch({
    singlesnp <- mr_singlesnp(dat)
    p3 <- mr_funnel_plot(singlesnp)
    ggsave("mr_output/Fig_MR_Funnel.pdf", p3[[1]], width = 7, height = 5)
    cat("  Funnel plot saved\n")
  }, error = function(e) cat("  Funnel plot failed:", e$message, "\n"))
}

# Leave-one-out
if (nrow(dat) >= 3) {
  tryCatch({
    loo <- mr_leaveoneout(dat)
    p4 <- mr_leaveoneout_plot(loo)
    ggsave("mr_output/Fig_MR_LeaveOneOut.pdf", p4[[1]], width = 8, height = 6)
    cat("  Leave-one-out plot saved\n")
  }, error = function(e) cat("  Leave-one-out plot failed:", e$message, "\n"))
}

# ================================================================
# Step 5: Save results
# ================================================================
write.csv(mr_res, "mr_output/mr_results.csv", row.names = FALSE)
write.csv(dat, "mr_output/mr_harmonised_data.csv", row.names = FALSE)

results$MR_IVW_beta <- round(mr_res$b[mr_res$method == "Inverse variance weighted"], 4)
results$MR_IVW_P <- signif(mr_res$pval[mr_res$method == "Inverse variance weighted"], 3)
results$MR_IVW_OR <- round(exp(mr_res$b[mr_res$method == "Inverse variance weighted"]), 3)
results$MR_nSNPs <- nrow(dat)
results$MR_Egger_intercept <- round(pleio$egger_intercept, 4)
results$MR_Egger_P <- signif(pleio$pval, 3)
results$MR_status <- "Complete"

write.csv(data.frame(Key = names(results), Value = unlist(results)),
          "mr_output/mr_status.csv", row.names = FALSE)

cat(sprintf("\n=== MR Complete ===\n"))
cat(sprintf("IVW: beta=%.4f, P=%.3f, OR=%.3f\n", results$MR_IVW_beta, results$MR_IVW_P, results$MR_IVW_OR))
