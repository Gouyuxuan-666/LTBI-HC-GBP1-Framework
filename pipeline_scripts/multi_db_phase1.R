###############################################################################
# GBP1 Multi-Database Phase 1: GTEx + HPA + GWAS Catalog
# Three API queries, no downloads needed (~1 min total)
# Reference repos: ropensci/gtexr, anhtr/HPAanalyze, ramiromagno/gwasrapidd
# Run AFTER run_all.R: Rscript multi_db_phase1.R
###############################################################################

options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))

# ---- Install & Load ----
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")

pkg_list <- c("ggplot2","dplyr","reshape2","ggpubr")
for (p in pkg_list) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
}

if (!requireNamespace("gtexr", quietly=TRUE)) install.packages("gtexr")
if (!requireNamespace("HPAanalyze", quietly=TRUE)) BiocManager::install("HPAanalyze", update=FALSE, ask=FALSE)
if (!requireNamespace("gwasrapidd", quietly=TRUE)) install.packages("gwasrapidd")

suppressMessages({
  library(ggplot2); library(dplyr); library(gtexr)
  library(HPAanalyze); library(gwasrapidd)
})

# ---- Output dirs ----
for (d in c("extval/db_GTEx","extval/db_HPA","extval/db_GWAS")) dir.create(d, recursive=TRUE, showWarnings=FALSE)

# ============================================================================
# 1. GTEx: GBP1 median TPM across 54 human tissues
# ============================================================================
cat("\n========== GTEx: GBP1 Tissue Expression ==========\n")

gbp1_ensg <- "ENSG00000166128"  # GBP1 Ensembl ID

tryCatch({
  gtex_expr <- get_median_gene_expression(gtexr:::gtex_db_to_dataset_ids("gtex_v8")[[1]],
    gencodeIds=gbp1_ensg)
  gtex_expr <- gtex_expr[order(gtex_expr$median, decreasing=TRUE), ]
  gtex_expr$tissueSiteDetailId <- factor(gtex_expr$tissueSiteDetailId,
    levels=rev(gtex_expr$tissueSiteDetailId))

  write.table(gtex_expr, "extval/db_GTEx/GBP1_GTEx_tissues.txt", sep="\t", quote=FALSE, row.names=FALSE)

  pdf("extval/db_GTEx/GBP1_GTEx_barplot.pdf", width=10, height=12)
  print(ggplot(gtex_expr, aes(x=tissueSiteDetailId, y=median)) +
    geom_bar(stat="identity", fill="#3182BD", alpha=0.85) + coord_flip() +
    labs(title="GBP1 Expression Across 54 GTEx Tissues",
         subtitle="GTEx v8 Median TPM", x="", y="Median TPM") +
    theme_bw(11) + theme(axis.text.y=element_text(size=7)))
  dev.off()
  cat(sprintf("GTEx: %d tissues, max=%.1f TPM (%s)\n",
    nrow(gtex_expr), max(gtex_expr$median), gtex_expr$tissueSiteDetailId[1]))
}, error=function(e) {
  cat(sprintf("GTEx failed: %s\n", e$message))
  cat("Trying fallback: gtexr utilities...\n")
  tryCatch({
    datasets <- gtexr::get_datasets()
    cat(sprintf("Available GTEx datasets: %d\n", nrow(datasets)))
  }, error=function(e2) cat(sprintf("GTEx fallback also failed: %s\n", e2$message)))
})

# ============================================================================
# 2. HPA: GBP1 protein expression + IHC evidence
# ============================================================================
cat("\n========== HPA: GBP1 Protein Expression ==========\n")

tryCatch({
  hpa_data <- hpaDownload(downloadList="all", version="latest", archive=FALSE)

  # RNA tissue expression
  hpa_rna <- hpaSubset(data=hpa_data$rna_tissue, targetGene="GBP1")
  if (nrow(hpa_rna) > 0) {
    hpa_rna <- hpa_rna[order(hpa_rna$nTPM, decreasing=TRUE), ]
    write.table(hpa_rna, "extval/db_HPA/GBP1_HPA_RNA.txt", sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("HPA RNA: %d tissues, max=%.1f nTPM\n", nrow(hpa_rna), max(hpa_rna$nTPM, na.rm=TRUE)))

    pdf("extval/db_HPA/GBP1_HPA_tissue_barplot.pdf", width=10, height=8)
    hpa_plot <- head(hpa_rna, 30)
    hpa_plot$Tissue <- factor(hpa_plot$Tissue, levels=rev(hpa_plot$Tissue))
    print(ggplot(hpa_plot, aes(x=Tissue, y=nTPM)) +
      geom_bar(stat="identity", fill="#DE2D26", alpha=0.85) + coord_flip() +
      labs(title="GBP1 RNA Expression (HPA Consensus)", x="", y="Normalized TPM") +
      theme_bw(11))
    dev.off()
  }

  # Normal tissue IHC
  hpa_ihc <- hpaSubset(data=hpa_data$normal_tissue, targetGene="GBP1")
  if (nrow(hpa_ihc) > 0) {
    write.table(hpa_ihc, "extval/db_HPA/GBP1_HPA_IHC.txt", sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("HPA IHC: %d tissue entries\n", nrow(hpa_ihc)))

    # Summarize staining levels
    stain_sum <- table(hpa_ihc$Level, hpa_ihc$Tissue)
    if (ncol(stain_sum) > 0) {
      stain_df <- as.data.frame.matrix(stain_sum)
      cat("Staining levels across tissues:\n")
      print(stain_df[, 1:min(8, ncol(stain_df))])
    }
  }

  # Pathology (cancer) IHC
  hpa_patho <- hpaSubset(data=hpa_data$pathology, targetGene="GBP1")
  if (nrow(hpa_patho) > 0) {
    write.table(hpa_patho, "extval/db_HPA/GBP1_HPA_pathology.txt", sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("HPA pathology: %d cancer samples\n", nrow(hpa_patho)))
  }

}, error=function(e) cat(sprintf("HPA failed: %s\n", e$message)))

# ============================================================================
# 3. GWAS Catalog: SNPs and traits near GBP1 locus
# ============================================================================
cat("\n========== GWAS Catalog: GBP1 Locus ==========\n")

tryCatch({
  # Query variants mapped to GBP1 gene
  gwas_vars <- get_variants(gene_name="GBP1")
  if (length(gwas_vars@variants) > 0 && nrow(gwas_vars@variants) > 0) {
    cat(sprintf("Variants near GBP1: %d\n", nrow(gwas_vars@variants)))

    # Get associations for these variants
    assoc_ids <- unique(unlist(gwas_vars@variant_associations))
    if (length(assoc_ids) > 0) {
      gwas_assoc <- get_associations(association_id=assoc_ids[1:min(100, length(assoc_ids))])

      if (length(gwas_assoc@associations) > 0) {
        gw_df <- gwas_assoc@associations
        gw_df <- gw_df[order(gw_df$pvalue), ]

        # Extract trait info
        trait_df <- gwas_assoc@traits
        risk_df <- gwas_assoc@risk_alleles

        # Merge key columns
        gw_out <- gw_df[, c("association_id","pvalue","pvalue_description"), drop=FALSE]
        if (nrow(trait_df) > 0) {
          trait_s <- sapply(gw_assoc@associations$association_id, function(aid) {
            rows <- trait_df[trait_df$association_id == aid, ]
            if (nrow(rows) > 0) paste(unique(rows$trait), collapse="; ") else ""
          })
          gw_out$Traits <- trait_s
        }

        gw_out <- gw_out[order(gw_out$pvalue), ]
        write.table(gw_out, "extval/db_GWAS/GBP1_GWAS_associations.txt", sep="\t", quote=FALSE, row.names=FALSE)
        cat(sprintf("GWAS associations: %d (min p=%.2e)\n", nrow(gw_out), min(gw_out$pvalue, na.rm=TRUE)))
        cat("Top 10 traits:\n")
        print(head(gw_out, 10))
      }
    }
  } else {
    cat("No GWAS variants found for GBP1.\n")
    cat("Trying genomic region query (chr1:89.05-89.07 Mb)...\n")
    gwas_vars2 <- get_variants(genomic_range=list(
      chromosome="1", start=89050000, end=89070000))
    if (length(gwas_vars2@variants) > 0) {
      cat(sprintf("Region variants: %d\n", nrow(gwas_vars2@variants)))
    }
  }
}, error=function(e) cat(sprintf("GWAS failed: %s\n", e$message)))

# ============================================================================
# Save Summary
# ============================================================================
cat("\n========== Multi-DB Phase 1 Complete ==========\n")
cat("Outputs: extval/db_GTEx/ extval/db_HPA/ extval/db_GWAS/\n")
