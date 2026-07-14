###############################################################################
# GBP1 Multi-Database Phase 1: GTEx + HPA + GWAS Catalog
# Direct REST API via httr — no dependency on gtexr/HPAanalyze/gwasrapidd
###############################################################################
setwd("F:/GBP1_pipeline_2hao")
options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))

for (p in c("httr","jsonlite","ggplot2","dplyr")) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
}
suppressMessages({ library(httr); library(jsonlite); library(ggplot2); library(dplyr) })

for (d in c("extval/db_GTEx","extval/db_HPA","extval/db_GWAS")) dir.create(d, recursive=TRUE, showWarnings=FALSE)

# ============================================================================
# 1. GTEx v8: median TPM per tissue via Portal API v2
# ============================================================================
cat("\n========== GTEx: GBP1 Tissue Expression ==========\n")
tryCatch({
  resp <- GET("https://gtexportal.org/rest/v2/expression/medianTranscriptExpression",
    query=list(gencodeId="ENSG00000166128.13", datasetId="gtex_v8", format="json"))
  if (status_code(resp) == 200) {
    gtex <- fromJSON(content(resp, "text", encoding="UTF-8"))
    gtex_df <- gtex$medianTranscriptExpression
    if (!is.null(gtex_df) && nrow(gtex_df) > 0) {
      gtex_df <- gtex_df[order(gtex_df$median, decreasing=TRUE), ]
      write.table(gtex_df, "extval/db_GTEx/GBP1_GTEx_tissues.txt", sep="\t", quote=FALSE, row.names=FALSE)
      cat(sprintf("GTEx: %d tissues, max=%.1f TPM (%s)\n", nrow(gtex_df),
        max(gtex_df$median), gtex_df$tissueSiteDetailId[1]))
      # Barplot
      top30 <- head(gtex_df, 30)
      top30$tissueSiteDetailId <- factor(top30$tissueSiteDetailId,
        levels=rev(top30$tissueSiteDetailId))
      pdf("extval/db_GTEx/GBP1_GTEx_barplot.pdf", width=10, height=8)
      print(ggplot(top30, aes(x=tissueSiteDetailId, y=median)) +
        geom_bar(stat="identity", fill="#3182BD") + coord_flip() +
        labs(title="GBP1 Expression Across GTEx Tissues (v8)", x="", y="Median TPM") + theme_bw(11))
      dev.off(); cat("  GTEx barplot done\n")
    }
  } else { cat(sprintf("GTEx API returned %d\n", status_code(resp))) }
}, error=function(e) cat(sprintf("GTEx failed: %s\n", e$message)))

# ============================================================================
# 2. HPA: protein expression via downloadable TSV
# ============================================================================
cat("\n========== HPA: GBP1 Protein Expression ==========\n")
tryCatch({
  hpa_url <- "https://www.proteinatlas.org/download/normal_tissue.tsv.zip"
  hpa_zip <- "extval/db_HPA/normal_tissue.tsv.zip"
  if (!file.exists(hpa_zip)) {
    cat("Downloading HPA normal_tissue.tsv.zip (~15 MB)...\n")
    download.file(hpa_url, hpa_zip, mode="wb")
  }
  hpa <- read.table(unz(hpa_zip, "normal_tissue.tsv"), header=TRUE, sep="\t")
  gbp1_hpa <- hpa[hpa$Gene.name == "GBP1" | hpa$Gene == "GBP1", ]
  if (nrow(gbp1_hpa) == 0) {
    gbp1_hpa <- hpa[grep("GBP1", hpa$Gene.name, ignore.case=TRUE), ]
  }
  if (nrow(gbp1_hpa) > 0) {
    write.table(gbp1_hpa, "extval/db_HPA/GBP1_HPA_tissues.txt", sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("HPA: %d tissues with GBP1 expression\n", nrow(gbp1_hpa)))
    print(gbp1_hpa[, c("Tissue","Cell.type","Level","Reliability")])
    # Barplot of nTPM
    if ("nTPM" %in% colnames(gbp1_hpa)) {
      gbp1_hpa <- gbp1_hpa[order(gbp1_hpa$nTPM, decreasing=TRUE), ]
      pdf("extval/db_HPA/GBP1_HPA_barplot.pdf", width=10, height=7)
      top20 <- head(gbp1_hpa, 20)
      top20$Tissue <- factor(top20$Tissue, levels=rev(top20$Tissue))
      print(ggplot(top20, aes(x=Tissue, y=nTPM)) +
        geom_bar(stat="identity", fill="#DE2D26") + coord_flip() +
        labs(title="GBP1 RNA Expression (Human Protein Atlas)", x="", y="nTPM") + theme_bw(11))
      dev.off(); cat("  HPA barplot done\n")
    }
  } else { cat("GBP1 not found in HPA tissue data\n") }
}, error=function(e) cat(sprintf("HPA failed: %s\n", e$message)))

# ============================================================================
# 3. GWAS Catalog: SNPs near GBP1 via REST API
# ============================================================================
cat("\n========== GWAS Catalog: GBP1 Locus ==========\n")
tryCatch({
  # Search associations by gene
  gw_url <- "https://www.ebi.ac.uk/gwas/rest/api/associations/search"
  resp <- GET(gw_url, query=list(geneName="GBP1", pageSize=50))
  if (status_code(resp) == 200) {
    gw <- fromJSON(content(resp, "text", encoding="UTF-8"))
    n_assoc <- gw$page$totalElements
    cat(sprintf("GWAS associations for GBP1: %d\n", n_assoc))

    if (n_assoc > 0) {
      assoc_list <- gw$`_embedded`$associations
      gw_df <- data.frame(
        pvalue=assoc_list$pvalue,
        trait=sapply(assoc_list$trait, function(x) paste(unique(unlist(x)), collapse="; ")),
        riskAllele=sapply(assoc_list$loci, function(x) {
          alleles <- unlist(x$strongestRiskAlleles)
          if (length(alleles) > 0) alleles[1] else ""
        }),
        stringsAsFactors=FALSE)
      gw_df <- gw_df[order(gw_df$pvalue), ]
      write.table(gw_df, "extval/db_GWAS/GBP1_GWAS_associations.txt", sep="\t", quote=FALSE, row.names=FALSE)
      cat(sprintf("Top trait: %s (P=%.2e)\n", gw_df$trait[1], gw_df$pvalue[1]))
    }
  } else {
    cat(sprintf("GWAS API returned %d\n", status_code(resp)))
  }
}, error=function(e) cat(sprintf("GWAS failed: %s\n", e$message)))

cat("\n========== Multi-DB Phase 1 Complete ==========\n")
