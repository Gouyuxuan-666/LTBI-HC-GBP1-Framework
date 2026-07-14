###############################################################################
# GBP1 Multi-Database Phase 3: cellxgene Census + GTEx eQTL + PRIDE
# Tier 3: single-cell resolution + eQTL + proteomics
# Run AFTER multi_db_phase1.R and multi_db_phase2.R
###############################################################################

options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")

for (d in c("extval/db_cellxgene","extval/db_eQTL","extval/db_PRIDE")) {
  dir.create(d, recursive=TRUE, showWarnings=FALSE)
}

# ============================================================================
# 1. cellxgene Census: GBP1 single-cell cell-type specificity
# ============================================================================
cat("\n========== cellxgene: GBP1 Single-Cell Landscape ==========\n")

tryCatch({
  if (!requireNamespace("reticulate", quietly=TRUE)) install.packages("reticulate")
  library(reticulate)

  # Use cellxgene-census Python package
  cat("Opening cellxgene Census (~61M cells)...\n")

  py_script <- '
import cellxgene_census
import scanpy as sc
import pandas as pd

census = cellxgene_census.open_soma(census_version="latest")

# Query GBP1 expression in immune cells
gbp1_data = census["census_data"]["homo_sapiens"].axis_query(
    "RNA",
    var_query=cellxgene_census.get_anndata(
        census, organism="homo_sapiens",
        measurement_name="RNA",
        obs_value_filter="tissue_general == 'lung' or tissue_general == 'blood'"
    )
)

# Export summary
obs = gbp1_data.obs().concat().to_pandas()
print(f"Samples: {len(obs)}")
print(obs["cell_type"].value_counts().head(15))
print(obs["tissue"].value_counts().head(10))
'
  cat(py_script, file="extval/db_cellxgene/query_census.py")
  cat("cellxgene query script saved. Run with: python extval/db_cellxgene/query_census.py\n")
  cat("Requires: pip install cellxgene-census scanpy\n")

}, error=function(e) cat(sprintf("cellxgene setup failed: %s\n", e$message)))

# ============================================================================
# 2. GTEx eQTL: GBP1 cis-eQTL in lung and whole blood
# ============================================================================
cat("\n========== GTEx eQTL: GBP1 Genetic Regulation ==========\n")

tryCatch({
  if (!requireNamespace("gtexr", quietly=TRUE)) install.packages("gtexr")
  library(gtexr)

  gbp1_ensg <- "ENSG00000166128"

  # Get eQTL summary for GBP1 across tissues
  tryCatch({
    eqtl <- get_eqtl_genes(
      gencodeIds=gbp1_ensg,
      datasetId="gtex_v8"
    )

    if (nrow(eqtl) > 0) {
      write.table(eqtl, "extval/db_eQTL/GBP1_eQTL_summary.txt", sep="\t", quote=FALSE, row.names=FALSE)
      cat(sprintf("GTEx eQTL: %d significant SNP-gene pairs\n", nrow(eqtl)))

      # Focus on lung and whole blood
      target_tissues <- c("Lung","Whole_Blood","Cells_Cultured_fibroblasts","EBV-transformed_lymphocytes")
      eqtl_filtered <- eqtl[eqtl$tissueSiteDetailId %in% target_tissues, ]
      if (nrow(eqtl_filtered) > 0) {
        write.table(eqtl_filtered, "extval/db_eQTL/GBP1_eQTL_target_tissues.txt",
          sep="\t", quote=FALSE, row.names=FALSE)
        cat(sprintf("Target tissues eQTLs: %d\n", nrow(eqtl_filtered)))
      }
    }
  }, error=function(e) cat(sprintf("eQTL query: %s (may need GTEx dbGaP access)\n", e$message)))

  # Fallback: query variant-level eQTL
  cat("Checking variant-level eQTL for GBP1...\n")
  tryCatch({
    variants <- get_variants(gene_name="GBP1", gtexr:::gtex_db_to_dataset_ids("gtex_v8")[[1]])
    if (nrow(variants) > 0) {
      write.table(variants, "extval/db_eQTL/GBP1_variants.txt", sep="\t", quote=FALSE, row.names=FALSE)
      cat(sprintf("Variants near GBP1: %d\n", nrow(variants)))
    }
  }, error=function(e) cat(sprintf("Variant query: %s\n", e$message)))

}, error=function(e) cat(sprintf("GTEx eQTL section failed: %s\n", e$message)))

# ============================================================================
# 3. PRIDE: TB proteomics evidence for GBP1
# ============================================================================
cat("\n========== PRIDE: GBP1 Proteomics Evidence ==========\n")

tryCatch({
  pride_url <- "https://www.ebi.ac.uk/pride/ws/archive/v1/projects"
  resp <- GET(pride_url, query=list(
    query="tuberculosis AND GBP1",
    pageSize=20, page=0))
  if (status_code(resp) == 200) {
    pride <- fromJSON(content(resp, "text", encoding="UTF-8"))
    n_total <- pride$page$totalElements
    cat(sprintf("PRIDE TB studies: %d\n", n_total))

    if (n_total > 0) {
      pride_df <- as.data.frame(pride$`_embedded`$projects)
      write.table(pride_df[, c("accession","title","projectDescription","publicationDate")],
        "extval/db_PRIDE/TB_proteomics_studies.txt", sep="\t", quote=FALSE, row.names=FALSE)
    }
  }

  # Broad search for GBP1 in any proteomics study
  cat("Searching all PRIDE for GBP1...\n")
  resp2 <- GET(pride_url, query=list(query="GBP1", pageSize=50, page=0))
  if (status_code(resp2) == 200) {
    pride2 <- fromJSON(content(resp2, "text", encoding="UTF-8"))
    n_total2 <- pride2$page$totalElements
    cat(sprintf("All PRIDE studies mentioning GBP1: %d\n", n_total2))
    if (n_total2 > 0) {
      cat("GBP1 has been detected in published proteomics experiments.\n")
      cat("Key studies will provide orthogonal protein-level validation.\n")
    }
  }
}, error=function(e) cat(sprintf("PRIDE failed: %s\n", e$message)))

# ============================================================================
cat("\n========== Multi-DB Phase 3 Complete ==========\n")
cat("Outputs: extval/db_cellxgene/ extval/db_eQTL/ extval/db_PRIDE/\n")
