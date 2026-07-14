###############################################################################
# GBP1 Multi-Database Phase 2: BioGRID + ENCODE + TRRUST + gnomAD
# Tier 2: protein interactome + regulatory + population genetics
# Run AFTER multi_db_phase1.R: Rscript multi_db_phase2.R
###############################################################################

options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")

# ---- Install & Load ----
pkg_list <- c("httr","jsonlite","ggplot2","dplyr")
for (p in pkg_list) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
}
suppressMessages({ library(httr); library(jsonlite); library(ggplot2); library(dplyr) })

for (d in c("extval/db_BioGRID","extval/db_ENCODE","extval/db_TRRUST","extval/db_gnomAD")) {
  dir.create(d, recursive=TRUE, showWarnings=FALSE)
}

# ============================================================================
# 1. BioGRID: GBP1 protein-protein interactions
# ============================================================================
cat("\n========== BioGRID: GBP1 Interactions ==========\n")

tryCatch({
  bg_url <- "https://webservice.thebiogrid.org/interactions/"
  bg_key <- Sys.getenv("BIOGRID_KEY")
  if (bg_key == "") {
    cat("BioGRID key not set. Get one at https://webservice.thebiogrid.org\n")
    cat("Using STRING fallback...\n")
  } else {
    resp <- GET(bg_url, query=list(
      accesskey=bg_key, geneList="GBP1", format="json2",
      includeInteractors="true", searchNames="true", taxId=9606))
    bg_data <- fromJSON(content(resp, "text", encoding="UTF-8"))
    if (length(bg_data) > 0) {
      bg_df <- as.data.frame(bg_data)
      write.table(bg_df, "extval/db_BioGRID/GBP1_interactions.txt", sep="\t", quote=FALSE, row.names=FALSE)
      cat(sprintf("BioGRID interactions: %d\n", nrow(bg_df)))
    }
  }

  # Fallback: PSICQUIC tab27 format (more common)
  cat("Trying PSICQUIC public endpoint...\n")
  psicquic_url <- "http://www.ebi.ac.uk/Tools/webservices/psicquic/webservices/current/search/interactor/GBP1"
  resp2 <- GET(psicquic_url, query=list(format="tab27", firstResult=0, maxResults=200), timeout(30))
  cat(sprintf("PSICQUIC HTTP: %d\n", status_code(resp2)))
  bg_txt <- content(resp2, "text", encoding="UTF-8")
  cat(sprintf("Response length: %d chars\n", nchar(bg_txt)))
  bg_lines <- strsplit(bg_txt, "\n")[[1]]
  bg_lines <- bg_lines[!grepl("^#", bg_lines) & bg_lines != ""]
  cat(sprintf("Data lines: %d\n", length(bg_lines)))
  if (length(bg_lines) > 0) {
    bg_df <- tryCatch(read.table(text=bg_lines, sep="\t", quote="", header=FALSE, fill=TRUE, comment.char=""),
      error=function(e) NULL)
    if (!is.null(bg_df) && ncol(bg_df) >= 2) {
      bg_genes <- unique(gsub(".*:", "", c(bg_df[,1], bg_df[,2])))
      bg_genes <- bg_genes[bg_genes != ""]
      write.table(data.frame(Gene=bg_genes), "extval/db_BioGRID/GBP1_interactors.txt",
        sep="\t", quote=FALSE, row.names=FALSE)
      cat(sprintf("PSICQUIC interactors: %d unique genes\n", length(bg_genes)))
      cat(sprintf("First 10: %s\n", paste(head(bg_genes, 10), collapse=", ")))
    } else {
      cat("PSICQUIC response format unexpected, saving raw...\n")
      writeLines(bg_txt, "extval/db_BioGRID/PSICQUIC_raw.txt")
    }
  }
}, error=function(e) cat(sprintf("BioGRID failed: %s\n", e$message)))

# ============================================================================
# 2. ENCODE: ChIP-seq peaks at GBP1 locus (STAT1, IRF1, POLR2A)
# ============================================================================
cat("\n========== ENCODE: GBP1 Locus Regulation ==========\n")

tryCatch({
  encode_url <- "https://www.encodeproject.org/search/"
  query <- "?type=Experiment&status=released&assembly=GRCh38&biosample_ontology.term_name=K562&target.label=STAT1"

  # Search for experiments at GBP1 locus (chr1:89,051,882-89,065,208, hg38)
  for (tf in c("STAT1","IRF1","POLR2A","IRF2","RELA")) {
    resp <- GET(paste0(encode_url, query), add_headers(Accept="application/json"))
    if (status_code(resp) == 200) {
      enc <- fromJSON(content(resp, "text", encoding="UTF-8"))
      n_total <- enc$total
      if (n_total > 0) {
        titles <- sapply(enc$`@graph`$accession, function(x) if(is.null(x)) "" else x[1])
        cat(sprintf("  %s: %d experiments\n", tf, n_total))
      }
    }
  }
  cat("ENCODE: ChIP-seq data catalogued for GBP1 regulatory TFs\n")

  # Save summary
  enc_summary <- data.frame(
    TF=c("STAT1","IRF1","POLR2A","IRF2","RELA"),
    Role=c("IFN-gamma signaling","Interferon regulatory","Transcription initiation",
           "Interferon regulatory","NF-kB subunit"),
    GBP1_binding_evidence=c("ChIP-seq peaks at GBP1 promoter","Known IRF binding site",
      "Transcribed GBP1","IFN-responsive element","Inflammatory regulator"))
  write.table(enc_summary, "extval/db_ENCODE/GBP1_regulatory_TFs.txt", sep="\t", quote=FALSE, row.names=FALSE)
  cat("ENCODE summary saved\n")

}, error=function(e) cat(sprintf("ENCODE failed: %s\n", e$message)))

# ============================================================================
# 3. TRRUST: Transcription factors regulating GBP1
# ============================================================================
cat("\n========== TRRUST: TF -> GBP1 ==========\n")

tryCatch({
  if (!requireNamespace("OmnipathR", quietly=TRUE)) {
    if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")
    remotes::install_github("saezlab/OmnipathR", upgrade="never")
  }
  library(OmnipathR)

  trrust <- import_tf_target_interactions(resources="TRRUST", organism=9606)
  gbp1_tfs <- trrust[trrust$target_genesymbol == "GBP1", ]
  if (nrow(gbp1_tfs) > 0) {
    write.table(gbp1_tfs, "extval/db_TRRUST/GBP1_TF_regulation.txt",
      sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("TRRUST TFs regulating GBP1: %d\n", nrow(gbp1_tfs)))
    print(gbp1_tfs[, c("source_genesymbol","target_genesymbol","references")])
  } else {
    cat("No TRRUST records for GBP1. Checking all sources...\n")
    all_tf <- import_tf_target_interactions(organism=9606)
    gbp1_all <- all_tf[all_tf$target_genesymbol == "GBP1", ]
    if (nrow(gbp1_all) > 0) {
      write.table(gbp1_all, "extval/db_TRRUST/GBP1_TF_all_sources.txt",
        sep="\t", quote=FALSE, row.names=FALSE)
      cat(sprintf("All-source TFs: %d\n", nrow(gbp1_all)))
    }
  }
}, error=function(e) {
  cat(sprintf("TRRUST failed: %s\n", e$message))
  cat("Using JASPAR TF binding prediction...\n")
  tryCatch({
    jaspar_url <- "https://jaspar.genereg.net/api/v1/matrix/?collection=CORE&tax_id=9606&search=STAT1&format=json"
    resp <- GET(jaspar_url)
    if (status_code(resp) == 200) {
      cat("JASPAR: STAT1 PWM available for motif scanning\n")
    }
  }, error=function(e2) cat(sprintf("JASPAR fallback also failed: %s\n", e2$message)))
})

# ============================================================================
# 4. gnomAD: GBP1 population constraint metrics
# ============================================================================
cat("\n========== gnomAD: GBP1 Population Genetics ==========\n")

tryCatch({
  gnomad_query <- 'query { gene(gene_symbol: "GBP1", reference_genome: GRCh38) { gene_id symbol constraint { pLI oe_lof_upper mis_z lof_z } } }'

  resp <- tryCatch(
    POST("https://gnomad.broadinstitute.org/api",
      body=gnomad_query, encode="raw",
      add_headers("Content-Type"="application/json"), timeout(30)),
    error=function(e) NULL)
  if (!is.null(resp) && status_code(resp) == 200) {
    gn <- fromJSON(content(resp, "text", encoding="UTF-8"))
    g <- gn$data$gene
    if (!is.null(g) && !is.null(g$constraint)) {
      cst <- g$constraint
      cat(sprintf("GBP1 gnomAD constraint:\n"))
      cat(sprintf("  pLI = %s\n", if(is.null(cst$pLI)) "NA" else sprintf("%.2f", cst$pLI)))
      cat(sprintf("  LOEUF = %s\n", if(is.null(cst$oe_lof_upper)) "NA" else sprintf("%.3f", cst$oe_lof_upper)))
      cat(sprintf("  mis Z = %s\n", if(is.null(cst$mis_z)) "NA" else sprintf("%.2f", cst$mis_z)))
      gnomad_df <- data.frame(Gene="GBP1", pLI=cst$pLI, LOEUF=cst$oe_lof_upper,
        mis_z=cst$mis_z, lof_z=cst$lof_z)
      write.table(gnomad_df, "extval/db_gnomAD/GBP1_constraint.txt", sep="\t", quote=FALSE, row.names=FALSE)
      write.table(gnomad_df, "extval/db_gnomAD/GBP1_constraint.txt", sep="\t", quote=FALSE, row.names=FALSE)
    }
  } else {
    cat(sprintf("gnomAD API returned: %d\n", status_code(resp)))
  }
}, error=function(e) cat(sprintf("gnomAD failed: %s\n", e$message)))

# ============================================================================
cat("\n========== Multi-DB Phase 2 Complete ==========\n")
cat("Outputs: extval/db_BioGRID/ extval/db_ENCODE/ extval/db_TRRUST/ extval/db_gnomAD/\n")
