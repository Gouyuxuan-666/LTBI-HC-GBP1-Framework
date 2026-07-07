# ============================================================
# LTBI vs HC 跨平台合并脚本
# GPL10558 (7 datasets) + GPL6947 (3 datasets)
# 探针→基因符号 → ComBat批次校正 → 合并矩阵
# ============================================================

library(GEOquery)
library(limma)
library(sva)
library(dplyr)

# ---- Paths ----
data_dir <- "C:/Users/1/Desktop/LTBI_HC_merge/data/"
out_dir  <- "C:/Users/1/Desktop/LTBI_HC_merge/output/"
dir.create(data_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

setwd(data_dir)

# ============================================================
# STEP 1: 定义所有数据集和平台
# ============================================================
cat("===== Step 1: Dataset configuration =====\n\n")

# GPL10558 (Illumina HT-12 v4) — 7 datasets
gpl10558_datasets <- list(
  GSE37250 = list(ltbi_pat = "LTBI", hc_pat = "Control|Healthy"),
  GSE39940 = list(ltbi_pat = "LTBI|Latent", hc_pat = "Control|Healthy"),
  GSE42834 = list(ltbi_pat = "LTBI|Latent", hc_pat = "Control|Healthy"),
  GSE83456 = list(ltbi_pat = "LTBI|Latent", hc_pat = "Control|Healthy"),
  GSE39939 = list(ltbi_pat = "LTBI|Latent", hc_pat = "Control|Healthy"),
  GSE74092 = list(ltbi_pat = "LTBI|Latent", hc_pat = "Control|Healthy"),
  GSE42825 = list(ltbi_pat = "LTBI|Latent", hc_pat = "Control|Healthy")
)

# GPL6947 (Illumina HT-12 v3) — 3 datasets
gpl6947_datasets <- list(
  GSE19491 = list(ltbi_pat = "LTBI|Latent", hc_pat = "Control|Healthy"),
  GSE19439 = list(ltbi_pat = "LTBI|Latent", hc_pat = "Control|Healthy"),
  GSE19444 = list(ltbi_pat = "LTBI|Latent", hc_pat = "Control|Healthy")
)

# ============================================================
# STEP 2: 下载平台注释文件（probe → gene symbol）
# ============================================================
cat("===== Step 2: Downloading platform annotations =====\n\n")

# GPL10558
if(!file.exists("GPL10558.annot.gz")) {
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz",
    "GPL10558.annot.gz")
}
gpl10558 <- read.table("GPL10558.annot.gz", header=TRUE, sep="\t",
                        comment.char="", quote="", fill=TRUE)
cat("GPL10558 probes:", nrow(gpl10558), "\n")

# GPL6947
if(!file.exists("GPL6947.annot.gz")) {
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6947/annot/GPL6947.annot.gz",
    "GPL6947.annot.gz")
}
gpl6947 <- read.table("GPL6947.annot.gz", header=TRUE, sep="\t",
                       comment.char="", quote="", fill=TRUE)
cat("GPL6947 probes:", nrow(gpl6947), "\n")

# ============================================================
# STEP 3: 处理 GPL10558 数据集
# ============================================================
cat("\n===== Step 3: Processing GPL10558 datasets =====\n\n")

gpl10558_all_expr <- list()
gpl10558_all_group <- list()

for(gse_name in names(gpl10558_datasets)) {
  cat("---", gse_name, "---\n")

  # Download series matrix
  num <- gsub("GSE", "", gse_name)
  prefix <- substr(num, 1, nchar(num) - 3)
  url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE", prefix, "nnn/",
                gse_name, "/matrix/", gse_name, "_series_matrix.txt.gz")

  fname <- paste0(gse_name, "_series_matrix.txt.gz")
  if(!file.exists(fname)) {
    tryCatch(download.file(url, fname), error=function(e) cat("  Download failed\n"))
  }

  if(!file.exists(fname)) next

  # Read the matrix
  gse <- getGEO(filename=fname, GSEMatrix=TRUE, getGPL=FALSE)
  expr <- exprs(gse)

  # Probe → gene symbol mapping
  probes <- rownames(expr)
  # Map: try Symbol column first, fall back to Gene Symbol
  probe_to_gene <- gpl10558[match(probes, gpl10558$ID), ]
  gene_col <- if("Symbol" %in% colnames(probe_to_gene)) "Symbol" else
              if("Gene.symbol" %in% colnames(probe_to_gene)) "Gene.symbol" else NULL

  if(is.null(gene_col)) {
    cat("  WARNING: Cannot find gene symbol column\n")
    next
  }

  genes <- probe_to_gene[[gene_col]]
  # Remove probes without gene symbols
  keep <- !is.na(genes) & genes != ""
  expr <- expr[keep, ]; genes <- genes[keep]

  # Deduplicate genes by max mean expression
  expr_agg <- aggregate(expr, by=list(gene=genes), FUN=mean)
  rownames(expr_agg) <- expr_agg$gene
  expr_agg$gene <- NULL
  cat("  Genes after mapping:", nrow(expr_agg), "\n")

  # ---- Extract sample groups ----
  pd <- pData(gse)
  # Find the characteristic column describing disease status
  # Common column names: characteristics_ch1, title, source_name_ch1
  char_cols <- grep("characteristics|title|source|description",
                    colnames(pd), value=TRUE, ignore.case=TRUE)

  # Combine all text columns for pattern matching
  group_text <- apply(pd[, char_cols, drop=FALSE], 1, paste, collapse=" ")
  ltbi_idx <- grepl(gpl10558_datasets[[gse_name]]$ltbi_pat, group_text, ignore.case=TRUE)
  hc_idx   <- grepl(gpl10558_datasets[[gse_name]]$hc_pat, group_text, ignore.case=TRUE)
  # Exclude active TB (ATB) from both groups
  atb_idx  <- grepl("active|ATB|tuberculosis.*active|TB patient|pulmonary TB",
                    group_text, ignore.case=TRUE)
  ltbi_idx <- ltbi_idx & !atb_idx
  hc_idx   <- hc_idx & !atb_idx

  n_ltbi <- sum(ltbi_idx); n_hc <- sum(hc_idx)
  cat(sprintf("  LTBI: %d | HC: %d | Other: %d\n",
              n_ltbi, n_hc, ncol(expr) - n_ltbi - n_hc))

  if(n_ltbi > 0 && n_hc > 0) {
    expr_sub <- expr_agg[, c(which(ltbi_idx), which(hc_idx)), drop=FALSE]
    group <- c(rep("LTBI", n_ltbi), rep("HC", n_hc))
    colnames(expr_sub) <- paste0(gse_name, "_", colnames(expr_sub))

    gpl10558_all_expr[[gse_name]] <- expr_sub
    gpl10558_all_group[[gse_name]] <- group
  } else {
    cat("  SKIPPED: insufficient LTBI or HC samples\n")
  }
}

# ============================================================
# STEP 4: 处理 GPL6947 数据集
# ============================================================
cat("\n===== Step 4: Processing GPL6947 datasets =====\n\n")

gpl6947_all_expr <- list()
gpl6947_all_group <- list()

for(gse_name in names(gpl6947_datasets)) {
  cat("---", gse_name, "---\n")

  num <- gsub("GSE", "", gse_name)
  prefix <- substr(num, 1, nchar(num) - 3)
  url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE", prefix, "nnn/",
                gse_name, "/matrix/", gse_name, "_series_matrix.txt.gz")

  fname <- paste0(gse_name, "_series_matrix.txt.gz")
  if(!file.exists(fname)) {
    tryCatch(download.file(url, fname), error=function(e) cat("  Download failed\n"))
  }

  if(!file.exists(fname)) next

  gse <- getGEO(filename=fname, GSEMatrix=TRUE, getGPL=FALSE)
  expr <- exprs(gse)

  probes <- rownames(expr)
  probe_to_gene <- gpl6947[match(probes, gpl6947$ID), ]
  gene_col <- if("Symbol" %in% colnames(probe_to_gene)) "Symbol" else
              if("Gene.symbol" %in% colnames(probe_to_gene)) "Gene.symbol" else NULL

  if(is.null(gene_col)) { cat("  WARNING: no gene column\n"); next }

  genes <- probe_to_gene[[gene_col]]
  keep <- !is.na(genes) & genes != ""
  expr <- expr[keep, ]; genes <- genes[keep]

  expr_agg <- aggregate(expr, by=list(gene=genes), FUN=mean)
  rownames(expr_agg) <- expr_agg$gene
  expr_agg$gene <- NULL
  cat("  Genes after mapping:", nrow(expr_agg), "\n")

  pd <- pData(gse)
  char_cols <- grep("characteristics|title|source|description",
                    colnames(pd), value=TRUE, ignore.case=TRUE)
  group_text <- apply(pd[, char_cols, drop=FALSE], 1, paste, collapse=" ")

  ltbi_idx <- grepl(gpl6947_datasets[[gse_name]]$ltbi_pat, group_text, ignore.case=TRUE)
  hc_idx   <- grepl(gpl6947_datasets[[gse_name]]$hc_pat, group_text, ignore.case=TRUE)
  atb_idx  <- grepl("active|ATB|tuberculosis.*active|TB patient|pulmonary TB",
                    group_text, ignore.case=TRUE)
  ltbi_idx <- ltbi_idx & !atb_idx
  hc_idx   <- hc_idx & !atb_idx

  n_ltbi <- sum(ltbi_idx); n_hc <- sum(hc_idx)
  cat(sprintf("  LTBI: %d | HC: %d | Other: %d\n",
              n_ltbi, n_hc, ncol(expr) - n_ltbi - n_hc))

  if(n_ltbi > 0 && n_hc > 0) {
    expr_sub <- expr_agg[, c(which(ltbi_idx), which(hc_idx)), drop=FALSE]
    group <- c(rep("LTBI", n_ltbi), rep("HC", n_hc))
    colnames(expr_sub) <- paste0(gse_name, "_", colnames(expr_sub))
    gpl6947_all_expr[[gse_name]] <- expr_sub
    gpl6947_all_group[[gse_name]] <- group
  } else {
    cat("  SKIPPED: insufficient samples\n")
  }
}

# ============================================================
# STEP 5: 同平台内合并
# ============================================================

cat("\n===== Step 5: Merging within platforms =====\n\n")

# --- GPL10558 within-platform merge ---
if(length(gpl10558_all_expr) > 0) {
  # Find common genes across all GPL10558 datasets
  common_genes_10558 <- Reduce(intersect, lapply(gpl10558_all_expr, rownames))
  cat("GPL10558 common genes:", length(common_genes_10558), "\n")

  # Subset to common genes and cbind
  gpl10558_merged <- do.call(cbind, lapply(gpl10558_all_expr,
    function(x) x[common_genes_10558, , drop=FALSE]))
  gpl10558_group <- unlist(gpl10558_all_group)

  cat("GPL10558 merged:", ncol(gpl10558_merged), "samples x",
      nrow(gpl10558_merged), "genes\n")
  cat("GPL10558 groups:\n"); print(table(gpl10558_group))
}

# --- GPL6947 within-platform merge ---
if(length(gpl6947_all_expr) > 0) {
  common_genes_6947 <- Reduce(intersect, lapply(gpl6947_all_expr, rownames))
  cat("\nGPL6947 common genes:", length(common_genes_6947), "\n")

  gpl6947_merged <- do.call(cbind, lapply(gpl6947_all_expr,
    function(x) x[common_genes_6947, , drop=FALSE]))
  gpl6947_group <- unlist(gpl6947_all_group)

  cat("GPL6947 merged:", ncol(gpl6947_merged), "samples x",
      nrow(gpl6947_merged), "genes\n")
  cat("GPL6947 groups:\n"); print(table(gpl6947_group))
}

# ============================================================
# STEP 6: 跨平台合并 (gene symbol intersection)
# ============================================================
cat("\n===== Step 6: Cross-platform merge =====\n\n")

data_list <- list()
batch_labels <- c()

if(exists("gpl10558_merged") && ncol(gpl10558_merged) > 0) {
  data_list[["GPL10558"]] <- gpl10558_merged
  batch_labels <- c(batch_labels, rep("GPL10558", ncol(gpl10558_merged)))
}
if(exists("gpl6947_merged") && ncol(gpl6947_merged) > 0) {
  data_list[["GPL6947"]] <- gpl6947_merged
  batch_labels <- c(batch_labels, rep("GPL6947", ncol(gpl6947_merged)))
}

cat("Platforms found:", length(data_list), "\n")

if(length(data_list) >= 2) {
  # Find genes common to ALL platforms
  cross_genes <- Reduce(intersect, lapply(data_list, rownames))
  cat("Cross-platform common genes:", length(cross_genes), "\n")

  # Subset and combine
  cross_matrices <- lapply(data_list, function(x) x[cross_genes, , drop=FALSE])
  merged_all <- do.call(cbind, cross_matrices)

  # Combine group labels
  all_groups <- c(gpl10558_group, gpl6947_group)
  names(all_groups) <- colnames(merged_all)

} else if(length(data_list) == 1) {
  merged_all <- data_list[[1]]
  cross_genes <- rownames(merged_all)
  all_groups <- if(exists("gpl10558_group")) gpl10558_group else gpl6947_group
}

cat("Final merged matrix:", ncol(merged_all), "samples x",
    nrow(merged_all), "genes\n")
cat("Final groups:\n"); print(table(all_groups))

# ============================================================
# STEP 7: ComBat batch correction
# ============================================================
cat("\n===== Step 7: ComBat batch correction =====\n\n")

if(length(unique(batch_labels)) > 1) {
  cat("Applying ComBat for", length(unique(batch_labels)), "platforms...\n")

  # ComBat needs: expression matrix, batch vector, model matrix
  mod <- model.matrix(~ all_groups)
  merged_combat <- ComBat(dat=as.matrix(merged_all), batch=batch_labels, mod=mod)

  cat("ComBat done. Checking batch effect removal...\n")
  cat("Before ComBat - platform variance explained:\n")
  # Simple check: PCA by batch
  pca_before <- prcomp(t(merged_all), scale.=TRUE)
  cat("  PC1 var by platform:", summary(aov(pca_before$x[,1] ~ batch_labels))[[1]][1,5], "\n")
  pca_after <- prcomp(t(merged_combat), scale.=TRUE)
  cat("  After - PC1 var by platform:", summary(aov(pca_after$x[,1] ~ batch_labels))[[1]][1,5], "\n")
} else {
  cat("Only one platform — skipping ComBat\n")
  merged_combat <- merged_all
}

# ============================================================
# STEP 8: Save outputs
# ============================================================
cat("\n===== Step 8: Saving outputs =====\n\n")

write.csv(merged_combat, paste0(out_dir, "merged_expression_combat.csv"))
write.csv(all_groups, paste0(out_dir, "sample_groups.csv"))
saveRDS(list(expr=merged_combat, group=all_groups, batch=batch_labels,
             genes=cross_genes),
        paste0(out_dir, "LTBI_HC_merged.rds"))

cat("\n===== Pipeline complete =====\n")
cat("Output files:\n")
cat("  merged_expression_combat.csv  —", ncol(merged_combat), "samples\n")
cat("  sample_groups.csv              — group labels\n")
cat("  LTBI_HC_merged.rds             — full R object\n")
cat("\nGroup summary:\n")
print(table(all_groups))
