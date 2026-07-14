###############################################################################
# Resume from AUC evaluation — models already trained & saved to 13/
# Run after partial completion: Rscript resume_from_auc.R
###############################################################################
setwd("F:/GBP1_pipeline_2hao")
options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))

suppressMessages({
  library(limma); library(sva); library(WGCNA)
  library(ggplot2); library(ggrepel); library(pheatmap)
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
  library(ggvenn); library(caret); library(pROC)
  library(randomForestSRC); library(plsRglm); library(gbm)
  library(mboost); library(e1071); library(xgboost); library(ComplexHeatmap)
})

source("refer.ML.R")

# ---- R 4.6 cor fix ----
cor <- function(x, y = NULL, use = "everything", method = c("pearson", "kendall", "spearman"), ...) {
  if (use == "p") use <- "pairwise.complete.obs"
  stats::cor(x, y, use = use, method = match.arg(method))
}

ML_SEED <- 123
logFCfilter <- 0.585; adj.P.Val.Filter <- 0.05
WGCNA_MIN_MOD <- 60; WGCNA_MERGE_H <- 0.25
classVar <- "Type"; min.selected.var <- 2

# ---- Load data ----
cat("Loading pipeline state...\n")
expr <- readRDS("6/merged_data.rds")$expr
batch_vec <- readRDS("6/merged_data.rds")$batch
type_vec <- readRDS("6/merged_data.rds")$type

# DEG
group <- factor(type_vec, levels=c("Control","Treat"))
design <- model.matrix(~0+group); colnames(design) <- c("Control","Treat")
fit <- lmFit(expr, design); fit2 <- contrasts.fit(fit, makeContrasts(Treat-Control, levels=design))
fit2 <- eBayes(fit2)
all_diff <- topTable(fit2, adjust="fdr", number=200000)
deg <- all_diff[abs(all_diff$logFC) > logFCfilter & all_diff$adj.P.Val < adj.P.Val.Filter, ]
deg_genes <- rownames(deg)

# WGCNA
dataW <- expr[apply(expr,1,sd) > 0.5, ]
if (ncol(dataW) > 150) {
  keep <- c()
  for(b in unique(batch_vec)) for(t in c("Control","Treat")) {
    idx <- which(batch_vec==b & type_vec==t)
    if(length(idx) > 15) idx <- sample(idx, 15)
    keep <- c(keep, idx)
  }
  dataW <- dataW[, sort(keep)]
}
TypeW <- type_vec[colnames(dataW)]
datExpr0 <- t(dataW)
gsg <- goodSamplesGenes(datExpr0, verbose=0)
datExpr <- datExpr0[, gsg$goodGenes]
set.seed(123)
sft <- pickSoftThreshold(datExpr, powerVector=1:20, verbose=0)
beta <- sft$powerEstimate; if(is.na(beta)) beta <- 6
net <- blockwiseModules(datExpr, power=beta, TOMType="unsigned",
  minModuleSize=WGCNA_MIN_MOD, mergeCutHeight=WGCNA_MERGE_H,
  numericLabels=TRUE, pamRespectsDendro=FALSE, verbose=0)
trait <- data.frame(Treat=ifelse(TypeW=="Treat",1,0))
modCor <- cor(net$MEs, trait, use="p")
best_mod <- names(which.max(abs(modCor[,1])))
best_genes <- colnames(datExpr)[net$colors == as.numeric(gsub("ME","",best_mod))]
disease_genes <- intersect(deg_genes, best_genes)
cat(sprintf("Disease genes: %d\n", length(disease_genes)))

# ---- ML data (recreate) ----
ml_genes <- disease_genes
if (length(ml_genes) < 5) { ml_genes <- rownames(deg) }
if (length(ml_genes) > 200) {
  top_deg <- head(rownames(deg)[order(abs(deg$logFC), decreasing=TRUE)], 200)
  ml_genes <- intersect(ml_genes, top_deg)
}
ml_expr <- expr[ml_genes, , drop=FALSE]
ml_label <- data.frame(Type=ifelse(type_vec=="Treat", 1, 0), row.names=colnames(ml_expr))
set.seed(ML_SEED)
train_idx <- createDataPartition(ml_label$Type, p=0.7, list=FALSE)
train_set <- t(ml_expr[, train_idx, drop=FALSE]); test_set <- t(ml_expr[, -train_idx, drop=FALSE])
train_lab <- data.frame(Type=ml_label$Type[train_idx], row.names=rownames(train_set))
test_lab <- data.frame(Type=ml_label$Type[-train_idx], Cohort="Test", row.names=rownames(test_set))
train_set <- scaleData(train_set, centerFlags=TRUE, scaleFlags=TRUE)
test_set <- scaleData(test_set, cohort=rep("Test", nrow(test_set)), centerFlags=TRUE, scaleFlags=TRUE)

# ---- Load trained models ----
model <- readRDS("13/model.MLmodel.rds")
cat(sprintf("Loaded %d trained models\n", length(model)))

# ---- Step 14: AUC Evaluation ----
cat("\n========== AUC Evaluation ==========\n")
AUC_list <- list()
for (m in names(model)) {
  AUC_list[[m]] <- tryCatch(
    RunEval(fit=model[[m]], Test_set=test_set, Test_label=test_lab,
            Train_set=train_set, Train_label=train_lab, Train_name="Train",
            cohortVar="Cohort", classVar=classVar),
    error=function(e) { cat(sprintf("  %s AUC failed: %s\n", m, e$message)); NULL }
  )
}
AUC_list <- AUC_list[!sapply(AUC_list, is.null)]
AUC_mat <- do.call(rbind, AUC_list)
write.table(cbind(Method=rownames(AUC_mat), AUC_mat), "14/model.AUCmatrix.txt", sep="\t", quote=FALSE, row.names=FALSE)

avg_AUC <- apply(AUC_mat, 1, mean); avg_AUC <- sort(avg_AUC, decreasing=TRUE)
AUC_mat <- AUC_mat[names(avg_AUC), ]; avg_AUC <- round(avg_AUC, 3)
best_method <- names(avg_AUC)[1]; best_genes <- ExtractVar(model[[best_method]])
cat(sprintf("Best model: %s (AUC=%.3f), %d genes\n", best_method, avg_AUC[1], length(best_genes)))

CohortCol <- RColorBrewer::brewer.pal(n=max(3, ncol(AUC_mat)), name="Paired")[1:ncol(AUC_mat)]
names(CohortCol) <- colnames(AUC_mat)
pdf("14/AUCheatmap.pdf", width=ncol(AUC_mat)*1.2+6, height=nrow(AUC_mat)*0.4)
hm <- SimpleHeatmap(Cindex_mat=AUC_mat, avg_Cindex=avg_AUC,
  CohortCol=CohortCol, barCol="steelblue", cellwidth=1.2, cellheight=0.4,
  cluster_columns=FALSE, cluster_rows=FALSE)
draw(hm, heatmap_legend_side="right", annotation_legend_side="right")
dev.off(); cat("AUCheatmap.pdf done\n")

# ---- Step 15: ROC ----
cat("\n========== ROC Curves ==========\n")
top5 <- names(avg_AUC)[1:min(5, length(avg_AUC))]
pdf("15/ROC_top5.pdf", width=8, height=7)
col_vec <- c("#DE2D26","#3182BD","#31A354","#756BB1","#E6550D"); names(col_vec) <- top5
plot.roc(test_lab$Type, CalPredictScore(model[[top5[1]]], rbind(train_set, test_set))[rownames(test_set)],
  col=col_vec[1], lwd=2, legacy.axes=TRUE, print.auc=TRUE, print.auc.cex=1.2,
  auc.polygon=TRUE, auc.polygon.col=adjustcolor(col_vec[1], 0.2))
for (i in 2:length(top5)) {
  plot.roc(test_lab$Type, CalPredictScore(model[[top5[i]]], rbind(train_set, test_set))[rownames(test_set)],
    col=col_vec[i], lwd=2, add=TRUE, print.auc=TRUE, print.auc.cex=1.2,
    print.auc.y=1-0.08*(i-1), print.auc.col=col_vec[i])
}
legend("bottomright", legend=top5, col=col_vec[1:length(top5)], lwd=2, cex=0.8)
dev.off(); cat("ROC_top5.pdf done\n")

# ---- Step 16: SHAP ----
cat("\n========== SHAP ==========\n")
shap_model <- model[[best_method]]; shap_genes <- ExtractVar(shap_model)
if (grepl("XGBoost", best_method)) {
  shp <- predict(shap_model, newdata=as.matrix(train_set[, shap_genes, drop=FALSE]), predcontrib=TRUE)
  shap_cols <- setdiff(colnames(shp), "BIAS")
  shap_imp <- sort(colMeans(abs(shp[, shap_cols, drop=FALSE])), decreasing=TRUE)
  write.table(data.frame(Gene=names(shap_imp), SHAP=shap_imp), "16/SHAP_importance.txt", sep="\t", quote=FALSE, row.names=FALSE)
  pdf("16/SHAP_bar.pdf", width=8, height=6)
  barplot(rev(shap_imp[1:min(20,length(shap_imp))]), horiz=TRUE, las=1, col="steelblue",
    xlab="mean(|SHAP|)", main=paste0("SHAP: ", best_method))
  dev.off()
  pdf("16/SHAP_beeswarm.pdf", width=10, height=8)
  top20 <- names(shap_imp)[1:min(20,length(shap_imp))]
  shp_sub <- shp[, top20, drop=FALSE]
  shap_long <- reshape2::melt(shp_sub, varnames=c("Sample","Gene"), value.name="SHAP")
  shap_long$FeatureValue <- as.vector(train_set[rep(1:nrow(train_set), length(top20)), top20, drop=FALSE])
  print(ggplot(shap_long, aes(x=SHAP, y=Gene, color=FeatureValue)) +
    geom_jitter(size=0.5, alpha=0.5, height=0.2) +
    scale_color_gradient2(low="#3182BD", mid="grey90", high="#DE2D26") +
    labs(title=paste0("SHAP Beeswarm: ", best_method)) + theme_bw(10))
  dev.off(); cat("SHAP done\n")
} else if (grepl("Lasso|Ridge|Enet|Stepglm|RF", best_method)) {
  imp <- rep(NA, length(shap_genes)); names(imp) <- shap_genes
  pdf("16/feature_bar.pdf", width=8, height=6)
  barplot(rep(1, min(20,length(shap_genes))), horiz=TRUE, las=1, col="steelblue",
    names.arg=shap_genes[1:min(20,length(shap_genes))],
    xlab="Selected by model", main=paste0("Features: ", best_method))
  dev.off(); cat("Feature list saved\n")
}

# ---- Step 17: Gene Consensus ----
cat("\n========== Gene Consensus ==========\n")
fea_list <- lapply(model, ExtractVar)
fea_df <- do.call(rbind, lapply(names(fea_list), function(m) data.frame(Gene=fea_list[[m]], Model=m)))
write.table(fea_df, "17/model.genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
gene_freq <- sort(table(fea_df$Gene), decreasing=TRUE)
write.table(data.frame(Gene=names(gene_freq), Frequency=as.integer(gene_freq)),
  "17/gene_frequency.txt", sep="\t", quote=FALSE, row.names=FALSE)
pdf("17/gene_freq_bar.pdf", width=10, height=6)
barplot(rev(head(gene_freq, 20)), horiz=TRUE, las=1, col="#756BB1",
  xlab="Times selected", main="Gene selection frequency")
dev.off(); cat("Gene consensus done\n")

# ---- Step 18: Nomogram ----
cat("\n========== Nomogram ==========\n")
if (requireNamespace("rms", quietly=TRUE)) {
  library(rms)
  nomo_genes <- head(names(gene_freq), min(5, length(gene_freq)))
  nomo_data <- as.data.frame(train_set[, nomo_genes, drop=FALSE]); nomo_data$Type <- train_lab$Type
  dd <- datadist(nomo_data); options(datadist="dd")
  nomo_fit <- lrm(as.formula(paste("Type ~", paste(nomo_genes, collapse=" + "))), data=nomo_data, x=TRUE, y=TRUE)
  pdf("18/nomogram.pdf", width=10, height=7)
  plot(nomogram(nomo_fit, fun=function(x)1/(1+exp(-x)), funlabel="Risk", conf.int=FALSE))
  dev.off(); cat("Nomogram done\n")
}

cat("\n========== ML Pipeline Complete ==========\n")
cat(sprintf("Best: %s AUC=%.3f | %d genes\n", best_method, avg_AUC[1], length(best_genes)))
