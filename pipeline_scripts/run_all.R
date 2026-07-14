###############################################################################
# GBP1-TB Full Pipeline for 2号机 (200GB RAM)
# Data: 6 GEO datasets, 13836 genes x 796 samples
# Steps: SVA merge → QC → DEG → WGCNA → Disease genes → ML → SHAP
#        TCGA immune infiltration (CIBERSORT + GBP1 correlation)
# Usage: Rscript run_all.R
###############################################################################

# ---- CONFIG (IRON RULES from skill) ----
# Set working directory to script location
args <- commandArgs(trailingOnly=FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0 && script_path != "") {
  setwd(dirname(script_path))
}
cat(sprintf("Working dir: %s\n", getwd()))

# Create output directories
for (d in c("6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25")) {
  dir.create(d, showWarnings=FALSE)
}

logFCfilter <- 0.585
adj.P.Val.Filter <- 0.05
WGCNA_MIN_MOD <- 60
WGCNA_MERGE_H <- 0.25
ML_SEED <- 123

# ---- Config ----
options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))

# ---- Install packages if needed ----
need_pkgs <- c("limma","sva","WGCNA","ggplot2","ggrepel","pheatmap",
               "caret","pROC","randomForest","glmnet","clusterProfiler",
               "org.Hs.eg.db","enrichplot","ggvenn","UpSetR",
               "reshape2","corrplot","ggpubr",
               "randomForestSRC","plsRglm","gbm","mboost","e1071","xgboost","ComplexHeatmap")
for (pkg in need_pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg)
}
# Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("illuminaHumanv4.db","impute","preprocessCore","sva","limma","WGCNA",
  "clusterProfiler","org.Hs.eg.db","enrichplot","ggvenn"), update=FALSE, ask=FALSE)

suppressMessages({
  library(limma); library(sva); library(WGCNA)
  library(ggplot2); library(ggrepel); library(pheatmap)
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
  library(ggvenn)
})

cat("Sequential mode (stable)\n")

# R 4.6 compatibility: translate WGCNA's use codes and swallow extra args
cor <- function(x, y = NULL, use = "everything", method = c("pearson", "kendall", "spearman"), ...) {
  if (use == "p") use <- "pairwise.complete.obs"
  if (use == "na.or.complete") use <- "na.or.complete"
  stats::cor(x, y, use = use, method = match.arg(method))
}

# ==================== STEP 6: MERGE ====================
cat("\n========== STEP 6: SVA Merge ==========\n")
gse_all <- c("GSE83456","GSE34608","GSE19491","GSE37250","GSE62525","GSE42830")

expr_list <- list(); batch_vec <- c(); type_vec <- c()
for (gse in gse_all) {
  f <- file.path("5", paste0(gse, ".normalize.txt"))
  if (!file.exists(f)) { cat(sprintf("SKIP %s\n", gse)); next }
  rt <- as.data.frame(data.table::fread(f, header=TRUE, sep="\t", check.names=FALSE))
  rownames(rt) <- rt[[1]]; rt[[1]] <- NULL
  rt <- as.matrix(rt); mode(rt) <- "numeric"
  rt <- rt[rowSums(is.na(rt)) == 0, , drop=FALSE]
  cnames <- colnames(rt)

  # Classify samples
  types <- rep(NA, length(cnames))
  types[grepl("_Control$|_control$", cnames)] <- "Control"
  types[grepl("_Treat$|_treat$|_TB$", cnames)] <- "Treat"

  # For GSM-only names (new datasets), use s1/s2
  if (any(is.na(types))) {
    s1f <- file.path("5", gse, "s1.txt"); s2f <- file.path("5", gse, "s2.txt")
    if (file.exists(s1f) && file.exists(s2f)) {
      s1c <- colnames(read.table(s1f, header=TRUE, sep="\t", check.names=FALSE, row.names=1))
      s2c <- colnames(read.table(s2f, header=TRUE, sep="\t", check.names=FALSE, row.names=1))
      types[is.na(types) & cnames %in% s1c] <- "Control"
      types[is.na(types) & cnames %in% s2c] <- "Treat"
    }
  }

  keep <- !is.na(types)
  if (sum(keep) == 0) { cat(sprintf("  %s: no typed samples, skip\n", gse)); next }
  rt <- rt[, keep, drop=FALSE]; types <- types[keep]
  expr_list[[gse]] <- rt
  batch_vec <- c(batch_vec, rep(gse, ncol(rt)))
  type_vec <- c(type_vec, types)
  cat(sprintf("  %s: %d genes x %d samples (C:%d T:%d)\n", gse, nrow(rt), ncol(rt),
              sum(types=="Control"), sum(types=="Treat")))
}

common <- Reduce(intersect, lapply(expr_list, rownames))
cat(sprintf("Common genes: %d\n", length(common)))

merged <- do.call(cbind, lapply(expr_list, function(x) x[common, , drop=FALSE]))
write.table(merged, "6/merge.preNorm.txt", sep="\t", quote=FALSE)
cat(sprintf("Pre-norm: %d x %d\n", nrow(merged), ncol(merged)))

merged_norm <- ComBat(merged, batch=batch_vec, par.prior=TRUE)
merged_norm <- normalizeBetweenArrays(merged_norm)
write.table(merged_norm, "6/merge.normalize.txt", sep="\t", quote=FALSE)
cat(sprintf("Post-norm: %d x %d\n", nrow(merged_norm), ncol(merged_norm)))
cat(sprintf("Batch: %s\n", paste(unique(batch_vec), collapse=", ")))

saveRDS(list(expr=merged_norm, batch=batch_vec, type=type_vec), "6/merged_data.rds")

# ==================== STEP 7-8: QC ====================
cat("\n========== QC ==========\n")
expr <- merged_norm
rt500 <- expr[order(-apply(expr,1,sd)), ][1:min(500,nrow(expr)), ]
data_t <- t(rt500)

df_b <- data.frame(Project=batch_vec, Type=type_vec, Median=apply(data_t, 1, median))
pdf("7/boxplot_postNorm.pdf", width=12, height=5)
print(ggplot(df_b, aes(x=Project, y=Median, fill=Type)) +
  geom_boxplot(outlier.size=0.3) + scale_fill_manual(values=c("Control"="#3182BD","Treat"="#DE2D26")) +
  labs(title="After ComBat Batch Correction") + theme_bw(10) +
  theme(axis.text.x=element_text(angle=30, hjust=1)))
dev.off()

pca <- prcomp(data_t, scale.=TRUE); pred <- predict(pca)
var1 <- round(summary(pca)$importance[2,1]*100,1); var2 <- round(summary(pca)$importance[2,2]*100,1)
df_p <- data.frame(PC1=pred[,1], PC2=pred[,2], Cohort=batch_vec, Type=type_vec)
pdf("8/PCA_postNorm.pdf", width=9, height=6)
print(ggplot(df_p, aes(x=PC1, y=PC2, color=Cohort, shape=Type)) +
  geom_point(size=1.0, alpha=0.4) + scale_shape_manual(values=c("Control"=16,"Treat"=17)) +
  scale_color_brewer(palette="Set1") +
  labs(x=paste0("PC1 (",var1,"%)"), y=paste0("PC2 (",var2,"%)"), title="After ComBat") +
  theme_bw(10))
dev.off()
cat("QC done\n")

# ==================== STEP 9: DEG ====================
cat("\n========== DEG ==========\n")
group <- factor(type_vec, levels=c("Control","Treat"))
design <- model.matrix(~0+group); colnames(design) <- c("Control","Treat")
fit <- lmFit(expr, design); fit2 <- contrasts.fit(fit, makeContrasts(Treat-Control, levels=design))
fit2 <- eBayes(fit2)
all_diff <- topTable(fit2, adjust="fdr", number=200000)

deg <- all_diff[abs(all_diff$logFC) > logFCfilter & all_diff$adj.P.Val < adj.P.Val.Filter, ]
write.table(deg, "9/diff.txt", sep="\t", quote=FALSE)
cat(sprintf("DEGs: %d (Up:%d Down:%d)\n", nrow(deg), sum(deg$logFC>0), sum(deg$logFC<0)))
top10 <- head(deg[order(-abs(deg$logFC)), ], 10)
cat(sprintf("Top10: %s\n", paste(rownames(top10), collapse=", ")))

# Volcano
all_diff$Sig <- with(all_diff,
  ifelse((adj.P.Val<adj.P.Val.Filter)&(abs(logFC)>logFCfilter),
    ifelse(logFC>logFCfilter,"Up","Down"),"Not"))
pdf("9/vol.pdf", width=6, height=5)
all_diff$ID <- rownames(all_diff)
top_lbl <- head(all_diff[order(all_diff$adj.P.Val), ], 10)
print(ggplot(all_diff, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(color=Sig), size=0.3, alpha=0.5) +
  scale_color_manual(values=c("Down"="#3182BD","Not"="grey85","Up"="#DE2D26")) +
  geom_label_repel(data=top_lbl, aes(label=ID), size=2.5, box.padding=0.3, max.overlaps=30) +
  labs(title=paste0("DEGs: ",sum(all_diff$Sig=="Up")," up, ",sum(all_diff$Sig=="Down")," down")) + theme_bw(11))
dev.off()

de <- expr[rownames(expr) %in% rownames(deg), ]
write.table(de, "9/diffGeneExp.txt", sep="\t", quote=FALSE)
write.table(all_diff, "9/all.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("DEG done\n")

# ==================== STEP 10: WGCNA ====================
cat("\n========== WGCNA ==========\n")
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
cat(sprintf("WGCNA input: %d genes x %d samples\n", nrow(dataW), ncol(dataW)))

datExpr0 <- t(dataW)
gsg <- goodSamplesGenes(datExpr0, verbose=0)
datExpr <- datExpr0[, gsg$goodGenes]

sft <- pickSoftThreshold(datExpr, powerVector=1:20, verbose=0)
pdf("10/3_scale_independence.pdf", width=9, height=5)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab="Soft Threshold", ylab="Scale Free Topology Fit", type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=1:20, col="red")
abline(h=0.85, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold", ylab="Mean Connectivity", type="n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=1:20, col="red")
dev.off()

beta <- sft$powerEstimate; if(is.na(beta)) beta <- 6
cat(sprintf("Soft threshold power: %d\n", beta))

net <- blockwiseModules(datExpr, power=beta, TOMType="unsigned",
  minModuleSize=WGCNA_MIN_MOD, mergeCutHeight=WGCNA_MERGE_H,
  numericLabels=TRUE, pamRespectsDendro=FALSE, verbose=1)

trait <- data.frame(Treat=ifelse(TypeW=="Treat",1,0))
modCor <- cor(net$MEs, trait, use="p")
modP <- corPvalueStudent(modCor, nrow(datExpr))

pdf("10/8_Module_trait.pdf", width=8, height=6)
labeledHeatmap(Matrix=modCor, xLabels="Treat", yLabels=names(net$MEs),
  colorLabels=FALSE, colors=blueWhiteRed(50),
  textMatrix=paste(signif(modCor,2),"\n(",signif(modP,1),")",sep=""),
  setStdMargins=FALSE, cex.text=0.7, zlim=c(-1,1), main="Module-Trait")
dev.off()

best_mod <- names(which.max(abs(modCor[,1])))
best_genes <- colnames(datExpr)[net$colors == as.numeric(gsub("ME","",best_mod))]
write.table(best_genes, file.path("10", paste0("module_", tolower(best_mod), ".txt")),
  row.names=FALSE, col.names=FALSE, quote=FALSE)
cat(sprintf("Best module %s: %d genes\n", best_mod, length(best_genes)))
for(m in unique(net$colors)) {
  write.table(colnames(datExpr)[net$colors==m],
    file.path("10", paste0("module_",m,".txt")), row.names=FALSE, col.names=FALSE, quote=FALSE)
}

# ==================== STEP 11: Disease Genes (DEG ∩ WGCNA) ====================
cat("\n========== Disease Genes ==========\n")
deg_genes <- rownames(deg)
wgcna_genes <- best_genes
disease_genes <- intersect(deg_genes, wgcna_genes)
write.table(disease_genes, "11/Disease.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
cat(sprintf("DEG ∩ WGCNA(%s): %d genes\n", best_mod, length(disease_genes)))

# ==================== STEP 12: ML Data Preparation ====================
cat("\n========== ML Data Prep ==========\n")
source("refer.ML.R")
suppressMessages({
  library(caret); library(randomForestSRC); library(plsRglm)
  library(gbm); library(mboost); library(e1071); library(xgboost)
  library(ComplexHeatmap); library(pROC)
})

cat(sprintf("Disease genes (DEG ∩ WGCNA): %d\n", length(disease_genes)))
if (length(disease_genes) < 5) {
  cat("WARNING: Too few disease genes! Using DEGs only.\n")
  ml_genes <- rownames(deg)
} else {
  ml_genes <- disease_genes
}
if (length(ml_genes) > 200) {
  top_deg <- head(rownames(deg)[order(abs(deg$logFC), decreasing=TRUE)], 200)
  ml_genes <- intersect(ml_genes, top_deg)
  cat(sprintf("Trimmed to top 200 by |logFC|: %d genes\n", length(ml_genes)))
}
cat(sprintf("ML features: %d genes\n", length(ml_genes)))

ml_expr <- expr[ml_genes, , drop=FALSE]
ml_label <- data.frame(Type=ifelse(type_vec=="Treat", 1, 0), row.names=colnames(ml_expr))

set.seed(ML_SEED)
train_idx <- createDataPartition(ml_label$Type, p=0.7, list=FALSE)
train_set <- t(ml_expr[, train_idx, drop=FALSE])
test_set  <- t(ml_expr[, -train_idx, drop=FALSE])
train_lab <- data.frame(Type=ml_label$Type[train_idx], row.names=rownames(train_set))
test_lab  <- data.frame(Type=ml_label$Type[-train_idx], Cohort="Test", row.names=rownames(test_set))

train_set <- scaleData(train_set, centerFlags=TRUE, scaleFlags=TRUE)
test_set  <- scaleData(test_set, cohort=rep("Test", nrow(test_set)), centerFlags=TRUE, scaleFlags=TRUE)

write.table(cbind(train_set, Type=train_lab$Type), "12/data.train.txt", sep="\t", quote=FALSE)
write.table(cbind(test_set, Type=test_lab$Type), "12/data.test.txt", sep="\t", quote=FALSE)
cat(sprintf("Train: %d, Test: %d\n", nrow(train_set), nrow(test_set)))

# ==================== STEP 13: ML Model Training ====================
cat("\n========== ML Training (113 models) ==========\n")
method_file <- "refer.methodLists.txt"
if (!file.exists(method_file)) {
  methods_all <- c("Lasso","Ridge","Enet[alpha=0.1]","Enet[alpha=0.2]","Enet[alpha=0.3]",
    "Enet[alpha=0.4]","Enet[alpha=0.5]","Enet[alpha=0.6]","Enet[alpha=0.7]","Enet[alpha=0.8]",
    "Enet[alpha=0.9]","Stepglm[direction=both]","Stepglm[direction=backward]",
    "Stepglm[direction=forward]","SVM","LDA","glmBoost","plsRglm","RF","GBM","XGBoost",
    "NaiveBayes","Lasso+RF","Lasso+GBM","Lasso+XGBoost","Lasso+glmBoost",
    "Ridge+RF","Ridge+GBM","Ridge+XGBoost","Ridge+glmBoost",
    "Enet[alpha=0.1]+RF","Enet[alpha=0.5]+RF","Enet[alpha=0.9]+RF",
    "Enet[alpha=0.1]+GBM","Enet[alpha=0.5]+GBM","Enet[alpha=0.9]+GBM",
    "Enet[alpha=0.1]+XGBoost","Enet[alpha=0.5]+XGBoost","Enet[alpha=0.9]+XGBoost",
    "Enet[alpha=0.1]+glmBoost","Enet[alpha=0.5]+glmBoost","Enet[alpha=0.9]+glmBoost",
    "Stepglm[direction=both]+RF","Stepglm[direction=both]+GBM","Stepglm[direction=both]+XGBoost",
    "Stepglm[direction=both]+glmBoost","Stepglm[direction=backward]+RF",
    "Stepglm[direction=backward]+GBM","Stepglm[direction=backward]+XGBoost",
    "Stepglm[direction=backward]+glmBoost","Stepglm[direction=forward]+RF",
    "Stepglm[direction=forward]+GBM","Stepglm[direction=forward]+XGBoost",
    "Stepglm[direction=forward]+glmBoost","SVM+RF","SVM+GBM","SVM+XGBoost",
    "SVM+glmBoost","LDA+RF","LDA+GBM","LDA+XGBoost","LDA+glmBoost",
    "glmBoost+RF","glmBoost+GBM","glmBoost+XGBoost","plsRglm+RF","plsRglm+GBM",
    "plsRglm+XGBoost","plsRglm+glmBoost","NaiveBayes+RF","NaiveBayes+GBM",
    "NaiveBayes+XGBoost","NaiveBayes+glmBoost")
  write.table(data.frame(Model=methods_all), method_file, sep="\t", quote=FALSE, row.names=FALSE)
}
methodRT <- read.table(method_file, header=TRUE, sep="\t")
methods <- methodRT$Model
methods <- gsub("-| ", "", methods)

classVar <- "Type"
min.selected.var <- 2
Variable <- colnames(train_set)

preTrain.method <- strsplit(methods, "\\+")
preTrain.method <- lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method <- unique(unlist(preTrain.method))

preTrain.var <- list()
set.seed(ML_SEED)
cat(sprintf("  Pre-training %d methods...\n", length(preTrain.method)))
preTrain.var <- lapply(preTrain.method, function(method) {
  tryCatch(
    RunML(method=method, Train_set=train_set, Train_label=train_lab, mode="Variable", classVar=classVar),
    error=function(e) { cat(sprintf("  %s var sel failed: %s\n", method, e$message)); colnames(train_set) }
  )
})
names(preTrain.var) <- preTrain.method
preTrain.var[["simple"]] <- colnames(train_set)

model <- list()
set.seed(ML_SEED)
train_set_bk <- train_set
cat(sprintf("  Training %d models...\n", length(methods)))
model <- list()
set.seed(ML_SEED)
train_set_bk <- train_set
for (i in seq_along(methods)) {
  method <- methods[i]
  cat(sprintf("  [%d/%d] %s\n", i, length(methods), method))
  method_parts <- strsplit(method, "\\+")[[1]]
  if (length(method_parts) == 1) method_parts <- c("simple", method_parts)
  Variable <- preTrain.var[[method_parts[1]]]
  train_set <- train_set_bk[, Variable, drop=FALSE]
  fit <- tryCatch(
    RunML(method=method_parts[2], Train_set=train_set, Train_label=train_lab, mode="Model", classVar=classVar),
    error=function(e) { cat(sprintf("  %s failed: %s\n", method, e$message)); NULL }
  )
  if (!is.null(fit) && length(ExtractVar(fit)) <= min.selected.var) fit <- NULL
  if (!is.null(fit)) model[[method]] <- fit
}
train_set <- train_set_bk; rm(train_set_bk)
saveRDS(model, "13/model.MLmodel.rds")
cat(sprintf("Trained: %d models\n", length(model)))

# ==================== STEP 14: AUC Evaluation ====================
cat("\n========== AUC Evaluation ==========\n")
AUC_list <- list()
cat(sprintf("  Evaluating %d models...\n", length(model)))
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

avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing=TRUE)
AUC_mat <- AUC_mat[names(avg_AUC), ]
avg_AUC <- round(avg_AUC, 3)
best_method <- names(avg_AUC)[1]
best_genes <- ExtractVar(model[[best_method]])
cat(sprintf("Best model: %s (AUC=%.3f), %d genes\n", best_method, avg_AUC[1], length(best_genes)))
cat(sprintf("Top 10:\n%s\n", paste(sprintf("  %s: %.3f", names(avg_AUC)[1:min(10,length(avg_AUC))],
  avg_AUC[1:min(10,length(avg_AUC))]), collapse="\n")))

CohortCol <- RColorBrewer::brewer.pal(n=max(3, ncol(AUC_mat)), name="Paired")[1:ncol(AUC_mat)]
names(CohortCol) <- colnames(AUC_mat)
pdf("14/AUCheatmap.pdf", width=ncol(AUC_mat)*1.2+6, height=nrow(AUC_mat)*0.4)
hm <- SimpleHeatmap(Cindex_mat=AUC_mat, avg_Cindex=avg_AUC,
  CohortCol=CohortCol, barCol="steelblue",
  cellwidth=1.2, cellheight=0.4, cluster_columns=FALSE, cluster_rows=FALSE)
draw(hm, heatmap_legend_side="right", annotation_legend_side="right")
dev.off(); cat("  AUCheatmap.pdf done\n")

# ==================== STEP 15: ROC Curves ====================
cat("\n========== ROC Curves ==========\n")
top5 <- names(avg_AUC)[1:min(5, length(avg_AUC))]
pdf("15/ROC_top5.pdf", width=8, height=7)
col_vec <- c("#DE2D26","#3182BD","#31A354","#756BB1","#E6550D")
names(col_vec) <- top5
plot.roc(test_lab$Type, CalPredictScore(model[[top5[1]]], rbind(train_set, test_set))[rownames(test_set)],
  col=col_vec[1], lwd=2, legacy.axes=TRUE, print.auc=TRUE, print.auc.cex=1.2,
  auc.polygon=TRUE, auc.polygon.col=adjustcolor(col_vec[1], 0.2))
for (i in 2:length(top5)) {
  plot.roc(test_lab$Type, CalPredictScore(model[[top5[i]]], rbind(train_set, test_set))[rownames(test_set)],
    col=col_vec[i], lwd=2, add=TRUE, print.auc=TRUE, print.auc.cex=1.2,
    print.auc.y=1-0.08*(i-1), print.auc.col=col_vec[i])
}
legend("bottomright", legend=top5, col=col_vec[1:length(top5)], lwd=2, cex=0.8)
dev.off(); cat("  ROC_top5.pdf done\n")

# ==================== STEP 16: SHAP Analysis ====================
cat("\n========== SHAP ==========\n")
shap_model <- model[[best_method]]
shap_genes <- ExtractVar(shap_model)
cat(sprintf("SHAP on %s: %d genes\n", best_method, length(shap_genes)))

if (grepl("XGBoost", best_method)) {
  cat("Using xgboost built-in SHAP contributions...\n")
  shp <- predict(shap_model, newdata=as.matrix(train_set[, shap_genes, drop=FALSE]), predcontrib=TRUE)
  shap_cols <- setdiff(colnames(shp), "BIAS")
  shap_imp <- colMeans(abs(shp[, shap_cols, drop=FALSE]))
  shap_imp <- sort(shap_imp, decreasing=TRUE)
  write.table(data.frame(Gene=names(shap_imp), SHAP=shap_imp),
    "16/SHAP_importance.txt", sep="\t", quote=FALSE, row.names=FALSE)

  pdf("16/SHAP_bar.pdf", width=8, height=6)
  top20 <- names(shap_imp)[1:min(20, length(shap_imp))]
  barplot(rev(shap_imp[top20]), horiz=TRUE, las=1, col="steelblue",
    xlab="mean(|SHAP|)", main=paste0("SHAP Feature Importance: ", best_method))
  dev.off(); cat("  SHAP_bar.pdf done\n")

  pdf("16/SHAP_beeswarm.pdf", width=10, height=8)
  shp_sub <- shp[, top20, drop=FALSE]
  shap_long <- reshape2::melt(shp_sub, varnames=c("Sample","Gene"), value.name="SHAP")
  shap_long$FeatureValue <- as.vector(train_set[rep(1:nrow(train_set), length(top20)), top20, drop=FALSE])
  print(ggplot(shap_long, aes(x=SHAP, y=Gene, color=FeatureValue)) +
    geom_jitter(size=0.5, alpha=0.5, height=0.2) +
    scale_color_gradient2(low="#3182BD", mid="grey90", high="#DE2D26") +
    labs(title=paste0("SHAP Beeswarm: ", best_method)) + theme_bw(10))
  dev.off(); cat("  SHAP_beeswarm.pdf done\n")

} else if (grepl("Lasso|Ridge|Enet|Stepglm", best_method)) {
  coefs <- coef(shap_model)[, 1]
  coefs <- coefs[coefs != 0]
  coefs <- coefs[names(coefs) != "(Intercept)"]
  coefs <- sort(abs(coefs), decreasing=TRUE)
  write.table(data.frame(Gene=names(coefs), AbsCoef=coefs),
    "16/model_coefficients.txt", sep="\t", quote=FALSE, row.names=FALSE)
  pdf("16/coef_bar.pdf", width=8, height=6)
  barplot(rev(coefs[1:min(20,length(coefs))]), horiz=TRUE, las=1, col="steelblue",
    xlab="|Coefficient|", main=paste0("Feature Coefficients: ", best_method))
  dev.off(); cat("  coef_bar.pdf done\n")
} else {
  cat(sprintf("Best model genes (%s): %s\n", best_method, paste(shap_genes, collapse=", ")))
}

# ==================== STEP 17: Gene Importance Consensus ====================
cat("\n========== Gene Consensus ==========\n")
fea_list <- list()
for (method in names(model)) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}
fea_df <- do.call(rbind, lapply(names(fea_list), function(m) {
  data.frame(Gene=fea_list[[m]], Model=m, stringsAsFactors=FALSE)
}))
write.table(fea_df, "17/model.genes.txt", sep="\t", quote=FALSE, row.names=FALSE)

gene_freq <- sort(table(fea_df$Gene), decreasing=TRUE)
cat(sprintf("Consensus genes (selected in >=3 models):\n"))
consensus <- names(gene_freq[gene_freq >= 3])
cat(sprintf("  %s\n", paste(consensus, collapse=", ")))
write.table(data.frame(Gene=names(gene_freq), Frequency=as.integer(gene_freq)),
  "17/gene_frequency.txt", sep="\t", quote=FALSE, row.names=FALSE)

pdf("17/gene_freq_bar.pdf", width=10, height=6)
top_freq <- head(gene_freq, 20)
barplot(rev(top_freq), horiz=TRUE, las=1, col="#756BB1",
  xlab="Times selected", main="Gene selection frequency across all ML models")
dev.off(); cat("  gene_freq_bar.pdf done\n")

# ==================== STEP 18: Nomogram ====================
cat("\n========== Nomogram ==========\n")
if (requireNamespace("rms", quietly=TRUE)) {
  library(rms)
  nomo_genes <- head(names(gene_freq), min(5, length(gene_freq)))
  nomo_data <- as.data.frame(train_set[, nomo_genes, drop=FALSE])
  nomo_data$Type <- train_lab$Type
  dd <- datadist(nomo_data); options(datadist="dd")
  f <- as.formula(paste("Type ~", paste(nomo_genes, collapse=" + ")))
  nomo_fit <- lrm(f, data=nomo_data, x=TRUE, y=TRUE)
  pdf("18/nomogram.pdf", width=10, height=7)
  plot(nomogram(nomo_fit, fun=function(x) 1/(1+exp(-x)), funlabel="Risk",
    conf.int=FALSE, abbrev=TRUE))
  dev.off(); cat("  nomogram.pdf done\n")

  pdf("18/calibration.pdf", width=6, height=6)
  cal <- calibrate(nomo_fit, method="boot", B=200)
  plot(cal, xlab="Predicted", ylab="Actual", main="Calibration")
  dev.off(); cat("  calibration.pdf done\n")
} else {
  cat("rms package not available, skipping nomogram\n")
}

# ==================== GO/KEGG Enrichment ====================
cat("\n========== GO/KEGG ==========\n")
tryCatch({
  eg <- bitr(disease_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  go <- enrichGO(eg$ENTREZID, OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.2)
  if (nrow(go) > 0) {
    pdf("19/GO_barplot.pdf", width=10, height=7)
    print(barplot(go, showCategory=20))
    dev.off()
    write.table(as.data.frame(go), "19/GO.txt", sep="\t", quote=FALSE)
  }
  kegg <- enrichKEGG(eg$ENTREZID, organism="hsa", pvalueCutoff=0.05, qvalueCutoff=0.2)
  if (nrow(kegg) > 0) {
    pdf("20/KEGG_barplot.pdf", width=10, height=7)
    print(barplot(kegg, showCategory=20))
    dev.off()
    write.table(as.data.frame(kegg), "20/KEGG.txt", sep="\t", quote=FALSE)
  }
}, error=function(e) cat(sprintf("GO/KEGG error: %s\n", e$message)))

# ==================== STEP 21: TCGA LUAD+LUSC Download ====================
cat("\n========== TCGA: Download ==========\n")
tcga_url <- "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz"
tcga_gz <- "21/tcga_RSEM_gene_tpm.gz"
if (!file.exists(tcga_gz)) {
  cat("Downloading TCGA TPM from UCSC Xena (~700 MB)...\n")
  download.file(tcga_url, tcga_gz, mode="wb")
}
cat("Reading TCGA expression matrix...\n")
tcga <- read.table(gzfile(tcga_gz), header=TRUE, sep="\t", check.names=FALSE, row.names=1, comment.char="")
cat(sprintf("TCGA raw: %d genes x %d samples\n", nrow(tcga), ncol(tcga)))

keep <- grepl("^TCGA-LUAD|^TCGA-LUSC", colnames(tcga))
tcga_lung <- tcga[, keep, drop=FALSE]
cat(sprintf("LUAD+LUSC: %d genes x %d samples\n", nrow(tcga_lung), ncol(tcga_lung)))

sample_type <- ifelse(grepl("-11[A-Z0-9]*-", colnames(tcga_lung)), "Normal", "Tumor")
cat(sprintf("Tumor: %d, Normal: %d\n", sum(sample_type=="Tumor"), sum(sample_type=="Normal")))
saveRDS(list(expr=tcga_lung, type=sample_type), "21/tcga_luad_lusc.rds")

# ==================== STEP 22: GBP1 High/Low ====================
cat("\n========== GBP1 Group Split ==========\n")
gbp1_row <- which(rownames(tcga_lung) == "GBP1")
if (length(gbp1_row) == 0) {
  gbp1_row <- grep("GBP1", rownames(tcga_lung), ignore.case=FALSE)[1]
}
if (length(gbp1_row) == 0 || is.na(gbp1_row)) stop("GBP1 not found in TCGA data")
gbp1_tpm <- as.numeric(tcga_lung[gbp1_row, ])
cat(sprintf("GBP1: row=%s, median=%.2f, range=[%.2f, %.2f]\n",
  rownames(tcga_lung)[gbp1_row], median(gbp1_tpm), min(gbp1_tpm), max(gbp1_tpm)))

gbp1_group <- ifelse(gbp1_tpm > median(gbp1_tpm), "GBP1-High", "GBP1-Low")
names(gbp1_group) <- colnames(tcga_lung)
write.table(data.frame(Sample=colnames(tcga_lung), GBP1_TPM=gbp1_tpm, Group=gbp1_group, Type=sample_type),
  "22/gbp1_groups.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("GBP1-High: %d, GBP1-Low: %d\n", sum(gbp1_group=="GBP1-High"), sum(gbp1_group=="GBP1-Low")))

# ==================== STEP 23: CIBERSORT Deconvolution ====================
cat("\n========== CIBERSORT ==========\n")
if (!requireNamespace("IOBR", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cran.r-project.org")
  BiocManager::install("IOBR", update=FALSE, ask=FALSE)
  library(IOBR)
} else {
  library(IOBR)
}

tcga_expr <- as.matrix(tcga_lung)
mode(tcga_expr) <- "numeric"

if (max(tcga_expr, na.rm=TRUE) < 50) {
  tcga_expr <- 2^tcga_expr - 0.001
  tcga_expr[tcga_expr < 0] <- 0
  cat("Anti-log applied (log2 detected)\n")
}

cib <- deconvo_tme(eset=tcga_expr, method="cibersort", perm=1000)

cib_filtered <- cib[cib$P.value < 0.05, ]
cell_cols <- c("B.cells.naive","B.cells.memory","Plasma.cells","T.cells.CD8","T.cells.CD4.naive",
  "T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper",
  "T.cells.regulatory.Tregs","T.cells.gamma.delta","NK.cells.resting","NK.cells.activated",
  "Monocytes","Macrophages.M0","Macrophages.M1","Macrophages.M2","Dendritic.cells.resting",
  "Dendritic.cells.activated","Mast.cells.resting","Mast.cells.activated","Eosinophils","Neutrophils")
cell_cols <- intersect(cell_cols, colnames(cib_filtered))
cib_mat <- cib_filtered[, cell_cols, drop=FALSE]
write.table(cib_filtered, "23/CIBERSORT-Results.txt", sep="\t", quote=FALSE)
cat(sprintf("CIBERSORT: %d samples pass P<0.05 (from %d input)\n", nrow(cib_mat), nrow(cib)))

# ==================== STEP 24: CIBERSORT Visualization ====================
cat("\n========== CIBERSORT Visualization ==========\n")
common_s <- intersect(rownames(cib_mat), names(gbp1_group))
cat(sprintf("Samples with both CIBERSORT + GBP1 group: %d\n", length(common_s)))

cib_plot <- cib_mat[common_s, , drop=FALSE]
grp <- gbp1_group[common_s]
grp_col <- c("GBP1-High"="#DE2D26", "GBP1-Low"="#3182BD")

# Plot 1: Stacked barplot
pdf("24/barplot.pdf", width=14, height=7)
bar_data <- cib_plot; bar_data$Sample <- rownames(bar_data); bar_data$Group <- grp
bar_melt <- reshape2::melt(bar_data, id.vars=c("Sample","Group"), variable.name="CellType", value.name="Fraction")
bar_melt <- bar_melt[order(bar_melt$Group, bar_melt$Sample), ]
print(ggplot(bar_melt, aes(x=Sample, y=Fraction, fill=CellType)) +
  geom_bar(stat="identity", width=1) +
  facet_grid(~Group, scales="free_x", space="free_x") +
  labs(title="CIBERSORT Immune Cell Proportions by GBP1 Group", y="Fraction", x="") +
  theme_bw(10) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        strip.text=element_text(size=11, face="bold"),
        panel.spacing=unit(0.1,"lines")))
dev.off(); cat("  barplot.pdf done\n")

# Plot 2: Group-diff boxplot (top 10 cell types by mean)
pdf("24/immune.diff.pdf", width=12, height=6)
cell_means <- colMeans(cib_plot); top_cells <- names(sort(cell_means, decreasing=TRUE)[1:10])
diff_data <- reshape2::melt(data.frame(Sample=common_s, Group=grp, cib_plot[, top_cells, drop=FALSE]),
  id.vars=c("Sample","Group"), variable.name="CellType", value.name="Fraction")
print(ggplot(diff_data, aes(x=Group, y=Fraction, fill=Group)) +
  geom_boxplot(outlier.size=0.3, alpha=0.7) +
  facet_wrap(~CellType, scales="free_y", nrow=2) +
  scale_fill_manual(values=grp_col) +
  ggpubr::stat_compare_means(aes(group=Group), label="p.signif", hide.ns=TRUE, size=3) +
  labs(title="Immune Cell Differences: GBP1-High vs Low") + theme_bw(10) +
  theme(axis.text.x=element_text(angle=30, hjust=1)))
dev.off(); cat("  immune.diff.pdf done\n")

# Plot 3: Correlation heatmap
pdf("24/corHeatmap.pdf", width=10, height=8)
cib_cor <- cor(cib_plot, method="spearman")
corrplot::corrplot(cib_cor, method="color", type="upper", tl.cex=0.6, tl.col="black",
  title="Immune Cell Spearman Correlation (TCGA Lung)", mar=c(0,0,2,0))
dev.off(); cat("  corHeatmap.pdf done\n")

# ==================== STEP 25: GBP1-Immune Scatter ====================
cat("\n========== GBP1-Immune Correlation ==========\n")
cib_common <- cib_mat[common_s, , drop=FALSE]
gbp1_vals <- gbp1_tpm[match(common_s, colnames(tcga_lung))]

gbp1_cor <- sapply(colnames(cib_common), function(ct) {
  cor.test(gbp1_vals, cib_common[, ct], method="spearman")
}, simplify=FALSE)
gbp1_r <- sapply(gbp1_cor, `[[`, "estimate")
gbp1_p <- sapply(gbp1_cor, `[[`, "p.value")
gbp1_stats <- data.frame(CellType=names(gbp1_r), Spearman_R=round(gbp1_r, 3), Pvalue=signif(gbp1_p, 3))
gbp1_stats <- gbp1_stats[order(-abs(gbp1_stats$Spearman_R)), ]
write.table(gbp1_stats, "25/gbp1_immune_cor.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("Top GBP1-correlated immune cells:\n")
print(head(gbp1_stats, 10))

pdf("25/gbp1_immune_scatter.pdf", width=12, height=9)
top6 <- head(gbp1_stats$CellType, 6)
par(mfrow=c(2,3), mar=c(4,4,3,1))
for (ct in top6) {
  plot(gbp1_vals, cib_common[, ct],
    xlab="GBP1 TPM", ylab=paste0(ct, " fraction"),
    main=sprintf("%s\nR=%.2f, P=%s", ct, gbp1_stats$Spearman_R[gbp1_stats$CellType==ct],
      format(gbp1_stats$Pvalue[gbp1_stats$CellType==ct], digits=2)),
    pch=19, col=rgb(0,0,0,0.3), cex=0.5)
  abline(lm(cib_common[, ct] ~ gbp1_vals), col="red", lwd=2)
  lines(lowess(gbp1_vals, cib_common[, ct]), col="blue", lwd=2)
}
dev.off(); cat("  gbp1_immune_scatter.pdf done\n")

# ==================== SAVE SUMMARY ====================
cat("\n========== PIPELINE COMPLETE ==========\n")
cat(sprintf("--- TB Pipeline ---\n"))
cat(sprintf("Datasets: %d\n", length(unique(batch_vec))))
cat(sprintf("Samples: %d (C:%d T:%d)\n", ncol(expr), sum(type_vec=="Control"), sum(type_vec=="Treat")))
cat(sprintf("Genes (common): %d\n", nrow(expr)))
cat(sprintf("DEGs: %d\n", nrow(deg)))
cat(sprintf("WGCNA best module (%s): %d genes\n", best_mod, length(best_genes)))
cat(sprintf("Disease genes (DEG ∩ WGCNA): %d\n", length(disease_genes)))
if (exists("top10")) cat(sprintf("Top DEGs: %s\n", paste(rownames(top10), collapse=", ")))
cat(sprintf("\n--- TCGA Immune ---\n"))
cat(sprintf("TCGA LUAD+LUSC: %d samples (T:%d N:%d)\n",
  ncol(tcga_lung), sum(sample_type=="Tumor"), sum(sample_type=="Normal")))
cat(sprintf("CIBERSORT pass P<0.05: %d samples\n", nrow(cib_mat)))
cat(sprintf("GBP1-High: %d, GBP1-Low: %d\n", sum(gbp1_group=="GBP1-High"), sum(gbp1_group=="GBP1-Low")))
if (exists("avg_AUC")) {
  cat(sprintf("\n--- ML Results ---\n"))
  cat(sprintf("Models trained: %d\n", length(model)))
  cat(sprintf("Best model: %s (Test AUC=%.3f)\n", best_method, avg_AUC[1]))
  cat(sprintf("Best model genes (%d): %s\n", length(best_genes), paste(best_genes, collapse=", ")))
}

saveRDS(list(deg=deg, wgcna_genes=best_genes, disease_genes=disease_genes,
             best_mod=best_mod, batch=batch_vec, type=type_vec,
             tcga_expr=tcga_lung, tcga_type=sample_type,
             cibersort=cib_filtered, gbp1_group=gbp1_group, gbp1_immune_cor=gbp1_stats,
             ml_models=if(exists("model")) model else NULL,
             ml_auc=if(exists("AUC_mat")) AUC_mat else NULL,
             ml_best=if(exists("best_method")) best_method else NULL,
             ml_genes=if(exists("best_genes")) best_genes else NULL),
        "pipeline_results.rds")
cat("Saved pipeline_results.rds\n")
