#------------------------------ LASSO基因数精准控制版 ---------------
rm(list = ls()); gc()

# 安装必要包（仅首次运行）
pkg_list <- c("survival", "glmnet", "ggsci", "patchwork", "randomForest", 
              "VennDiagram", "ggpubr", "limma", "data.table", "tibble", "e1071", "caret")
install.packages(pkg_list[!pkg_list %in% installed.packages()[,1]])

# 加载包
library(survival); library(glmnet); library(ggplot2); library(ggsci); library(patchwork)
library(limma); library(data.table); library(tibble); library(randomForest); library(VennDiagram)
library(ggpubr); library(e1071); library(caret)

# 工作路径（修改为你的实际路径）
setwd("YOUR_WORKING_DIR")  # 改成你存放数据的目录

# ========== 步骤1：数据读取与预处理（核心整合） ==========
# 1.1 读取表达矩阵（GeneSymbol为行名）
exp_df <- read.csv("02normal_exp_uniq_gene.csv", check.names = F)
rownames(exp_df) <- exp_df$GeneSymbol
exp_df <- exp_df[, -which(colnames(exp_df) == "GeneSymbol"), drop = F]

# 1.2 读取共表达基因列表并精准匹配（统一基因名格式）
cogene_list <- read.csv("cogene.csv", header = T)[,1]
clean_gene <- function(gene_vec) {  # 基因名清洗函数
  gsub("///.*|\\-.*|\\_.*", "", tolower(gene_vec))
}
exp_gene_clean <- clean_gene(rownames(exp_df))
cogene_clean <- clean_gene(cogene_list)
match_idx <- match(cogene_clean, exp_gene_clean) %>% na.omit()
exp_mat <- exp_df[match_idx, ]  # 仅保留匹配的基因
rownames(exp_mat) <- rownames(exp_df)[match_idx]  # 恢复原始基因名

# 1.3 读取分组文件并匹配样本
group_df <- read.csv("group.csv", header = T, check.names = F)
group_df <- group_df[order(group_df$group == "Normal", decreasing = T), ]
matched_samples <- intersect(group_df$ID, colnames(exp_mat))
if(length(matched_samples) == 0) stop("无匹配样本名！请检查分组文件和表达矩阵的样本名")
exp_mat <- exp_mat[, matched_samples]
group_df <- group_df[match(matched_samples, group_df$ID), ]

# 1.4 数据格式转换（样本为行，基因为列）+ 修复基因名特殊字符
exp_mat_t <- as.data.frame(t(exp_mat))
colnames(exp_mat_t) <- gsub('-', "_", colnames(exp_mat_t))  # 避免模型报错

# 1.5 保存标准化文件（供后续模型调用）
write.table(group_df, 'sample.txt', sep='\t', row.names = F, col.names = F, quote = F)
write.table(exp_mat_t, 'cogenedata.txt', sep = '\t', row.names = T, col.names = NA, quote = F)

# ========== 步骤2：模型输入数据统一准备 ==========
# 读取标准化表达矩阵
rt <- read.table("cogenedata.txt", header=T, sep="\t", check.names=F, row.names=1, fill=T)
data_mat <- matrix(as.numeric(as.matrix(rt)), nrow=nrow(rt), dimnames=dimnames(rt))
data_mat <- avereps(data_mat)  # 合并重复样本的表达值

# 读取分组标签并构建二分类因子
sample_group <- read.table("sample.txt", sep='\t', header=F, check.names=F)
colnames(sample_group) <- c("ID","group")
sample_group <- sample_group[match(rownames(data_mat), sample_group$ID), ]
y_bin <- ifelse(sample_group$group == "Normal", 0, 1)  # 数值型标签
y_factor <- factor(y_bin, levels = c(0,1), labels = c("Normal","EM"))  # 因子型标签

# 缺失值/无穷值填充（保证模型运行）
x_imputed <- apply(data_mat, 2, function(col) {
  col[is.na(col) | is.infinite(col)] <- mean(col, na.rm = TRUE)
  return(col)
}) %>% as.matrix()

# ========== 步骤3：LASSO回归筛选基因（精准控数+稳定性验证） ==========
set.seed(1221)
# 动态调整交叉验证折数（适配小样本）
nfolds <- ifelse(min(table(y_factor)) >=10, 10, min(table(y_factor)))
nfolds <- max(nfolds, 2)  # 至少2折

# 【核心调整1：收紧惩罚强度（恢复纯LASSO+适度λ范围）】
cvfit <- cv.glmnet(
  x = x_imputed,
  y = y_factor,
  family = "binomial", 
  nlambda=150,                  # 减少λ搜索数量（从200→150）
  lambda.min.ratio = 1e-4,      # 恢复默认λ下限（从1e-5→1e-4）
  alpha=1,                      # 恢复纯LASSO（α=1），惩罚更严格
  nfolds = nfolds,     
  grouped = FALSE,              # 适配样本数不均的情况
  type.measure = "deviance"     # 小样本更稳定
)
fit <- glmnet(x_imputed, y_factor, family = "binomial", lambda = cvfit$lambda)

# 【核心调整2：选择适中的λ（lambda.1se，平衡基因数和泛化能力）】
# lambda.1se：误差在最小误差1个标准误内的最大λ，基因数更少、模型更稳健
custom_lambda <- cvfit$lambda.1se  
# 备选：若仍多，用lambda.min（比1se宽松，但比×0.5严格）
# custom_lambda <- cvfit$lambda.min  
lasso_coef <- coef(fit, s = custom_lambda)

# 提取非零系数基因（移除截距项）
lasso_gene_raw <- row.names(lasso_coef)[which(lasso_coef != 0)]
lasso_gene_raw <- lasso_gene_raw[!lasso_gene_raw %in% "(Intercept)"]

# 【核心调整3：收紧稳定性验证阈值（保留≥70%频率的基因）】
boot_rep <- 100  # 抽样次数
boot_gene_list <- replicate(boot_rep, {
  idx <- sample(1:nrow(x_imputed), replace = TRUE)  # 自助法抽样
  cvfit_boot <- cv.glmnet(
    x_imputed[idx,], y_factor[idx], family = "binomial",
    nlambda=150, lambda.min.ratio=1e-4, alpha=1, nfolds=nfolds
  )
  coef_boot <- coef(cvfit_boot, s = cvfit_boot$lambda.1se)  # 同步λ选择
  row.names(coef_boot)[which(coef_boot != 0) & row.names(coef_boot)!="(Intercept)"]
})
# 统计基因出现频率，保留≥70%频率的核心基因（原50%→70%）
gene_freq <- table(unlist(boot_gene_list))
stable_lasso_gene <- names(gene_freq)[gene_freq >= boot_rep * 0.7]
lasso_gene <- stable_lasso_gene  # 最终核心基因列表

# 保存LASSO结果
write.table(
  data.frame(Gene=lasso_gene), 
  "lassoset.txt", 
  sep="\t", quote=F, col.names=F, row.names=F
)

# 【可视化：LASSO曲线（标注lambda.1se）】
tiff("01lasso_optimized.tif", height = 5, width = 10, units = "in", res = 300)
par(mfrow = c(1, 2))  # 1行2列布局

# 图1：系数变化曲线（标注lambda.1se）
plot(
  fit, xvar = 'lambda', 
  main = "LASSO系数变化曲线（纯L1惩罚）",
  xlab = "log(λ)", ylab = "系数值", lwd = 1.5
)
abline(v = log(custom_lambda), col = "red", lty = 2, lwd = 2)
text(
  x = log(custom_lambda), y = max(fit$beta), 
  labels = paste0("λ.1se\n(logλ=", round(log(custom_lambda), 2), ")"), 
  col = "red", pos = 4, cex = 0.8
)

# 图2：交叉验证曲线（突出lambda.1se）
plot(
  cvfit, main = "LASSO交叉验证误差曲线",
  xlab = "log(λ)", ylab = "交叉验证偏差", lwd = 1.5
)
abline(v = log(cvfit$lambda.min), col = "black", lty = 2, lwd = 2)
abline(v = log(cvfit$lambda.1se), col = "red", lty = 1, lwd = 2)
# 添加图例
legend(
  "topright", 
  legend = c(
    paste0("λ.min (", round(cvfit$lambda.min, 4), ")"),
    paste0("λ.1se (", round(cvfit$lambda.1se, 4), ")")
  ),
  col = c("black", "red"), lty = c(2,1), lwd = 2, cex = 0.8
)
dev.off()

# 【验证输出：LASSO基因数控制效果】
cat("===== LASSO基因数控制验证 =====\n")
cat("纯LASSOα值：", cvfit$alpha, "\n")
cat("λ选择：", ifelse(custom_lambda == cvfit$lambda.1se, "lambda.1se", "lambda.min"), "\n")
cat("原始筛选基因数：", length(lasso_gene_raw), "\n")
cat("核心稳定基因数：", length(lasso_gene), "\n")

# ========== 步骤4：随机森林（RF）筛选基因 ==========
set.seed(1222)
# 训练RF模型（计算基因重要性）
rf_model <- randomForest(
  x = data_mat, 
  y = factor(ifelse(sample_group$group=="Normal","con","treat")),
  ntree=500,
  importance=TRUE
)
# 最优树数重训练
rf_opt <- randomForest(
  x = data_mat, 
  y = factor(ifelse(sample_group$group=="Normal","con","treat")),
  ntree=which.min(rf_model$err.rate[,1]),
  importance=TRUE
)

# 提取并整理基因重要性
rf_importance <- as.data.frame(importance(rf_opt)) %>% 
  rownames_to_column("Gene") %>% 
  select(Gene, importance=MeanDecreaseGini) %>% 
  arrange(desc(importance))
rf_gene <- subset(rf_importance, importance>0.3)$Gene  # 收紧RF阈值（0.25→0.3）

# 保存RF结果
write.table(rf_importance, "rfGenes.xls", sep="\t", quote=F, col.names=T, row.names=F)
write.table(rf_gene, "rfset.txt", sep="\t", quote=F, col.names=F, row.names=F)

# RF可视化1：误差曲线
tiff("02RF_forest.tif", width=6, height=6, units = "in", res = 300)
plot(rf_model, main="Random Forest 误差曲线", lwd=2)
dev.off()

# RF可视化2：TOP20重要基因
rf_top20 <- head(rf_importance, 20)
p_rf <- ggdotchart(
  rf_top20, x = "Gene", y = "importance", color = "importance",
  sorting = "descending", add = "segments", add.params = list(color = "lightgray", size = 2),
  dot.size = 6, font.label = list(color = "white", size = 9, vjust = 0.5),
  ggtheme = theme_bw(), rotate=TRUE
) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  gradient_color(palette = c(pal_npg()(2)[2], pal_npg()(2)[1])) + 
  grids()
tiff("03RF_importance.tif", width=8, height=6, units = "in", res = 300) 
print(p_rf); dev.off()

# ========== 步骤5：SVM筛选基因（置换检验计算重要性） ==========
set.seed(1223)
# 数据标准化（SVM必须）
data_scaled <- preProcess(data_mat, method = c("center", "scale")) %>% predict(data_mat)
svm_df <- as.data.frame(data_scaled) %>% 
  mutate(group = factor(ifelse(sample_group$group=="Normal","con","treat")))

# 训练SVM模型（径向基核）
svm_model <- svm(
  group ~ ., 
  data = svm_df, 
  kernel = "radial", 
  cost = 1,
  gamma = 1/ncol(data_scaled),
  probability = TRUE, 
  scale = FALSE
)

# 置换检验计算特征重要性（核心修正）
original_acc <- mean(predict(svm_model, svm_df) == svm_df$group)
svm_importance <- lapply(colnames(data_scaled), function(gene) {
  df_perm <- svm_df; df_perm[[gene]] <- sample(df_perm[[gene]])  # 打乱特征
  perm_acc <- mean(predict(svm_model, df_perm) == df_perm$group)
  data.frame(Gene = gene, importance = original_acc - perm_acc)
}) %>% bind_rows() %>% arrange(desc(importance))

# 收紧SVM筛选阈值（前40%→原均值/前50%）
svm_gene <- subset(svm_importance, importance > quantile(svm_importance$importance, 0.6))$Gene

# 保存SVM结果
write.table(svm_importance, "svmGenes.xls", sep="\t", quote=F, col.names=T, row.names=F)
write.table(svm_gene, "svmset.txt", sep="\t", quote=F, col.names=F, row.names=F)

# SVM可视化：TOP20重要基因
svm_top20 <- head(svm_importance, 20)
p_svm <- ggdotchart(
  svm_top20, x = "Gene", y = "importance", color = "importance",
  sorting = "descending", add = "segments", add.params = list(color = "lightgray", size = 2),
  dot.size = 6, font.label = list(color = "white", size = 9, vjust = 0.5),
  ggtheme = theme_bw(), rotate=TRUE
) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  gradient_color(palette = c(pal_aaas()(2)[2], pal_aaas()(2)[1])) + 
  grids()
tiff("05SVM_importance.tif", width=8, height=6, units = "in", res = 300) 
print(p_svm); dev.off()

# ========== 步骤6：三方法交集分析 + 可视化 ==========
# 读取各方法筛选的基因列表
g1 <- read.table('lassoset.txt', sep='\t', header = F)[,1]  # LASSO（控数后）
g2 <- read.table('rfset.txt', sep='\t', header = F)[,1]     # RF
g3 <- read.table('svmset.txt', sep='\t', header = F)[,1]    # SVM

# 计算三方法交集基因
common_gene <- Reduce(intersect, list(g1, g2, g3))

# 绘制三元韦恩图
tiff('04rf_lasso_svm_venn.tif', width = 6, height = 6, units = "in", res = 300) 
venn.plot <- venn.diagram(
  x = list(Random_Forest = g2, Lasso = g1, SVM = g3),
  filename = NULL,
  col = "black",
  fill = c("lightblue", "orangeRed", "lightgreen"),  # 三色区分
  alpha = 0.5,
  label.col = "black",
  cex = 1.2,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.col = 'black',
  cat.pos = c(150, 210, 360),  # 标签位置
  cat.dist = c(0.05, 0.05, 0.05)
)
grid.draw(venn.plot); dev.off()

# 保存交集基因并输出最终统计
write.csv(data.frame(common_gene = common_gene), 'rf_lasso_svm_gene.csv', row.names = F)
cat("\n===== 最终分析结果（基因数控制版）=====\n")
cat("LASSO核心基因数：", length(g1), "\n")
cat("RF核心基因数（importance>0.3）：", length(g2), "\n")
cat("SVM核心基因数（前40%）：", length(g3), "\n")
cat("三方法交集基因数：", length(common_gene), "\n")
if(length(common_gene) > 0) cat("交集核心基因：", paste(common_gene, collapse = ", "), "\n")

# ===================== 附加：10折CV 性能曲线（Fig 1G 右）=====================
# 说明：上面 LASSO/RF/SVM/三法交集 一个字没动，基因结果完全不变。
#       这里只读取已算好的 svm_importance / data_scaled，画 accuracy/errors 曲线。
set.seed(1234)   # 仅让 CV 曲线可复现，不影响上面的基因筛选
ranked_genes <- svm_importance$Gene          # 置换重要性排序（最重要在前）
X_cv <- data_scaled
y_cv <- factor(ifelse(sample_group$group == "Normal", "con", "treat"), levels = c("con", "treat"))
top_n <- 20

fold_idx <- sample(rep(1:10, length.out = nrow(X_cv)))
cv_errors <- c()
for (k in 1:top_n) {
  cur <- ranked_genes[1:k]
  fold_err <- c()
  for (fold in 1:10) {
    tr <- which(fold_idx != fold); te <- which(fold_idx == fold)
    m <- svm(x = X_cv[tr, cur, drop = FALSE], y = y_cv[tr], kernel = "radial", cost = 1, gamma = 1/length(cur))
    pred <- predict(m, X_cv[te, cur, drop = FALSE])
    fold_err <- c(fold_err, sum(pred != y_cv[te]) / length(te))
  }
  cv_errors[k] <- mean(fold_err)
  cat("  CV 特征数=", k, "，误差=", round(cv_errors[k], 3), "\n")
}
best_k <- which.min(cv_errors)
error_data <- data.frame(Feature_Num = 1:top_n, CV_Error = cv_errors)

pdf("errors.pdf", width = 8, height = 6)
print(ggplot(error_data, aes(Feature_Num, CV_Error)) +
  geom_line(color = "#2E86AB", linewidth = 1) +
  geom_point(color = "#A23B72", size = 2) +
  geom_vline(xintercept = best_k, color = "red", linewidth = 0.8) +
  labs(x = "Number of Features", y = "Cross-Validation Error",
       title = "10-Fold CV Error vs Feature Number") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                     axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  xlim(1, top_n))
dev.off()

pdf("accuracy.pdf", width = 8, height = 6)
print(ggplot(error_data, aes(Feature_Num, 1 - CV_Error)) +
  geom_line(color = "#2E86AB", linewidth = 1) +
  geom_point(color = "#A23B72", size = 2) +
  geom_vline(xintercept = best_k, color = "red", linewidth = 0.8) +
  labs(x = "Number of Features", y = "Cross-Validation Accuracy",
       title = "10-Fold CV Accuracy vs Feature Number") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                     axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  xlim(1, top_n))
dev.off()

write.table(error_data, "SVM-CV.cv_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("最优特征数（最小CV误差）：", best_k, "\n")
