# ---------------------- 环境初始化 & 包安装/加载 ----------------------
rm(list = ls()); gc()
setwd("YOUR_WORKING_DIR")  # 改成你存放数据的目录

# 安装必需包（仅首次运行，增加判断避免重复安装）
pkg_list <- c("rms", "pROC", "ggplot2", "dplyr", "tidyr", "plotROC", "ResourceSelection")
new_pkg <- pkg_list[!pkg_list %in% installed.packages()[,1]]
if (length(new_pkg) > 0) {
  install.packages(new_pkg, dependencies = TRUE, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}

# 加载包
library(rms)                
library(pROC)               
library(ggplot2)            
library(dplyr)              
library(tidyr)              
library(plotROC)            
library(ResourceSelection)  

# ---------------------- 数据读取 & 预处理（核心优化） ----------------------
# 1. 读取核心基因列表（增加异常判断）
if (!file.exists("rf_lasso_svm_gene.csv")) {
  stop("核心基因文件 rf_lasso_svm_gene.csv 不存在！请检查文件路径")
}
hubgene <- read.csv("rf_lasso_svm_gene.csv", header = TRUE)[,1]
hubgene <- hubgene[!is.na(hubgene)]  # 去除NA值
if (length(hubgene) == 0) stop("核心基因列表为空！")
cat("核心基因列表：", paste(hubgene, collapse = ","), "\n")

# 2. 读取表达矩阵和分组文件（增加文件存在性检查）
if (!file.exists("02normal_exp_uniq_gene.csv")) stop("表达矩阵文件不存在！")
if (!file.exists("group.csv")) stop("分组文件不存在！")

rt_raw <- read.csv("02normal_exp_uniq_gene.csv", header = TRUE, check.names = FALSE)
pdata <- read.csv("group.csv", header = TRUE, check.names = FALSE)  

# 3. 构建分析数据集（样本×核心基因+结局）
gene_col <- colnames(rt_raw)[1]  # 基因名列（第一列）
# 提取核心基因行（优化匹配逻辑）
rt_filtered <- rt_raw[rt_raw[[gene_col]] %in% hubgene, ]  
if (nrow(rt_filtered) == 0) stop("无匹配的核心基因！请检查基因名格式")
rownames(rt_filtered) <- rt_filtered[[gene_col]]          
rt_filtered <- rt_filtered[, -which(colnames(rt_filtered) == gene_col), drop = FALSE]  
exp_matrix <- t(rt_filtered)  # 转置：样本×基因
exp_df <- as.data.frame(exp_matrix, stringsAsFactors = FALSE)

# 4. 构建二元结局变量（修复因子转换bug）
sample_ids <- intersect(rownames(exp_df), pdata$ID)
if (length(sample_ids) == 0) stop("无匹配的样本ID！请检查表达矩阵和分组文件的样本名")
if (length(sample_ids) < 10) warning("样本量过少（<10），模型结果可能不可靠！")

# 对齐样本顺序
exp_df <- exp_df[sample_ids, , drop = FALSE]
pdata <- pdata[match(sample_ids, pdata$ID), ]
# 构建结局变量（统一转换为数值型后再转因子）
exp_df$y <- ifelse(pdata$group == "EM", 1, 0)
exp_df$y <- factor(exp_df$y, levels = c(0, 1), labels = c("Normal", "EM"))  # 增加标签便于解读

# 5. 数据清洗（优化标准化逻辑）
exp_df <- exp_df %>% na.omit()  # 去除缺失值
# 标准化核心基因表达量（增加存在性判断）
model_genes <- hubgene[hubgene %in% colnames(exp_df)]
if (length(model_genes) == 0) stop("无有效核心基因！")

for (gene in model_genes) {
  exp_df[[gene]] <- as.numeric(exp_df[[gene]])  # 确保为数值型
  exp_df[[gene]] <- scale(exp_df[[gene]])[,1]  # Z-score标准化
}

# 输出数据基本信息
cat("分析数据集维度：", dim(exp_df), "\n")
cat("结局分布（0=Normal,1=EM）：", table(exp_df$y), "\n")

# ---------------------- 构建rms环境（lrm必需） ----------------------
dd <- datadist(exp_df)
options(datadist = "dd")

# ---------------------- 1. 多因素Logistic回归模型（增强稳定性） ----------------------
# 构建公式（避免基因名特殊字符问题）
formula_str <- paste("y ~", paste(model_genes, collapse = "+"))
formula_obj <- as.formula(formula_str)

# 训练lrm模型（增加权重优化小样本）
model <- lrm(
  formula = formula_obj, 
  data = exp_df, 
  x = TRUE, 
  y = TRUE,
  weights = rep(1, nrow(exp_df))  # 权重，小样本可调整
)
# 输出模型结果（美化格式）
cat("\n===== 多因素Logistic回归模型结果 =====", "\n")
print(summary(model))  

# ---------------------- 2. Nomogram（列线图，优化可视化） ----------------------
# 构建Nomogram（调整参数提升美观度）
nom <- nomogram(
  model,
  fun = plogis,  # 逻辑回归概率转换
  lp = FALSE,
  funlabel = "Probability of Endometriosis (EM)",
  maxscale = 100,
  abbrev = FALSE,  # 不缩写基因名
  fun.at = seq(0, 1, 0.2)  # 概率刻度
)

# 保存Nomogram为TIF（高清，优化尺寸）
tiff("Nomogram.tif", width = 14, height = 9, units = "in", res = 300)
plot(nom, 
     xfrac = 0.15,        # 刻度标签宽度
     cex.axis = 0.9,      # 轴字体大小
     cex.main = 1.2,      # 标题字体大小
     col.grid = gray(c(0.85, 0.95)),  # 网格颜色更柔和
     lwd.grid = 0.5)      # 网格线宽度
dev.off()
cat("Nomogram已保存为Nomogram.tif\n")

# ---------------------- 3. ROC曲线（单基因+多基因模型，修复AUC计算） ----------------------
# ---------------------- 3. ROC曲线（单基因+多基因模型，彻底修复AUC+二分类问题） ----------------------
# 3.1 重新处理结局变量（核心：确保二分类且无NA）
# 从原始分组文件重建结局变量（绕过因子转换坑）
exp_df$y_num <- ifelse(pdata$group[match(rownames(exp_df), pdata$ID)] == "EM", 1, 0)
# 强制清洗：仅保留0/1，去除NA
exp_df$y_num <- as.numeric(exp_df$y_num)
exp_df$y_num[is.na(exp_df$y_num)] <- 0  # 可根据研究设计改为：exp_df <- exp_df[!is.na(exp_df$y_num), ]
exp_df <- exp_df[exp_df$y_num %in% c(0,1), ]  # 过滤非0/1的异常值

# 验证二分类（提前终止错误）
if (length(unique(exp_df$y_num)) < 2) {
  stop("结局变量不是二分类！请检查分组文件中是否同时包含Normal和EM样本")
}

# 3.2 多基因模型预测概率（确保与清洗后数据匹配）
exp_df$pred_prob <- predict(model, newdata = exp_df, type = "fitted")  # rms包预测概率
exp_df <- exp_df[!is.na(exp_df$pred_prob), ]  # 过滤预测概率的NA值

# 3.3 多基因ROC曲线（强制使用pROC包，修复AUC计算）
# 核心：显式调用pROC::roc + pROC::auc，避免函数冲突
roc_obj <- pROC::roc(
  response = exp_df$y_num,
  predictor = exp_df$pred_prob,
  levels = c(0, 1),        # 明确：0=阴性(Normal)，1=阳性(EM)
  direction = "<",         # 预测值越大，阳性概率越高（相反则改">"）
  ci = TRUE,               # 计算95%CI
  quiet = TRUE             # 关闭冗余输出
)

# 计算AUC（两种方式任选，均无报错）
auc_value <- roc_obj$auc        # 推荐：直接提取属性，完全避开auc()函数
# auc_value <- pROC::auc(roc_obj)  # 备选：显式调用pROC的auc函数

ci_auc <- pROC::ci.auc(roc_obj)  # 计算AUC的95%CI

# 保存ROC曲线
tiff("ROC_Curve.tif", width = 9, height = 8, units = "in", res = 300)
pROC::plot.roc(
  roc_obj,
  main = paste0("ROC Curve of Multi-Gene Model\nAUC = ", round(auc_value, 3), 
                " (95%CI: ", round(ci_auc[1],3), "-", round(ci_auc[3],3), ")"),
  col = "#E74C3C", lwd = 2.5,
  print.auc = TRUE,  
  print.auc.x = 0.6, print.auc.y = 0.2, 
  print.auc.cex = 1.2,
  ci = TRUE,  
  ci.col = "#BDC3C7",
  legacy.axes = TRUE  # 传统坐标轴（1-Specificity为x轴）
)
abline(a = 0, b = 1, lty = 2, col = "#7F8C8D", lwd = 1.5)
legend("bottomright", legend = "Multi-gene model", col = "#E74C3C", lwd = 2.5, bty = "n")
dev.off()

# 3.4 单基因ROC曲线（动态布局+异常处理）
model_genes <- model_genes[model_genes %in% colnames(exp_df)]  # 验证基因存在性
n_plot <- min(length(model_genes), 4)  # 最多显示4个基因
if (n_plot > 0) {
  tiff("Single_Gene_ROC.tif", width = 10, height = 8, units = "in", res = 300)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))  # 调整边距
  for (i in 1:n_plot) {
    gene <- model_genes[i]
    # 清洗单基因数据：去除NA
    gene_data <- exp_df[!is.na(exp_df[[gene]]), ]
    # 验证二分类样本充足
    if (length(unique(gene_data$y_num)) < 2) {
      warning(paste(gene, "无足够二分类样本，跳过该基因ROC绘制"))
      next
    }
    # 单基因ROC（强制使用pROC包）
    roc_single <- pROC::roc(
      response = gene_data$y_num,
      predictor = gene_data[[gene]],
      levels = c(0, 1),
      ci = TRUE,
      quiet = TRUE
    )
    auc_single <- roc_single$auc  # 直接提取AUC，避免函数冲突
    # 绘制单基因ROC
    pROC::plot.roc(
      roc_single,
      main = paste(gene, "(AUC =", round(auc_single, 3), ")"),
      col = "#3498DB", lwd = 2, 
      print.auc = TRUE, 
      legacy.axes = TRUE
    )
    abline(a = 0, b = 1, lty = 2, col = "#7F8C8D")
  }
  dev.off()
}
cat("ROC曲线已保存为ROC_Curve.tif和Single_Gene_ROC.tif\n")
cat("多基因模型AUC：", round(auc_value, 3), "(95%CI: ", round(ci_auc[1],3), "-", round(ci_auc[3],3), ")\n")
# ---------------------- 4. DCA（决策曲线分析，优化计算逻辑） ----------------------
# 4.1 定义DCA计算函数（修复分母为0问题）
calculate_dca <- function(obs, pred, thresholds = seq(0, 1, 0.01)) {
  obs <- as.numeric(obs)
  dca_result <- data.frame(
    threshold = thresholds,
    net_benefit = 0,
    nb_all = 0,
    nb_none = 0
  )
  n <- length(obs)
  pos <- sum(obs == 1)  # 阳性样本数
  
  for (i in 1:length(thresholds)) {
    p <- thresholds[i]
    if (p %in% c(0, 1)) next  # 跳过极端值
    
    pred_pos <- sum(pred >= p)
    tp <- sum(obs == 1 & pred >= p)
    fp <- sum(obs == 0 & pred >= p)
    
    # 计算净获益（增加防除零保护）
    if (p >= 1 - 1e-6) {
      nb <- 0
    } else {
      nb <- (tp / n) - (fp / n) * (p / (1 - p))
    }
    dca_result$net_benefit[i] <- nb
    # 所有干预的净获益
    dca_result$nb_all[i] <- (pos / n) - ((n - pos) / n) * (p / (1 - p))
  }
  return(dca_result)
}

# 4.2 计算DCA数据
dca_data <- calculate_dca(
  obs = exp_df$y_num,
  pred = exp_df$pred_prob,
  thresholds = seq(0, 1, by = 0.01)
)

# 4.3 绘制DCA曲线（优化可视化）
tiff("DCA_Curve.tif", width = 9, height = 8, units = "in", res = 300)
ggplot(dca_data, aes(x = threshold)) +
  # 模型净获益曲线
  geom_line(aes(y = net_benefit), color = "#E74C3C", linewidth = 1.5) +
  # 所有干预曲线
  geom_line(aes(y = nb_all), color = "#95A5A6", linewidth = 1, linetype = 1) +
  # 无干预曲线
  geom_line(aes(y = nb_none), color = "#7F8C8D", linewidth = 1, linetype = 2) +
  # 美化
  labs(title = "Decision Curve Analysis (DCA)",
       x = "Threshold Probability of EM",
       y = "Net Benefit") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(min(dca_data$net_benefit, na.rm = T) - 0.05, 
                                max(dca_data$nb_all, na.rm = T) + 0.05)) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_blank())
dev.off()
cat("DCA曲线已保存为DCA_Curve.tif\n")

# ---------------------- 5. Calibration Curve（校准曲线，增强可靠性） ----------------------
# 5.1 简单版校准曲线
tiff("Calibration_Curve.tif", width = 9, height = 8, units = "in", res = 300)
ggplot(data.frame(pred = exp_df$pred_prob, obs = exp_df$y_num), 
       aes(x = pred, y = obs)) +
  # 平滑校准线（逻辑回归）
  geom_smooth(method = "glm", method.args = list(family = binomial), 
              se = TRUE, color = "#E74C3C", linewidth = 1.2, fill = "#FADBD8") +
  # 理想参考线
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#7F8C8D", linewidth = 1) +
  labs(title = "Calibration Curve of Multi-Gene Model",
       x = "Predicted Probability of EM",
       y = "Observed Probability of EM") +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
dev.off()

# 5.2 详细版校准曲线（含Hosmer-Lemeshow检验，优化分组）
# 动态调整分组数（适配样本量）
n_groups <- ifelse(nrow(exp_df) < 50, 5, 10)
hl_test <- hoslem.test(exp_df$y_num, exp_df$pred_prob, g = n_groups)  
cat("Hosmer-Lemeshow检验P值：", round(hl_test$p.value, 3), "\n")

# 1. 生成分组
exp_df$pred_group <- cut(exp_df$pred_prob, breaks = n_groups, include.lowest = TRUE)

# 2. 计算分组均值（原生aggregate）
mean_data <- aggregate(
  cbind(pred_prob, y_num) ~ pred_group,
  data = exp_df,
  FUN = function(x) mean(x, na.rm = TRUE)
)
colnames(mean_data)[colnames(mean_data) == "y_num"] <- "obs_prob"  # 重命名

# 3. 统计分组样本数（原生table）
count_data <- as.data.frame(table(exp_df$pred_group, useNA = "no"))
colnames(count_data) <- c("pred_group", "n")

# 4. 合并数据
cal_data <- merge(mean_data, count_data, by = "pred_group", all.x = TRUE) %>%
  na.omit() %>%
  filter(n > 0)  # 过滤空分组
# 绘制详细校准曲线
tiff("Calibration_Curve_Detailed.tif", width = 9, height = 8, units = "in", res = 300)
ggplot(cal_data, aes(x = pred_prob, y = obs_prob)) +
  geom_point(size = 4, color = "#E74C3C", alpha = 0.8) +  # 分组点
  geom_errorbar(aes(ymin = obs_prob - sqrt(obs_prob*(1-obs_prob)/n), 
                    ymax = obs_prob + sqrt(obs_prob*(1-obs_prob)/n)),
                width = 0.02, color = "#E74C3C") +  # 误差棒
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", linewidth = 1.2) +  # 校准线
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#7F8C8D") +  # 理想线
  labs(title = paste0("Calibration Curve (HL test P = ", round(hl_test$p.value, 3), ")"),
       x = "Mean Predicted Probability",
       y = "Mean Observed Probability") +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
dev.off()
cat("校准曲线已保存为Calibration_Curve.tif和Calibration_Curve_Detailed.tif\n")

# ---------------------- 结果汇总（美化输出） ----------------------
cat("\n==================== 核心基因预测模型分析结果汇总 ====================", "\n")
cat("1. 纳入核心基因数：", length(model_genes), "\n")
cat("2. 分析样本量：", nrow(exp_df), " (Normal =", sum(exp_df$y_num == 0), ", EM =", sum(exp_df$y_num == 1), ")\n")
cat("3. 多基因模型AUC：", round(auc_value, 3), "(95%CI: ", round(ci_auc[1],3), "-", round(ci_auc[3],3), ")\n")
cat("4. Hosmer-Lemeshow检验P值：", round(hl_test$p.value, 3), "\n")
cat("   → P>0.05说明模型校准度良好，P<0.05提示模型校准不足\n")
cat("5. 输出文件清单：\n")
cat("   - Nomogram.tif（列线图）\n")
cat("   - ROC_Curve.tif（多基因ROC曲线）\n")
cat("   - Single_Gene_ROC.tif（单基因ROC曲线）\n")
cat("   - DCA_Curve.tif（决策曲线分析）\n")
cat("   - Calibration_Curve.tif（简易校准曲线）\n")
cat("   - Calibration_Curve_Detailed.tif（详细校准曲线）\n")
cat("6. 文件保存路径：", getwd(), "\n")