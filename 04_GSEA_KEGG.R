# 安装依赖包（按需取消注释）
if (!requireNamespace("BiocManager", quietly = TRUE))		
  install.packages("BiocManager")		
# BiocManager::install("limma")		
# BiocManager::install("org.Hs.eg.db")		
# BiocManager::install("DOSE")		
# BiocManager::install("clusterProfiler")		
# BiocManager::install("enrichplot")		
# BiocManager::install("dplyr")  

# 加载所有必需包
library(limma)		
library(org.Hs.eg.db)		
library(clusterProfiler)		
library(enrichplot)		
library(dplyr)  

# ===================== 配置文件路径 =====================
expFile="merge.normalize.txt"           # 表达数据文件（第一列是基因SYMBOL）
gmtFile="c2.cp.kegg.Hs.symbols.gmt"    # KEGG基因集（SYMBOL格式）
setwd("YOUR_WORKING_DIR")  # 改成你存放数据的目录
target_gene <- "PDE2A"                  # 替换为PDE2A基因

# ===================== 读取并整理表达数据 =====================
rt <- read.table(
  expFile, 
  header = TRUE, 
  sep = "\t", 
  check.names = FALSE, 
  stringsAsFactors = FALSE,
  fill = TRUE  
)		

rt <- as.matrix(rt)		
rownames(rt) <- rt[, 1]		       
exp <- rt[, 2:ncol(rt)]		

dimnames <- list(rownames(exp), colnames(exp))		
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)		
data <- avereps(data)		

data <- data[rowMeans(data) > 0, ]   

# ===================== 筛选处理组样品 =====================
treat_flag <- grepl("_Treat$", colnames(data))  # 匹配列名以"_Treat"结尾的样品
data_treat <- data[, treat_flag, drop = FALSE]  # 筛选处理组，drop=FALSE避免单列时转为向量

if (ncol(data_treat) >= 2) {
  data <- data_treat  
  cat("✅ 成功筛选出", ncol(data), "个处理组样品！\n")
} else {
  # 修正：特殊符号用引号包裹，字符串拼接补全
  warning(paste0("⚠️ 处理组样品数仅", ncol(data_treat), "个，将使用所有样品进行分析！"))
  # 保留所有样品（不筛选）
  data <- data
}

# ===================== 检查目标基因 =====================
if (!target_gene %in% rownames(data)) {
  stop(paste0("❌ 目标基因", target_gene, "不在表达矩阵中，请检查基因名拼写！"))
}
cat(paste0("✅ 确认", target_gene, "存在于表达矩阵中，开始分析...\n"))

# ===================== 高低表达分组 + logFC =====================
gene_exp <- data[target_gene, ]          
median_exp <- median(gene_exp)           

dataL <- data[, gene_exp < median_exp, drop = FALSE]  
dataH <- data[, gene_exp >= median_exp, drop = FALSE] 

if (ncol(dataL) < 1 | ncol(dataH) < 1) {
  stop(paste0("❌ 基因", target_gene, "高低表达组样品数不足，无法分析！"))
}
cat(paste0("✅ ", target_gene, "高低表达分组完成：低表达组", ncol(dataL), "个样品，高表达组", ncol(dataH), "个样品\n"))

meanL <- rowMeans(dataL)	
meanH <- rowMeans(dataH)	
meanL[meanL < 1e-5] <- 1e-5	
meanH[meanH < 1e-5] <- 1e-5	
logFC <- meanH - meanL	       

# ===================== 整理GSEA输入 =====================
logFC <- sort(logFC, decreasing = TRUE)  
logFC <- logFC[!is.na(logFC)]  

if (length(logFC) < 50) {
  stop(paste0("❌ 基因", target_gene, "有效logFC基因数<50，GSEA分析无意义！"))
}
cat(paste0("✅ 有效分析基因数：", length(logFC), "个\n"))

# ===================== 读取基因集 =====================
gmt <- read.gmt(gmtFile)	
if (nrow(gmt) == 0) {
  stop("❌ 基因集文件c2.cp.kegg.Hs.symbols.gmt为空！")
}
cat(paste0("✅ 成功读取KEGG基因集，共", nrow(gmt), "个通路\n"))

# ===================== GSEA分析 =====================
kk <- GSEA(
  logFC,                
  TERM2GENE = gmt,      
  pvalueCutoff = 1,     
  pAdjustMethod = "fdr",
  minGSSize = 10,       
  maxGSSize = 500,      
  verbose = FALSE       
)

if (is.null(kk)) {
  stop(paste0("❌ 基因", target_gene, "无GSEA富集结果！"))
}

# ===================== 筛选显著通路（严格筛选p<0.05） =====================
kkTab <- as.data.frame(kk)	
kkTab_filtered <- kkTab[kkTab$pvalue < 0.05, ]  

write.table(
  kkTab_filtered,
  file = paste0(target_gene, ".result.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)	
if (nrow(kkTab_filtered) > 0) {
  cat(paste0("✅ ", target_gene, "显著富集通路表已保存（p<0.05）：", target_gene, ".result.txt\n"))
} else {
  warning(paste0("⚠️ ", target_gene, "无显著富集通路（p<0.05），输出空文件！"))
}

# ===================== 调试信息 =====================
cat("\n===== 绘图条件调试信息 =====\n")
cat("显著富集通路总数（p<0.05）：", nrow(kkTab_filtered), "\n")  
kkUp <- kkTab_filtered[kkTab_filtered$NES > 0, ]
kkDown <- kkTab_filtered[kkTab_filtered$NES < 0, ]
cat("高表达组（NES>0）显著通路数：", nrow(kkUp), "\n")
cat("低表达组（NES<0）显著通路数：", nrow(kkDown), "\n")

# ===================== 绘制高表达组图（强制绘图） =====================
if (nrow(kkUp) >= 1) {  
  showTerm <- kkUp$ID[1:min(10, nrow(kkUp))]  
  # 生成绘图对象
  p <- gseaplot2(kk, geneSetID = showTerm, base_size = 8, 
                 title = paste0("Enriched in high ", target_gene, " group (", nrow(kkUp), " pathways)"),
                 pvalue_table = TRUE, color = "red")
  
  # 1. 导出预览PNG
  png(file = "highExp_preview.png", width = 1600, height = 1000, res = 100)
  print(p)
  dev.off()
  
  # 2. 导出PDF
  pdf(file = paste0(target_gene, ".highExp.pdf"), width = 15, height = 12, onefile = FALSE)
  print(p)
  dev.off()
  
  cat("✅ 高表达组PDF重新生成完成（宽度增大至12英寸，右侧文字完整显示）\n")
}
# 重新绘制低表达组图（增大宽度，完整显示右侧表格）
if (nrow(kkDown) >= 1) {  
  showTerm <- kkDown$ID[1:min(7, nrow(kkDown))]
  p <- gseaplot2(
    kk, 
    geneSetID = showTerm, 
    base_size = 8,  
    title = paste0("Enriched in low ", target_gene, " group (", nrow(kkDown), " pathways)"),
    pvalue_table = TRUE,
    color = "blue"
  )
  
  # 1. 导出预览PNG
  png(file = "lowExp_preview.png", width = 1600, height = 700, res = 100)
  print(p)
  dev.off()
  
  # 2. 导出PDF
  pdf(file = paste0(target_gene, ".lowExp.pdf"), width = 15, height = 10, onefile = FALSE)
  print(p)
  dev.off()
  
  cat("✅ 低表达组图已重新生成（宽度增大至10英寸，右侧文字完整显示）\n")
}