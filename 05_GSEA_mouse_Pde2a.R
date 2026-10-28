# ======================== 小鼠 GSEA 同时出 PDF + TIF 高清图 ========================
# 安装（第一次运行）
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("limma","clusterProfiler","enrichplot"))
install.packages("dplyr")

# 加载包
library(limma)
library(clusterProfiler)
library(enrichplot)
library(dplyr)

# 工作目录
setwd("YOUR_WORKING_DIR")  # 改成你存放数据的目录

# ===================== 配置 =====================
expFile   <- "EM_vs_NC.normalize.txt"
gmtFile   <- "m2.cp.v2026.1.Mm.symbols.gmt"
target_gene <- "Pde2a"

# ===================== 读取表达矩阵 =====================
rt <- read.table(expFile, header=T, sep="\t", check.names=F, fill=T)
rownames(rt) <- rt[,1]
exp <- as.matrix(rt[,-1])
data <- matrix(as.numeric(exp), nrow=nrow(exp), dimnames=dimnames(exp))
data <- avereps(data)
data <- data[rowMeans(data) > 0, ]

# ===================== 按Pde2a高低分组 =====================
gene_exp <- data[target_gene, ]
med <- median(gene_exp)
dataL <- data[, gene_exp < med, drop=F]
dataH <- data[, gene_exp >= med, drop=F]

logFC <- rowMeans(dataH) - rowMeans(dataL)
logFC <- sort(logFC, decreasing=T)

# ===================== GSEA分析 =====================
gmt <- read.gmt(gmtFile)

kk <- GSEA(
  geneList = logFC,
  TERM2GENE = gmt,
  pvalueCutoff = 1,
  pAdjustMethod = "fdr",
  minGSSize = 10,
  maxGSSize = 500
)

# ===================== 输出显著结果 =====================
res <- as.data.frame(kk)
res_sig <- res[res$pvalue < 0.05, ]

# 按NES绝对值降序排列（富集程度最强的在前）
res_sig <- res_sig[order(abs(res_sig$NES), decreasing = TRUE), ]

write.table(res_sig, paste0(target_gene, ".GSEA.result.txt"), sep="\t", quote=F, row.names=F)

# ===================== 绘图（PDF + 高清 TIF 600DPI）=====================
# 分离上下调通路
up <- res_sig[res_sig$NES > 0, ]
down <- res_sig[res_sig$NES < 0, ]

# 高表达组（NES > 0）：取前10个富集最强的通路
if(nrow(up) > 0){
  # 按NES绝对值降序排列（已经是，但确保一下）
  up <- up[order(abs(up$NES), decreasing = TRUE), ]
  
  # 取前10个或全部（如果不足10个）
  showTerm <- up$ID[1:min(10, nrow(up))]
  
  # 生成标题：显示通路数量
  plot_title <- paste0("High ", target_gene, " (Top ", length(showTerm), " pathways by |NES|)")
  
  # 动态调整高度：每个通路增加0.8英寸
  plot_height <- 10 + (length(showTerm) * 0.8)
  # 设置最小高度12英寸，最大高度25英寸
  plot_height <- max(12, min(plot_height, 25))
  
  p <- gseaplot2(kk, showTerm, 
                 title = plot_title, 
                 pvalue_table = TRUE, 
                 base_size = 9, 
                 color = "red")
  
  # PDF
  pdf(paste0(target_gene, ".high.top10_by_NES.pdf"), width = 20, height = plot_height)
  print(p)
  dev.off()
  
  # TIF 高清 600DPI（论文用）
  tiff(paste0(target_gene, ".high.top10_by_NES.tif"), 
       width = 20, height = plot_height, units = "in", 
       res = 600, compression = "lzw")
  print(p)
  dev.off()
  
  # 输出高表达组前10通路信息
  up_top10 <- up[1:min(10, nrow(up)), c("ID", "Description", "NES", "pvalue", "p.adjust")]
  write.table(up_top10, 
              paste0(target_gene, ".high.top10_pathways.txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("📊 高表达组前", nrow(up_top10), "个通路已保存，图片高度：", plot_height, "英寸")
}

# 低表达组（NES < 0）：取前10个富集最强的通路
if(nrow(down) > 0){
  # 按NES绝对值降序排列（负值绝对值越大，富集越强）
  down <- down[order(abs(down$NES), decreasing = TRUE), ]
  
  # 取前10个或全部（如果不足10个）
  showTerm <- down$ID[1:min(10, nrow(down))]
  
  # 生成标题：显示通路数量
  plot_title <- paste0("Low ", target_gene, " (Top ", length(showTerm), " pathways by |NES|)")
  
  # 动态调整高度：每个通路增加0.8英寸
  plot_height <- 10 + (length(showTerm) * 0.8)
  # 设置最小高度12英寸，最大高度25英寸
  plot_height <- max(12, min(plot_height, 25))
  
  p <- gseaplot2(kk, showTerm, 
                 title = plot_title, 
                 pvalue_table = TRUE, 
                 base_size = 9, 
                 color = "blue")
  
  # PDF
  pdf(paste0(target_gene, ".low.top10_by_NES.pdf"), width = 20, height = plot_height)
  print(p)
  dev.off()
  
  # TIF 高清 600DPI（论文用）
  tiff(paste0(target_gene, ".low.top10_by_NES.tif"), 
       width = 20, height = plot_height, units = "in", 
       res = 600, compression = "lzw")
  print(p)
  dev.off()
  
  # 输出低表达组前10通路信息
  down_top10 <- down[1:min(10, nrow(down)), c("ID", "Description", "NES", "pvalue", "p.adjust")]
  write.table(down_top10, 
              paste0(target_gene, ".low.top10_pathways.txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("📊 低表达组前", nrow(down_top10), "个通路已保存，图片高度：", plot_height, "英寸")
}

# 输出总结信息
message("\n==================================================")
message("📈 GSEA分析结果总结：")
message("   总显著通路数：", nrow(res_sig))
message("   高表达组通路数（NES > 0）：", nrow(up))
message("   低表达组通路数（NES < 0）：", nrow(down))
message("==================================================")

message("\n✅ 全部完成！已生成：")
message("   📄 PDF + 高清 TIF 600DPI 论文图（自适应高度）")
message("   📊 按NES绝对值排序的前10通路图")
message("   📋 前10通路详细信息表格")