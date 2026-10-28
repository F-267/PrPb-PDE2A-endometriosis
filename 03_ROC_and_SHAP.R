# ==============================================================================
# 终极修复版：全模型训练+ROC绘制+SHAP分析（含PLS/XGBoost适配）
# 核心：修复PLS/XGBoost训练逻辑，确保所有模型计入ROC曲线
# ==============================================================================

# 1. 安装包（清华镜像）
install.packages(
  c("caret", "DALEX", "ggplot2", "randomForest", "kernlab", 
    "kernelshap", "pROC", "shapviz", "xgboost", "klaR", "pls", "gbm", "nnet"),
  repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
  type = "binary",
  dependencies = TRUE,
  timeout = 120
)

# 2. 加载包
required_packages <- c("caret", "DALEX", "ggplot2", "randomForest", 
                       "kernlab", "kernelshap", "pROC", "shapviz", 
                       "xgboost", "klaR", "pls", "gbm", "nnet")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("包", pkg, "未安装！"))
  }
  library(pkg, character.only = TRUE)
}

# 3. 基础配置
setwd("YOUR_WORKING_DIR")  # 改成你存放数据的目录
set.seed(10086)
inputFile <- "merge.normalize.txt"
geneFile <- "interGenes.txt"
methodFile <- "refer.methodLists.txt"

# 4. 读取表达数据
expr_data <- read.table(
  inputFile, header = T, sep = "\t", check.names = F, row.names = 1,
  skipNul = TRUE, quote = ""
)
cat("✅ 表达数据：", nrow(expr_data), "基因 ×", ncol(expr_data), "样本\n")

# 5. 筛选目标基因
geneRT <- read.table(geneFile, header = F, sep = "\t", check.names = F, stringsAsFactors = F)
target_genes <- toupper(as.vector(geneRT[,1]))
rownames(expr_data) <- toupper(rownames(expr_data))
keep_genes <- intersect(target_genes, rownames(expr_data))
expr_filtered <- expr_data[keep_genes, , drop = FALSE]
rownames(expr_filtered) <- gsub("-", "_", rownames(expr_filtered))
cat("✅ 筛选后：", nrow(expr_filtered), "基因 ×", ncol(expr_filtered), "样本\n")

# 6. 核心修复：强制合并为Control/Treat二分类
expr_t <- t(expr_filtered)
sample_names <- row.names(expr_t)
group_raw <- sapply(strsplit(sample_names, " "), function(x) x[length(x)])
group_binary <- ifelse(
  grepl("Control", group_raw, ignore.case = TRUE),
  "Control", "Treat"
)
data_df <- as.data.frame(expr_t)
data_df$Type <- factor(group_binary, levels = c("Control", "Treat"), ordered = FALSE)

# 验证二分类
cat("\n✅ 最终二分类分布：\n")
print(table(data_df$Type))
if (length(unique(data_df$Type)) != 2) {
  stop("❌ 未成功合并为二分类！")
}

# 7. 分层抽样划分训练/测试集（确保测试集有两类）
inTrain <- createDataPartition(
  y = data_df$Type, 
  p = 0.7, 
  list = FALSE,
  times = 1,
  groups = min(table(data_df$Type))
)
train_df <- data_df[inTrain, ]
test_df <- data_df[-inTrain, ]
test_df$Type <- factor(test_df$Type, levels = c("Control", "Treat"))

# 验证划分
cat("\n✅ 训练集分布：\n")
print(table(train_df$Type))
cat("\n✅ 测试集分布：\n")
print(table(test_df$Type))
if (length(unique(test_df$Type)) != 2) {
  cat("⚠️ 调整训练集比例到50%...\n")
  inTrain <- createDataPartition(y = data_df$Type, p = 0.5, list = F)
  train_df <- data_df[inTrain, ]
  test_df <- data_df[-inTrain, ]
  test_df$Type <- factor(test_df$Type, levels = c("Control", "Treat"))
  cat("调整后测试集分布：\n")
  print(table(test_df$Type))
}

# 8. 提取测试集特征和标签
test_features <- test_df[, -ncol(test_df), drop = FALSE]
yTestClass <- test_df$Type
yTest <- ifelse(yTestClass == "Control", 0, 1)
cat("\n✅ 测试集标签类别数：", length(unique(yTest)), "\n")

# 9. 读取机器学习方法列表
methodRT <- read.table(methodFile, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
cat("✅ 读取", nrow(methodRT), "种方法\n")

# 10. 配置交叉验证（二分类专用）
control <- trainControl(
  method = "repeatedcv", 
  number = 5, 
  repeats = 3,
  savePredictions = "final",
  classProbs = TRUE,          
  summaryFunction = twoClassSummary  
)

# 11. 自定义模型参数网格（解决XGBoost/PLS训练失败）
# 11.1 XGBoost参数网格
xgb_grid <- expand.grid(
  nrounds = 50,
  max_depth = 3,
  eta = 0.1,
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)
# 11.2 PLS参数网格
pls_grid <- expand.grid(ncomp = 1:5)

# 12. 循环训练模型+绘制ROC（含PLS/XGBoost适配）
pdf(file = "ROC.pdf", width = 5.5, height = 5)
modelList <- list()
AUC <- c()
ROCcolor <- rainbow(nrow(methodRT))
plot_initialized <- FALSE

for (i in 1:nrow(methodRT)) {
  name <- methodRT[i, "Name"]
  method <- methodRT[i, "Method"]
  cat("\n🔹 训练模型：", name, "(", method, ")\n")
  
  # 适配不同模型的训练逻辑
  model <- tryCatch({
    switch(
      name,
      # PLS模型适配
      "PLS" = train(
        Type ~ ., 
        data = train_df, 
        method = "pls", 
        trControl = control,
        metric = "ROC",
        tuneGrid = pls_grid,  # 自定义PLS参数
        preProcess = c("center", "scale")  # PLS必须标准化
      ),
      # XGBoost模型适配（标签转0/1）
      "XGBoost" = {
        train_df_xgb <- train_df
        train_df_xgb$Type <- ifelse(train_df_xgb$Type == "Control", 0, 1)
        train(
          Type ~ ., 
          data = train_df_xgb, 
          method = "xgbTree", 
          trControl = control,
          metric = "ROC",
          tuneGrid = xgb_grid,  # 自定义XGBoost参数
          verbose = FALSE
        )
      },
      # SVM模型适配
      "SVM" = train(
        Type ~ ., 
        data = train_df, 
        method = "svmRadial", 
        prob.model = TRUE,
        trControl = control,
        metric = "ROC",
        preProcess = c("center", "scale")
      ),
      # 其他模型通用训练
      train(
        Type ~ ., 
        data = train_df, 
        method = method, 
        trControl = control,
        metric = "ROC"
      )
    )
  }, error = function(e) {
    cat("❌", name, "训练失败：", substr(e$message, 1, 50), "\n")
    return(NULL)
  })
  
  if (is.null(model)) next
  
  # 预测概率（适配不同模型的输出格式）
  pred <- tryCatch({
    # XGBoost预测后转换概率格式
    if (name == "XGBoost") {
      pred_prob <- predict(model, newdata = test_features, type = "prob")
      colnames(pred_prob) <- c("Control", "Treat")  # 强制列名
    } else {
      pred_prob <- predict(model, newdata = test_features, type = "prob")
    }
    # 统一概率列名
    if (all(colnames(pred_prob) %in% c("Control", "Treat"))) {
      pred_prob
    } else {
      colnames(pred_prob) <- c("Control", "Treat")
      pred_prob
    }
  }, error = function(e) {
    cat("❌", name, "预测失败：", substr(e$message, 1, 50), "\n")
    return(NULL)
  })
  
  if (is.null(pred) || ncol(pred) != 2) {
    cat("❌", name, "概率列异常，跳过\n")
    next
  }
  
  # 计算ROC（统一逻辑）
  roc_obj <- tryCatch({
    roc(
      response = yTestClass,
      predictor = pred[, "Treat"],
      levels = c("Control", "Treat"),
      direction = "<"
    )
  }, error = function(e) {
    cat("❌", name, "ROC失败：", substr(e$message, 1, 50), "\n")
    return(NULL)
  })
  
  if (is.null(roc_obj)) next
  
  # 绘制ROC曲线
  auc_val <- sprintf("%.03f", roc_obj$auc)
  AUC <- c(AUC, paste0(name, ': ', auc_val))
  modelList[[method]] <- roc_obj$auc
  
  if (!plot_initialized) {
    plot(roc_obj, print.auc = F, legacy.axes = T, main = "", col = ROCcolor[i], lwd = 3)
    plot_initialized <- TRUE
  } else {
    plot(roc_obj, print.auc = F, legacy.axes = T, main = "", col = ROCcolor[i], lwd = 3, add = T)
  }
}

# 添加图例
if (plot_initialized) {
  legend('bottomright', AUC, col = ROCcolor[1:length(AUC)], lwd = 3, bty = 'n', cex = 0.9)
}
dev.off()

# 13. 输出未计入的模型（便于排查）
all_models <- methodRT$Name
plot_models <- sapply(strsplit(AUC, ":"), function(x) x[1])
missing_models <- setdiff(all_models, plot_models)
cat("\n📋 方法列表所有模型：", paste(all_models, collapse = ", "), "\n")
cat("✅ 已计入ROC的模型：", paste(plot_models, collapse = ", "), "\n")
if (length(missing_models) > 0) {
  cat("❌ 未计入ROC的模型：", paste(missing_models, collapse = ", "), "\n")
} else {
  cat("🎉 所有模型均成功计入ROC曲线！\n")
}

# 14. 选择最优模型+SHAP分析（适配SVM模型）
if (length(modelList) > 0) {
  aucValue <- unlist(modelList)
  bestMethod <- names(which(aucValue == max(aucValue)))
  bestName <- methodRT[methodRT$Method == bestMethod, "Name"]
  cat("\n🏆 最优模型：", bestName, "(", bestMethod, ") | AUC：", max(aucValue), "\n")
  
  # 重新训练最优模型（统一格式）
  train_df_best <- train_df
  if (bestMethod == "xgbTree") {
    train_df_best$Type <- ifelse(train_df_best$Type == "Control", 0, 1)
  }
  best_model <- train(
    Type ~ ., 
    data = train_df_best, 
    method = bestMethod, 
    trControl = control,
    metric = "ROC"
  )
  
  # ========== 核心修复：适配SVM的SHAP计算（强制预测为数值型） ==========
  train_features <- train_df_best[, -ncol(train_df_best)]
  
  # 定义自定义预测函数（强制输出数值型）
  custom_pred_fun <- function(model, newdata) {
    # 1. 若为分类模型，输出正类概率（数值型）
    if (bestMethod %in% c("svmRadial", "rf", "glm", "knn", "pls")) {
      pred_prob <- predict(model, newdata = newdata, type = "prob")
      # 强制返回Treat类的概率（数值型向量）
      return(as.numeric(pred_prob[, "Treat"]))
    } 
    # 2. 若为XGBoost（数值标签），直接返回预测值
    else if (bestMethod == "xgbTree") {
      pred <- predict(model, newdata = newdata)
      return(as.numeric(pred))
    }
    # 3. 其他模型默认返回数值型预测
    else {
      pred <- predict(model, newdata = newdata)
      return(as.numeric(pred))
    }
  }
  
  # 计算SHAP值（使用自定义预测函数）
  cat("⚠️ 计算SHAP值（适配", bestName, "模型）\n")
  if (ncol(train_features) > 20) {
    # kernelshap适配
    fit <- kernelshap(
      best_model, 
      X = train_features,
      pred_fun = custom_pred_fun  # 使用自定义预测函数
    )
  } else {
    # permshap适配
    fit <- permshap(
      best_model, 
      X = train_features,
      pred_fun = custom_pred_fun  # 使用自定义预测函数
    )
  }
  
  # SHAP可视化（5类图，兼容所有模型）
  shp <- shapviz(fit, X_pred = train_features, X = train_features, interactions = TRUE)
  
  # 14.1 特征重要性柱状图
  pdf(file = "barplot.pdf", width = 6, height = 6)
  print(sv_importance(shp, kind = "bar", show_numbers = TRUE) + theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
  
  # 14.2 蜂群图
  pdf(file = "bee.pdf", width = 7, height = 6)
  print(sv_importance(shp, kind = "bee", show_numbers = TRUE) + theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
  
  # 14.3 依赖图
  pdf(file = "dependence.pdf", width = 9, height = 6)
  print(sv_dependence(shp, v = names(sort(colMeans(abs(shp$S)), TRUE))) + theme_bw())
  dev.off()
  
  # 14.4 瀑布图
  pdf(file = "waterfall.pdf", width = 7, height = 5)
  print(sv_waterfall(shp, row_id = 1) + theme_bw())
  dev.off()
  
  # 14.5 力图
  pdf(file = "force.pdf", width = 9, height = 5)
  print(sv_force(shp, row_id = 1))
  dev.off()
  
  cat("\n🎉 分析完成！生成：\n")
  cat("  - ROC.pdf（所有模型ROC曲线）\n")
  cat("  - barplot.pdf/bee.pdf/dependence.pdf/waterfall.pdf/force.pdf（SHAP可视化）\n")
} else {
  cat("\n❌ 无模型训练成功，请检查方法列表或数据！\n")
}