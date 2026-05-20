# ============================================================
# M5 生存分析 — R 包安装脚本
# 运行方式：在 RStudio 里打开此文件，全选 Ctrl+A，运行 Ctrl+Enter
# 预计时间：首次安装约 10-20 分钟（国内建议挂代理或使用清华镜像）
# ============================================================

# 1. 设置清华镜像（国内加速）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# 2. 安装 BiocManager（Bioconductor 包管理器）
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 3. 安装核心 Bioconductor 包
#    TCGAbiolinks：从 GDC 数据库自动下载 TCGA 数据
#    SummarizedExperiment：储存基因表达矩阵的标准容器
BiocManager::install(c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "BiocParallel"
), ask = FALSE, update = FALSE)

# 4. 安装 CRAN 包
#    survival：经典生存分析（K-M 曲线、Cox 回归）
#    survminer：ggplot2 风格的生存曲线可视化
#    tidyverse：数据整理（dplyr、tidyr、ggplot2 等）
#    ggpubr：排版多图
#    forestplot：绘制 Cox 森林图
packages_cran <- c(
  "survival",
  "survminer",
  "tidyverse",
  "ggpubr",
  "dplyr",
  "tibble",
  "ggplot2",
  "scales"
)

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
    message(paste("已安装:", pkg))
  } else {
    message(paste("已存在，跳过:", pkg))
  }
}

# 5. 验证安装
cat("\n========== 安装验证 ==========\n")
required <- c("TCGAbiolinks", "survival", "survminer", "tidyverse", "SummarizedExperiment")
for (pkg in required) {
  status <- ifelse(requireNamespace(pkg, quietly = TRUE), "✅ 已安装", "❌ 未安装")
  cat(pkg, ":", status, "\n")
}
cat("================================\n")
cat("全部 ✅ 后，运行 M5_Survival_Analysis.Rmd\n")
