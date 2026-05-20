# ============================================================
# M5 数据下载脚本 v5 — 直接调用 GDC REST API 下载
# 完全绕开 TCGAbiolinks 的下载/解压逻辑
# 在 RStudio Console 全选运行
# ============================================================

library(TCGAbiolinks)   # 只用它查询文件 ID 列表
library(httr)           # 直接 HTTP 下载
library(SummarizedExperiment)

GDCdata_dir <- "D:/Bio-Informatics Case Study/M5_Survival/GDC_data"
dir.create(GDCdata_dir, showWarnings = FALSE, recursive = TRUE)

# ── Step 1：查询获取文件 ID 列表 ──
cat("🔍 查询 TCGA-LAML 文件列表...\n")
query_exp <- GDCquery(
  project               = "TCGA-LAML",
  data.category         = "Transcriptome Profiling",
  data.type             = "Gene Expression Quantification",
  workflow.type         = "STAR - Counts",
  experimental.strategy = "RNA-Seq"
)

results   <- getResults(query_exp)
file_ids  <- results$id
file_names <- results$file_name
n_total   <- length(file_ids)
cat("✅ 共", n_total, "个文件\n\n")

# ── Step 2：逐个下载（直接 GET，R 内置解压，不用 tar.exe）──
cat("⬇️  开始逐个下载文件...\n")
cat("   每个文件约 4 MB，合计约 600 MB\n")
cat("   已有文件自动跳过\n\n")

success <- 0
failed  <- c()

for (i in seq_along(file_ids)) {
  fid   <- file_ids[i]
  fname <- file_names[i]

  # 每个文件存到以 file_id 命名的子目录（TCGAbiolinks 的标准结构）
  out_subdir <- file.path(GDCdata_dir, "TCGA-LAML", fid)
  out_file   <- file.path(out_subdir, fname)

  # 已下载则跳过
  if (file.exists(out_file) && file.size(out_file) > 1000) {
    success <- success + 1
    if (i %% 20 == 0) cat(sprintf("  [%d/%d] %d 个完成，跳过已有文件...\n", i, n_total, success))
    next
  }

  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)

  # GDC 单文件下载 URL
  url <- paste0("https://api.gdc.cancer.gov/data/", fid)

  tryCatch({
    resp <- GET(
      url,
      write_disk(out_file, overwrite = TRUE),
      timeout(120),
      progress()
    )

    if (status_code(resp) == 200 && file.size(out_file) > 1000) {
      success <- success + 1
    } else {
      file.remove(out_file)
      failed <- c(failed, fid)
      cat(sprintf("  ⚠️  [%d] 下载异常（状态码 %d）: %s\n", i, status_code(resp), fname))
    }
  }, error = function(e) {
    if (file.exists(out_file)) file.remove(out_file)
    failed <<- c(failed, fid)
    cat(sprintf("  ❌ [%d] 错误: %s\n", i, conditionMessage(e)))
  })

  # 每 10 个报一次进度
  if (i %% 10 == 0) {
    cat(sprintf("  [%d/%d] 成功: %d  失败: %d\n", i, n_total, success, length(failed)))
  }
}

cat(sprintf("\n📊 下载结果：成功 %d / %d，失败 %d\n", success, n_total, length(failed)))

if (length(failed) > 0) {
  cat("失败的文件 ID（可重新运行脚本补下）:\n")
  cat(paste(" ", head(failed, 5)), sep = "\n")
}

# ── Step 3：GDCprepare 整理 ──
if (success >= n_total * 0.9) {
  cat("\n📦 整理数据（GDCprepare）...\n")

  se_exp <- tryCatch({
    GDCprepare(
      query     = query_exp,
      directory = GDCdata_dir,
      save      = TRUE,
      save.filename = file.path(GDCdata_dir, "TCGA_LAML_expr.rda")
    )
  }, error = function(e) {
    cat("❌ GDCprepare 错误:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(se_exp)) {
    cat("\n🎉 全部完成！\n")
    cat("  基因数:", nrow(se_exp), "\n")
    cat("  样本数:", ncol(se_exp), "\n")
    cat("\n✅ 现在可以 Knit M5_Survival_Analysis.Rmd 了\n")
  }
} else {
  cat("\n⚠️  成功率不足 90%，建议重新运行本脚本补下失败文件\n")
}
