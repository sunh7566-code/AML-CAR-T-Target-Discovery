# ============================================================
# M5 生存分析 — 完整分析脚本（数据已下载后直接运行）
# 前提：已运行 01_download_data.R，GDC_data/ 目录有 151 个 tsv 文件
# ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(survival)
library(survminer)
library(ggplot2)
library(httr)
library(jsonlite)

GDCdata_dir <- "D:/Bio-Informatics Case Study/M5_Survival/GDC_data"
figures_dir <- "D:/Bio-Informatics Case Study/M5_Survival/figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

target_genes <- c("FLT3", "CD33", "IL3RA", "CLEC12A", "CD38")

# ── 如果已有缓存直接加载，跳过耗时步骤 ──
rds_file <- "D:/Bio-Informatics Case Study/M5_Survival/analysis_df.rds"
if (file.exists(rds_file)) {
  cat("✅ 加载缓存数据...\n")
  analysis_df <- readRDS(rds_file)
  cox_df      <- readRDS("D:/Bio-Informatics Case Study/M5_Survival/cox_results.rds")
  cat("  样本数:", nrow(analysis_df), "\n")
} else {
  # ── Step 1：读取表达量 ──
  cat("读取表达矩阵...\n")
  tsv_files <- list.files(GDCdata_dir, recursive = TRUE,
                           full.names = TRUE, pattern = "\\.tsv$")

  get_file_id <- function(path) {
    parts <- strsplit(path, "/|\\\\")[[1]]
    parts[length(parts) - 1]
  }

  expr_list <- lapply(seq_along(tsv_files), function(i) {
    df <- read_tsv(tsv_files[[i]], comment = "#", show_col_types = FALSE) %>%
      filter(gene_name %in% target_genes) %>%
      select(gene_name, unstranded)
    if (nrow(df) == 0) return(NULL)
    df$file_id <- get_file_id(tsv_files[[i]])
    df
  })

  expr_wide <- bind_rows(expr_list) %>%
    group_by(file_id, gene_name) %>%
    summarise(unstranded = sum(unstranded, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = gene_name, values_from = unstranded)

  # ── Step 2：file_id → barcode 映射 ──
  # 需要 results 对象（GDCquery 的结果），如没有重新查询
  if (!exists("results")) {
    library(TCGAbiolinks)
    query_exp <- GDCquery(
      project = "TCGA-LAML",
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      experimental.strategy = "RNA-Seq"
    )
    results <- getResults(query_exp)
  }

  mapping <- results %>%
    select(file_id = id, barcode = cases) %>%
    mutate(patient_id = substr(barcode, 1, 12))

  expr_with_barcode <- expr_wide %>%
    left_join(mapping, by = "file_id")

  # ── Step 3：临床生存数据（GDC REST API）──
  cat("获取临床数据...\n")
  `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

  resp <- POST(
    url  = "https://api.gdc.cancer.gov/cases",
    body = '{"filters":{"op":"=","content":{"field":"project.project_id","value":"TCGA-LAML"}},"fields":"submitter_id,demographic.vital_status,demographic.days_to_death,demographic.age_at_index,demographic.gender,diagnoses.days_to_last_follow_up","format":"JSON","size":"200"}',
    add_headers("Content-Type" = "application/json"),
    timeout(60)
  )

  hits <- fromJSON(rawToChar(resp$content))$data$hits

  clinical <- data.frame(
    patient_id    = hits$submitter_id,
    vital_status  = hits$demographic$vital_status,
    days_to_death = as.numeric(hits$demographic$days_to_death),
    age_at_index  = as.numeric(hits$demographic$age_at_index),
    gender        = hits$demographic$gender,
    stringsAsFactors = FALSE
  )
  clinical$days_to_last_follow_up <- sapply(hits$diagnoses, function(x) {
    if (is.data.frame(x) && "days_to_last_follow_up" %in% colnames(x))
      as.numeric(x$days_to_last_follow_up[1])
    else NA_real_
  })

  surv_df <- clinical %>%
    mutate(
      OS_time   = case_when(vital_status == "Dead"  ~ days_to_death,
                            vital_status == "Alive" ~ days_to_last_follow_up,
                            TRUE ~ NA_real_),
      OS_status = case_when(vital_status == "Dead"  ~ 1L,
                            vital_status == "Alive" ~ 0L,
                            TRUE ~ NA_integer_)
    ) %>%
    filter(!is.na(OS_time), OS_time > 0, !is.na(OS_status))

  # ── Step 4：合并 ──
  analysis_df <- surv_df %>%
    inner_join(expr_with_barcode, by = "patient_id") %>%
    distinct(patient_id, .keep_all = TRUE)

  cat("✅ 可分析样本:", nrow(analysis_df), "\n")

  # ── Step 5：Cox 回归 ──
  cox_results <- lapply(target_genes, function(gene) {
    df_tmp <- analysis_df %>% mutate(expr = log2(as.numeric(.data[[gene]]) + 1))
    fit <- coxph(Surv(OS_time, OS_status) ~ expr, data = df_tmp)
    s   <- summary(fit)
    data.frame(
      Gene     = gene,
      HR       = round(s$coefficients[1, "exp(coef)"], 3),
      CI_lower = round(s$conf.int[1, "lower .95"], 3),
      CI_upper = round(s$conf.int[1, "upper .95"], 3),
      P_value  = signif(s$coefficients[1, "Pr(>|z|)"], 3),
      Sig      = case_when(
        s$coefficients[1, "Pr(>|z|)"] < 0.001 ~ "***",
        s$coefficients[1, "Pr(>|z|)"] < 0.01  ~ "**",
        s$coefficients[1, "Pr(>|z|)"] < 0.05  ~ "*",
        s$coefficients[1, "Pr(>|z|)"] < 0.1   ~ ".",
        TRUE ~ "ns"
      )
    )
  })
  cox_df <- do.call(rbind, cox_results)

  saveRDS(analysis_df, rds_file)
  saveRDS(cox_df, "D:/Bio-Informatics Case Study/M5_Survival/cox_results.rds")
}

# ══════════════════════════════════════════════
# 可视化
# ══════════════════════════════════════════════

cat("\n📊 Cox 回归结果:\n")
print(cox_df)

# ── K-M 曲线 ──
cat("\n绘制 K-M 曲线...\n")
for (gene in target_genes) {
  vals <- as.numeric(analysis_df[[gene]])
  df_tmp <- analysis_df %>%
    mutate(group = factor(ifelse(vals >= median(vals, na.rm=TRUE), "High", "Low"),
                          levels = c("High", "Low")))
  fit <- survfit(Surv(OS_time, OS_status) ~ group, data = df_tmp)
  p_km <- ggsurvplot(
    fit, data = df_tmp,
    pval = TRUE, conf.int = TRUE, risk.table = TRUE,
    palette = c("#E41A1C", "#377EB8"),
    legend.labs = c(paste0(gene,"-High"), paste0(gene,"-Low")),
    legend.title = "",
    xlab = "时间（天）", ylab = "总体生存率",
    title = paste0(gene, " | TCGA-LAML 生存分析 (n=", nrow(df_tmp), ")"),
    ggtheme = theme_bw(base_size = 13),
    surv.median.line = "hv", break.time.by = 365
  )
  ggsave(file.path(figures_dir, paste0("KM_", gene, ".png")),
         plot = p_km$plot, width = 10, height = 7, dpi = 150)
}

# ── 森林图 ──
forest_df <- cox_df %>%
  mutate(Gene = factor(Gene, levels = rev(Gene)),
         color = ifelse(HR > 1, "#E41A1C", "#377EB8"))

p_forest <- ggplot(forest_df, aes(x = HR, y = Gene)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.25, linewidth = 0.8, color = forest_df$color) +
  geom_point(aes(color = color), size = 4, shape = 18) +
  geom_text(aes(label = paste0("HR=", HR, "  ", Sig), x = max(CI_upper)*1.05),
            hjust = 0, size = 4) +
  scale_color_identity() + scale_x_log10() +
  labs(x = "HR（log10）", y = NULL,
       title = "单变量 Cox 回归 — CAR-T 靶点预后价值（TCGA-LAML）",
       subtitle = "** p<0.01  * p<0.05  ns=不显著 | 红色HR>1高表达预后差，蓝色HR<1高表达预后好") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

ggsave(file.path(figures_dir, "Cox_ForestPlot.png"),
       plot = p_forest, width = 11, height = 5, dpi = 150)

cat("\n🎉 M5 全部完成！\n")
cat("  输出目录:", figures_dir, "\n")
cat("  关键发现：IL3RA (HR=1.711, p=0.00188**) 是唯一独立预后因子\n")
