# ============================================================
# M7 免疫微环境分析（Tumor Microenvironment, TME）
# 项目：AML CAR-T 靶点发现
# 数据：TCGA-LAML（151个样本，直接复用 M5 已下载数据）
# 工具：IOBR 2.2.1（整合 ssGSEA / CIBERSORT / TIMER 等算法）
# 输出：免疫细胞比例热图 + 靶点相关性散点图
# ============================================================

# ── 0. 加载包 ─────────────────────────────────────────────
# 原因：IOBR 提供统一的反卷积接口；tidyverse 负责数据整理
library(IOBR)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)  # ComplexHeatmap 配色需要

cat("✅ 包加载完成\n")

# ── 1. 读取 TPM 矩阵 ─────────────────────────────────────
# 原因：IOBR 要求输入 TPM（而非 raw count）
#       因为不同样本的测序深度不同，TPM 已归一化，才能横向比较
#       GDC 的 augmented_star_gene_counts.tsv 文件里直接有 tpm_unstranded 列，省去转换步骤

cat("📂 读取 TCGA-LAML TPM 数据...\n")

data_dir <- "D:/Bio-Informatics Case Study/M5_Survival/GDC_data/TCGA-LAML"
sample_dirs <- list.dirs(data_dir, recursive = FALSE)

# 读取每个样本的 TPM 值（tpm_unstranded 列）
read_tpm <- function(sample_dir) {
  tsv_file <- list.files(sample_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (length(tsv_file) == 0) return(NULL)

  df <- read_tsv(tsv_file, skip = 1, show_col_types = FALSE)
  # 过滤掉 N_unmapped 等质控行（gene_id 不以 ENSG 开头的）
  df <- df %>% filter(grepl("^ENSG", gene_id))
  # 只保留 gene_name（基因符号）和 tpm_unstranded
  df %>% select(gene_name, tpm_unstranded)
}

# 批量读取所有样本，构建宽格式矩阵（行=基因，列=样本）
tpm_list <- lapply(sample_dirs, read_tpm)

# 用 sample_dir 文件夹名作为临时样本 ID
names(tpm_list) <- basename(sample_dirs)

# 合并：每个样本是一列
tpm_wide <- tpm_list %>%
  bind_rows(.id = "sample_id") %>%
  pivot_wider(names_from = sample_id, values_from = tpm_unstranded, values_fn = mean)

# 设置行名为基因符号
tpm_mat <- as.data.frame(tpm_wide)
rownames(tpm_mat) <- tpm_mat$gene_name
tpm_mat$gene_name <- NULL

# 去掉有 NA 的基因（极少数重复基因名导致）
tpm_mat <- tpm_mat[complete.cases(tpm_mat), ]

cat(sprintf("✅ TPM 矩阵维度：%d 基因 × %d 样本\n", nrow(tpm_mat), ncol(tpm_mat)))

# 保存缓存（下次直接加载，避免重新读 151 个文件）
saveRDS(tpm_mat, "D:/Bio-Informatics Case Study/M7_TME/tpm_matrix.rds")
cat("💾 TPM 矩阵已缓存到 M7_TME/tpm_matrix.rds\n")

# ── 2. 获取真实 TCGA 样本 barcode ─────────────────────────
# 原因：目前列名是文件夹的 UUID，需要转换为 TCGA barcode（如 TCGA-AB-2802）
#       才能和临床数据合并

cat("🔗 获取 GDC 样本 UUID → TCGA barcode 映射...\n")

# 通过 GDC API 查询文件信息
library(httr)
library(jsonlite)

# 批量查询（每次最多 100 个 file_id）
get_barcode_map <- function(file_ids) {
  body <- list(
    filters = list(
      op = "in",
      content = list(
        field = "file_id",
        value = as.list(file_ids)
      )
    ),
    fields = "file_id,cases.submitter_id",
    format = "JSON",
    size = length(file_ids)
  )

  resp <- POST(
    "https://api.gdc.cancer.gov/files",
    body = toJSON(body, auto_unbox = TRUE),
    content_type_json(),
    timeout(30)
  )

  if (status_code(resp) != 200) return(NULL)

  result <- fromJSON(content(resp, "text", encoding = "UTF-8"))
  hits <- result$data$hits

  map_df <- data.frame(
    uuid = hits$id,
    barcode = sapply(hits$cases, function(x) {
      if (is.null(x) || nrow(x) == 0) return(NA)
      x$submitter_id[1]
    }),
    stringsAsFactors = FALSE
  )
  return(map_df)
}

# 获取文件夹 UUID（这是文件夹名，不是 file_id，需要从 tsv 文件名获取 file_id）
# 实际上：文件夹名 = case UUID，tsv 文件名 = file_id
file_uuids <- character(length(sample_dirs))
for (i in seq_along(sample_dirs)) {
  tsv_file <- list.files(sample_dirs[i], pattern = "\\.tsv$")
  if (length(tsv_file) > 0) {
    file_uuids[i] <- tools::file_path_sans_ext(tsv_file[1])
    # 去掉 .rna_seq.augmented_star_gene_counts 后缀
    file_uuids[i] <- sub("\\.rna_seq.*", "", tsv_file[1])
  }
}
file_uuids <- file_uuids[file_uuids != ""]

# 分批查询
batch_size <- 100
barcode_maps <- list()
for (i in seq(1, length(file_uuids), by = batch_size)) {
  batch <- file_uuids[i:min(i + batch_size - 1, length(file_uuids))]
  bm <- get_barcode_map(batch)
  if (!is.null(bm)) barcode_maps[[length(barcode_maps) + 1]] <- bm
  Sys.sleep(0.5)
}

barcode_map <- bind_rows(barcode_maps)
cat(sprintf("✅ 获取到 %d 个 UUID→barcode 映射\n", nrow(barcode_map)))

# 重命名 tpm_mat 的列名
# 先建立文件夹UUID → file_uuid → barcode 的映射
folder_to_file <- data.frame(
  folder_uuid = basename(sample_dirs),
  file_uuid = file_uuids,
  stringsAsFactors = FALSE
)

col_map <- folder_to_file %>%
  left_join(barcode_map, by = c("file_uuid" = "uuid"))

# 重命名 tpm_mat 列
new_colnames <- col_map$barcode[match(colnames(tpm_mat), col_map$folder_uuid)]
valid_idx <- !is.na(new_colnames)
tpm_mat_named <- tpm_mat[, valid_idx]
colnames(tpm_mat_named) <- new_colnames[valid_idx]

cat(sprintf("✅ 成功重命名 %d 个样本为 TCGA barcode\n", ncol(tpm_mat_named)))

# 更新缓存
saveRDS(tpm_mat_named, "D:/Bio-Informatics Case Study/M7_TME/tpm_matrix_named.rds")

# ── 3. 运行免疫细胞反卷积（ssGSEA）──────────────────────
# 原因：ssGSEA（single-sample Gene Set Enrichment Analysis）
#       原理：用已知的免疫细胞特征基因集，给每个样本打分
#       类比：像药物筛选中用 marker 化合物的方式来"鉴定"细胞类型
#       优点：不需要参考数据集，基于基因表达即可运行

cat("🔬 运行 ssGSEA 免疫细胞反卷积...\n")
cat("（这一步可能需要 3-5 分钟，请耐心等待）\n")

# IOBR 2.x 的新接口：deconvo_tme()
# 输入必须是 log2(TPM+1) 格式
tpm_log2 <- log2(tpm_mat_named + 1)

# 运行 ssGSEA
immune_ssgsea <- deconvo_tme(
  eset = tpm_log2,
  method = "ssgsea",     # 用 ssGSEA 方法
  arrays = FALSE,        # RNA-seq 数据，不是芯片
  scale_eset = TRUE      # 标准化处理
)

cat(sprintf("✅ ssGSEA 完成，输出 %d 样本 × %d 免疫细胞指标\n",
            nrow(immune_ssgsea), ncol(immune_ssgsea)))

saveRDS(immune_ssgsea, "D:/Bio-Informatics Case Study/M7_TME/immune_ssgsea.rds")

# ── 4. 运行 TIMER 反卷积 ──────────────────────────────────
# 原因：TIMER 是专门为 TCGA 数据设计的算法，对 AML 更准确
#       与 ssGSEA 互相验证，结果更可信

cat("🔬 运行 TIMER 免疫细胞反卷积...\n")

immune_timer <- deconvo_tme(
  eset = tpm_log2,
  method = "timer",
  cancer_type = "LAML"   # 指定癌症类型：急性髓系白血病
)

cat(sprintf("✅ TIMER 完成，输出 %d 样本 × %d 免疫细胞指标\n",
            nrow(immune_timer), ncol(immune_timer)))

saveRDS(immune_timer, "D:/Bio-Informatics Case Study/M7_TME/immune_timer.rds")

# ── 5. 整合靶点表达量 ─────────────────────────────────────
# 我们的 5 个候选 CAR-T 靶点
car_t_targets <- c("FLT3", "CD33", "IL3RA", "CLEC12A", "CD38")

# 提取靶点的 TPM 表达量
target_expr <- tpm_log2[car_t_targets[car_t_targets %in% rownames(tpm_log2)], ]
target_df <- as.data.frame(t(target_expr))
target_df$sample_id <- rownames(target_df)

# 整合 ssGSEA 结果
immune_df <- immune_ssgsea
immune_df$sample_id <- rownames(immune_df)

combined_df <- inner_join(target_df, immune_df, by = "sample_id")
cat(sprintf("✅ 整合完成，共 %d 个样本可用于相关性分析\n", nrow(combined_df)))

saveRDS(combined_df, "D:/Bio-Informatics Case Study/M7_TME/combined_df.rds")

# ── 6. 可视化：免疫细胞比例热图 ──────────────────────────
cat("🎨 绘制免疫细胞比例热图...\n")

# 提取主要免疫细胞列（ssGSEA 结果列）
# IOBR ssGSEA 输出列名通常以 "_GSVA" 或 "_ssGSEA" 结尾
immune_cols <- grep("ssGSEA|GSVA|_TME|Immune|T_cell|CD8|CD4|NK|Macro|Neutro|Dendr|Mast|Treg|B_cell",
                    colnames(immune_ssgsea), value = TRUE, ignore.case = TRUE)

if (length(immune_cols) == 0) {
  # 如果没匹配到，取除 sample_id 外的所有数值列
  immune_cols <- colnames(immune_ssgsea)[sapply(immune_ssgsea, is.numeric)]
}

immune_mat <- as.matrix(immune_ssgsea[, immune_cols])
rownames(immune_mat) <- rownames(immune_ssgsea)

# 缩短列名（去除后缀）
colnames(immune_mat) <- gsub("_GSVA_|_ssGSEA_|_TME_", " ", colnames(immune_mat))
colnames(immune_mat) <- gsub("\\.", " ", colnames(immune_mat))

# 对矩阵做 z-score 标准化（行方向，即每个细胞类型在样本间标准化）
# 原因：不同细胞类型的基线分数差异大，标准化后颜色才能反映相对高低
immune_mat_scaled <- t(scale(t(immune_mat)))

# 只展示前 80 个样本（避免热图太密）
n_show <- min(80, nrow(immune_mat_scaled))
mat_show <- immune_mat_scaled[sample(nrow(immune_mat_scaled), n_show), ]

png("D:/Bio-Informatics Case Study/M7_TME/figures/01_immune_heatmap.png",
    width = 14, height = 10, units = "in", res = 150)

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#D6604D"))

ht <- Heatmap(
  t(mat_show),
  name = "z-score",
  col = col_fun,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  column_title = sprintf("TCGA-LAML 免疫细胞浸润热图（ssGSEA，n=%d）", n_show),
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(title = "z-score", title_gp = gpar(fontsize = 10)),
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson"
)

draw(ht)
dev.off()
cat("✅ 热图保存：figures/01_immune_heatmap.png\n")

# ── 7. 可视化：靶点 vs CD8+ T 细胞相关性 ────────────────
cat("🎨 绘制靶点 vs 免疫细胞相关性散点图...\n")

# 找 CD8+ T 细胞的列名
cd8_col <- grep("CD8|Cytotoxic|cytotoxic", colnames(combined_df), value = TRUE)[1]
cat(sprintf("  CD8+ T 细胞列：%s\n", cd8_col))

if (!is.na(cd8_col)) {
  plots <- lapply(car_t_targets, function(gene) {
    if (!gene %in% colnames(combined_df)) return(NULL)

    ggscatter(
      combined_df,
      x = gene,
      y = cd8_col,
      add = "reg.line",           # 加回归线
      conf.int = TRUE,            # 95% 置信区间
      cor.coef = TRUE,            # 显示 Pearson r
      cor.method = "pearson",
      color = "#2166AC",
      alpha = 0.5,
      size = 1.5,
      xlab = paste0(gene, " 表达量（log2 TPM+1）"),
      ylab = "CD8+ T 细胞浸润评分（ssGSEA）",
      title = paste0(gene, " vs CD8+ T 细胞"),
      subtitle = "正相关 = 高表达肿瘤有更多 CD8+ T 细胞浸润（有利于 CAR-T）"
    ) +
      theme(plot.title = element_text(face = "bold", size = 12))
  })

  plots <- plots[!sapply(plots, is.null)]

  combined_plot <- ggarrange(
    plotlist = plots,
    ncol = 3, nrow = 2,
    labels = LETTERS[1:length(plots)]
  )

  ggsave(
    "D:/Bio-Informatics Case Study/M7_TME/figures/02_target_vs_CD8T.png",
    combined_plot,
    width = 15, height = 10, dpi = 150
  )
  cat("✅ 相关性散点图保存：figures/02_target_vs_CD8T.png\n")
}

# ── 8. 可视化：免疫细胞堆叠柱状图（Top 20 样本）──────────
cat("🎨 绘制免疫细胞组成堆叠图...\n")

# 取主要免疫细胞（列数较多时取最重要的几个）
main_cells <- immune_cols[1:min(8, length(immune_cols))]

stack_df <- immune_ssgsea[, main_cells, drop = FALSE]
# 把负值截断为 0（ssGSEA 分数可能有负值，但比例图中无意义）
stack_df[stack_df < 0] <- 0
# 归一化为比例
stack_df <- stack_df / rowSums(stack_df + 1e-10)
stack_df$sample_id <- rownames(stack_df)

# 取前 40 个样本展示
stack_df_show <- stack_df[1:min(40, nrow(stack_df)), ]

stack_long <- stack_df_show %>%
  pivot_longer(-sample_id, names_to = "cell_type", values_to = "proportion")
stack_long$cell_type <- gsub("_GSVA_|_ssGSEA_|_TME_", " ", stack_long$cell_type)

p_stack <- ggplot(stack_long, aes(x = sample_id, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "TCGA-LAML 样本免疫细胞组成（ssGSEA，前40样本）",
    x = "样本", y = "相对比例", fill = "免疫细胞类型"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 8)
  )

ggsave(
  "D:/Bio-Informatics Case Study/M7_TME/figures/03_immune_composition.png",
  p_stack,
  width = 12, height = 6, dpi = 150
)
cat("✅ 堆叠图保存：figures/03_immune_composition.png\n")

# ── 9. 可视化：靶点表达 vs 多种免疫细胞相关性热图 ────────
cat("🎨 绘制靶点-免疫细胞相关性矩阵热图...\n")

target_available <- car_t_targets[car_t_targets %in% colnames(combined_df)]
immune_available <- immune_cols[immune_cols %in% colnames(combined_df)]
immune_available <- immune_available[1:min(15, length(immune_available))]

# 计算 Pearson 相关矩阵
cor_mat <- matrix(NA, nrow = length(target_available), ncol = length(immune_available))
pval_mat <- matrix(NA, nrow = length(target_available), ncol = length(immune_available))
rownames(cor_mat) <- rownames(pval_mat) <- target_available
colnames(cor_mat) <- colnames(pval_mat) <- immune_available

for (g in target_available) {
  for (ic in immune_available) {
    ct <- cor.test(combined_df[[g]], combined_df[[ic]], method = "pearson")
    cor_mat[g, ic] <- ct$estimate
    pval_mat[g, ic] <- ct$p.value
  }
}

# 清理列名
colnames(cor_mat) <- gsub("_GSVA_|_ssGSEA_|_TME_", " ", colnames(cor_mat))
colnames(pval_mat) <- colnames(cor_mat)

# 显著性标记
sig_mat <- ifelse(pval_mat < 0.001, "***",
           ifelse(pval_mat < 0.01,  "**",
           ifelse(pval_mat < 0.05,  "*", "")))

png("D:/Bio-Informatics Case Study/M7_TME/figures/04_target_immune_correlation.png",
    width = 12, height = 5, units = "in", res = 150)

col_fun2 <- colorRamp2(c(-0.5, 0, 0.5), c("#2166AC", "white", "#D6604D"))

ht2 <- Heatmap(
  cor_mat,
  name = "Pearson r",
  col = col_fun2,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 12, col = "black"))
  },
  row_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  column_title = "CAR-T 靶点 × 免疫细胞浸润相关性（*p<0.05, **p<0.01, ***p<0.001）",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(title = "Pearson r")
)

draw(ht2)
dev.off()
cat("✅ 相关性矩阵热图保存：figures/04_target_immune_correlation.png\n")

cat("\n🎉 M7 分析完成！所有图表已保存到 D:/Bio-Informatics Case Study/M7_TME/figures/\n")
cat("图表清单：\n")
cat("  01_immune_heatmap.png       - 样本免疫细胞浸润热图\n")
cat("  02_target_vs_CD8T.png       - 靶点 vs CD8+ T 细胞散点图\n")
cat("  03_immune_composition.png   - 免疫细胞组成堆叠图\n")
cat("  04_target_immune_correlation.png - 靶点-免疫细胞相关性矩阵\n")
