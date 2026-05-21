# ============================================================
# M7 第二阶段：免疫细胞反卷积 + 可视化
# 直接加载缓存的 TPM 矩阵，跳过 barcode 转换
# ============================================================

library(IOBR)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

cat("✅ 包加载完成\n")

# ── 1. 加载缓存的 TPM 矩阵 ────────────────────────────────
# 原因：M5 数据已经读好缓存，不需要重新读 151 个文件
cat("📂 加载 TPM 矩阵缓存...\n")
tpm_mat <- readRDS("D:/Bio-Informatics Case Study/M7_TME/tpm_matrix.rds")
cat(sprintf("✅ TPM 矩阵：%d 基因 × %d 样本\n", nrow(tpm_mat), ncol(tpm_mat)))

# ── 2. 准备 log2(TPM+1) 矩阵 ─────────────────────────────
# 原因：IOBR 要求 log2 转换后的输入
#       log2(TPM+1) 是生信标准做法：+1 是为了避免 log(0) = -Inf
tpm_log2 <- log2(tpm_mat + 1)

# 确认靶点基因都在矩阵里
car_t_targets <- c("FLT3", "CD33", "IL3RA", "CLEC12A", "CD38")
found <- car_t_targets[car_t_targets %in% rownames(tpm_log2)]
cat(sprintf("✅ 找到 CAR-T 靶点：%s\n", paste(found, collapse = ", ")))

# ── 3. 运行 ssGSEA 免疫细胞反卷积 ────────────────────────
# 原理：ssGSEA 用已知免疫细胞的"特征基因集"给每个样本打分
#       分数越高 = 该细胞类型在肿瘤中浸润越多
#       药学类比：就像用一套抗体 panel 做流式，从基因表达推算细胞组成
cat("\n🔬 运行 ssGSEA 反卷积（约需 3-8 分钟，请耐心等待）...\n")

immune_ssgsea <- deconvo_tme(
  eset    = tpm_log2,
  method  = "ssgsea",
  arrays  = FALSE,      # RNA-seq 数据
  scale_eset = TRUE     # 跨样本标准化
)

cat(sprintf("✅ ssGSEA 完成：%d 样本 × %d 指标\n",
            nrow(immune_ssgsea), ncol(immune_ssgsea)))
print(head(colnames(immune_ssgsea), 20))  # 看看列名长什么样

saveRDS(immune_ssgsea, "D:/Bio-Informatics Case Study/M7_TME/immune_ssgsea.rds")

# ── 4. 整合靶点表达 + 免疫细胞分数 ──────────────────────
cat("\n🔗 整合靶点表达与免疫细胞分数...\n")

# 提取靶点表达（转置：行=样本，列=基因）
target_df <- as.data.frame(t(tpm_log2[found, , drop = FALSE]))
target_df$ID <- rownames(target_df)

# ssGSEA 结果的 ID 列
immune_df <- immune_ssgsea
# IOBR 2.x 输出里行名就是样本 ID
immune_df$ID <- rownames(immune_df)

combined_df <- inner_join(target_df, immune_df, by = "ID")
cat(sprintf("✅ 整合后：%d 样本可用于相关性分析\n", nrow(combined_df)))

saveRDS(combined_df, "D:/Bio-Informatics Case Study/M7_TME/combined_df.rds")

# ── 5. 确认免疫细胞列名 ──────────────────────────────────
# IOBR 2.x ssGSEA 输出的列名格式：Cell_Type_ssGSEA
all_cols <- colnames(immune_df)
immune_cols <- all_cols[all_cols != "ID"]

cat(sprintf("\n📋 共 %d 个免疫细胞指标：\n", length(immune_cols)))
print(immune_cols)

# 找 CD8+ T 细胞相关列
cd8_col <- grep("CD8|Cytotoxic|cytotoxic|T_cell_CD8|TCD8",
                immune_cols, value = TRUE, ignore.case = TRUE)
cat(sprintf("\n🔍 CD8+ T 细胞相关列：%s\n", paste(cd8_col, collapse = ", ")))

# ── 6. 图1：免疫细胞热图 ─────────────────────────────────
cat("\n🎨 绘制图1：免疫细胞浸润热图...\n")

immune_mat <- as.matrix(immune_df[, immune_cols])
rownames(immune_mat) <- immune_df$ID

# z-score 标准化（按列，即每种细胞类型在样本间标准化）
immune_mat_z <- scale(immune_mat)
# 截断极端值
immune_mat_z[immune_mat_z > 3]  <-  3
immune_mat_z[immune_mat_z < -3] <- -3

# 随机取 60 个样本展示（全部 151 个会太密）
set.seed(42)
n_show <- min(60, nrow(immune_mat_z))
idx <- sample(nrow(immune_mat_z), n_show)
mat_show <- t(immune_mat_z[idx, ])

# 清理列名
colnames(mat_show) <- paste0("S", seq_len(ncol(mat_show)))
rownames(mat_show) <- gsub("_ssGSEA$|_TME$|_GSVA$", "", rownames(mat_show))
rownames(mat_show) <- gsub("_", " ", rownames(mat_show))

png("D:/Bio-Informatics Case Study/M7_TME/figures/01_immune_heatmap.png",
    width = 14, height = 10, units = "in", res = 150)

col_fun <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#D6604D"))

ht <- Heatmap(
  mat_show,
  name            = "z-score",
  col             = col_fun,
  show_column_names = FALSE,
  show_row_names  = TRUE,
  row_names_gp    = gpar(fontsize = 8),
  column_title    = sprintf("TCGA-LAML 免疫细胞浸润热图（ssGSEA，n=%d）", n_show),
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  clustering_distance_rows    = "pearson",
  clustering_distance_columns = "pearson",
  heatmap_legend_param = list(title = "z-score")
)
draw(ht)
dev.off()
cat("✅ 图1 保存：figures/01_immune_heatmap.png\n")

# ── 7. 图2：靶点 vs CD8+ T 细胞散点图 ───────────────────
cat("\n🎨 绘制图2：靶点 vs CD8+ T 细胞散点图...\n")

if (length(cd8_col) > 0) {
  use_cd8 <- cd8_col[1]

  plots <- lapply(found, function(gene) {
    ggscatter(
      combined_df,
      x = gene, y = use_cd8,
      add = "reg.line",
      conf.int = TRUE,
      cor.coef = TRUE,
      cor.method = "pearson",
      color = "#2166AC", alpha = 0.5, size = 1.5,
      xlab = paste0(gene, "\n(log2 TPM+1)"),
      ylab = "CD8+ T 细胞评分（ssGSEA）",
      title = gene
    ) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
  })

  p2 <- ggarrange(plotlist = plots, ncol = 3, nrow = 2,
                  labels = LETTERS[1:length(plots)])
  p2 <- annotate_figure(p2,
    top = text_grob("CAR-T 靶点表达 vs CD8+ T 细胞浸润（TCGA-LAML）",
                    face = "bold", size = 14))

  ggsave("D:/Bio-Informatics Case Study/M7_TME/figures/02_target_vs_CD8T.png",
         p2, width = 15, height = 10, dpi = 150)
  cat("✅ 图2 保存：figures/02_target_vs_CD8T.png\n")
} else {
  cat("⚠️ 未找到 CD8+ T 细胞列，跳过图2\n")
  # 如果真的找不到，用第一个免疫细胞列代替
  use_cd8 <- immune_cols[1]
}

# ── 8. 图3：免疫细胞组成堆叠图 ───────────────────────────
cat("\n🎨 绘制图3：免疫细胞组成堆叠图...\n")

# 取前 8 个免疫细胞类型（分数最高的）
top_cells <- immune_cols[1:min(8, length(immune_cols))]

stack_df <- as.data.frame(immune_mat[, top_cells, drop = FALSE])
stack_df[stack_df < 0] <- 0  # 负分截为 0
row_sums <- rowSums(stack_df)
stack_df <- stack_df / ifelse(row_sums == 0, 1, row_sums)  # 归一化为比例
stack_df$sample_id <- rownames(stack_df)

# 取前 50 个样本
stack_show <- stack_df[1:min(50, nrow(stack_df)), ]
stack_long  <- pivot_longer(stack_show, -sample_id,
                             names_to = "cell_type", values_to = "proportion")
stack_long$cell_type <- gsub("_ssGSEA$|_TME$", "", stack_long$cell_type)
stack_long$cell_type <- gsub("_", " ", stack_long$cell_type)

p3 <- ggplot(stack_long, aes(x = sample_id, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "TCGA-LAML 免疫细胞组成（ssGSEA，前50样本）",
    subtitle = "每列为一个样本，颜色表示不同免疫细胞类型的相对比例",
    x = "样本", y = "相对比例", fill = "免疫细胞类型"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title   = element_text(face = "bold", size = 12)
  )

ggsave("D:/Bio-Informatics Case Study/M7_TME/figures/03_immune_composition.png",
       p3, width = 12, height = 6, dpi = 150)
cat("✅ 图3 保存：figures/03_immune_composition.png\n")

# ── 9. 图4：靶点 × 免疫细胞相关性矩阵热图 ───────────────
cat("\n🎨 绘制图4：靶点-免疫细胞相关性矩阵...\n")

# 最多用 15 个免疫细胞类型
use_immune <- immune_cols[1:min(15, length(immune_cols))]

cor_mat  <- matrix(NA, nrow = length(found), ncol = length(use_immune),
                   dimnames = list(found, use_immune))
pval_mat <- cor_mat

for (g in found) {
  for (ic in use_immune) {
    if (g %in% colnames(combined_df) && ic %in% colnames(combined_df)) {
      ct <- cor.test(combined_df[[g]], combined_df[[ic]], method = "pearson")
      cor_mat[g, ic]  <- ct$estimate
      pval_mat[g, ic] <- ct$p.value
    }
  }
}

# 去掉全 NA 的列
valid_cols <- colSums(!is.na(cor_mat)) > 0
cor_mat  <- cor_mat[,  valid_cols, drop = FALSE]
pval_mat <- pval_mat[, valid_cols, drop = FALSE]

# 清理列名
colnames(cor_mat) <- gsub("_ssGSEA$|_TME$|_GSVA$", "", colnames(cor_mat))
colnames(cor_mat) <- gsub("_", " ", colnames(cor_mat))
colnames(pval_mat) <- colnames(cor_mat)

# 显著性标记
sig_mat <- ifelse(pval_mat < 0.001, "***",
           ifelse(pval_mat < 0.01,  "**",
           ifelse(pval_mat < 0.05,  "*", "")))

png("D:/Bio-Informatics Case Study/M7_TME/figures/04_target_immune_correlation.png",
    width = 12, height = 4, units = "in", res = 150)

col_fun2 <- colorRamp2(c(-0.5, 0, 0.5), c("#2166AC", "white", "#D6604D"))

ht2 <- Heatmap(
  cor_mat,
  name  = "Pearson r",
  col   = col_fun2,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 14, fontface = "bold"))
  },
  row_names_gp    = gpar(fontsize = 12, fontface = "bold"),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  column_title    = "CAR-T 靶点 × 免疫细胞浸润相关性（*p<0.05  **p<0.01  ***p<0.001）",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(title = "Pearson r")
)
draw(ht2)
dev.off()
cat("✅ 图4 保存：figures/04_target_immune_correlation.png\n")

# ── 10. 汇总结果 ─────────────────────────────────────────
cat("\n", strrep("=", 50), "\n")
cat("🎉 M7 分析全部完成！\n")
cat(strrep("=", 50), "\n\n")

cat("📊 各靶点与 CD8+ T 细胞的相关性：\n")
if (use_cd8 %in% colnames(cor_mat)) {
  cd8_clean <- gsub("_ssGSEA$|_TME$|_GSVA$", "", use_cd8)
  cd8_clean <- gsub("_", " ", cd8_clean)
  for (g in found) {
    r   <- round(cor_mat[g, cd8_clean], 3)
    sig <- sig_mat[g, cd8_clean]
    cat(sprintf("  %-10s: r = %6.3f %s\n", g, r, sig))
  }
}

cat("\n📁 输出文件：\n")
cat("  M7_TME/figures/01_immune_heatmap.png\n")
cat("  M7_TME/figures/02_target_vs_CD8T.png\n")
cat("  M7_TME/figures/03_immune_composition.png\n")
cat("  M7_TME/figures/04_target_immune_correlation.png\n")
