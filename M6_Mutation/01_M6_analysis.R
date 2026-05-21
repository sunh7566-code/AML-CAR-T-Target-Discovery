# ============================================================
# M6: 基因组变异分析（WES/MAF）
# 项目：AML CAR-T 靶点发现 pipeline
# 数据：TCGA-LAML somatic mutation (131 samples, 3900 mutations)
# 作者：Hao
# 日期：2026-05-21
# ============================================================

# ── 0. 加载包 ────────────────────────────────────────────────
library(maftools)
library(TCGAbiolinks)
library(httr)
library(ggplot2)

# 创建输出目录
dir.create("D:/Bio-Informatics Case Study/M6_Mutation/data/GDCdata", recursive = TRUE, showWarnings = FALSE)
dir.create("D:/Bio-Informatics Case Study/M6_Mutation/figures", showWarnings = FALSE)

# ── 1. 下载 TCGA-LAML 体细胞突变数据 ─────────────────────────
# 说明：GDC 按病人分开存储 MAF 文件，需要用 TCGAbiolinks 批量下载并合并
# 不能直接用单个 file_id 下载，那样只会得到单样本文件

query_mut <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)
# 共 153 个文件，对应 131 个唯一病人（部分病人有多个测序批次）

GDCdownload(
  query = query_mut,
  directory = "D:/Bio-Informatics Case Study/M6_Mutation/data/GDCdata"
)

# 合并所有样本为一个数据框
maf_all <- GDCprepare(
  query = query_mut,
  directory = "D:/Bio-Informatics Case Study/M6_Mutation/data/GDCdata"
)
# 结果：3900 条突变记录，140 个字段，131 个唯一样本

# ── 2. 读取为 maftools MAF 对象 ───────────────────────────────
# MAF 对象是 maftools 的核心数据结构，后续所有分析都基于它
laml_full <- read.maf(maf = maf_all)

# 基本统计
cat("=== TCGA-LAML 突变数据摘要 ===\n")
cat("样本数：", nrow(getSampleSummary(laml_full)), "\n")
cat("突变基因数：", nrow(getGeneSummary(laml_full)), "\n")
cat("\nTop 10 最常突变基因：\n")
print(head(getGeneSummary(laml_full)[, c("Hugo_Symbol", "MutatedSamples")], 10))

# ── 3. 瀑布图 / Oncoprint（Top 20 突变基因）─────────────────
# 每列=一个病人，每行=一个基因，颜色=突变类型
# 右侧条形图=该基因在多少%病人中发生突变
png("D:/Bio-Informatics Case Study/M6_Mutation/figures/01_waterfall_top20.png",
    width = 3000, height = 2000, res = 200)
oncoplot(maf = laml_full, top = 20)
dev.off()

pdf("D:/Bio-Informatics Case Study/M6_Mutation/figures/01_waterfall_top20.pdf",
    width = 14, height = 10)
oncoplot(maf = laml_full, top = 20)
dev.off()

cat("✅ 图1：瀑布图已保存\n")

# ── 4. Lollipop 图（突变位点在蛋白质上的分布）──────────────
# 横轴=蛋白质氨基酸序列，棒棒糖高度=突变频率，颜色=突变类型
# 可以看出突变是否集中在特定功能域（hotspot）

# FLT3：突变集中在 PKc_like 激酶域（位置 600-800），全为错义突变
# 临床意义：激酶域突变导致持续激活，是 midostaurin 的靶点
png("D:/Bio-Informatics Case Study/M6_Mutation/figures/02_lollipop_FLT3.png",
    width = 2400, height = 1200, res = 200)
lollipopPlot(maf = laml_full, gene = "FLT3")
dev.off()

# NPM1：10 个突变高度聚集在蛋白质末端（第 294 位），移码插入为主
# 临床意义：破坏核定位信号，蛋白质从细胞核"逃逸"到细胞质，AML WHO 独立亚型
png("D:/Bio-Informatics Case Study/M6_Mutation/figures/03_lollipop_NPM1.png",
    width = 2400, height = 1200, res = 200)
lollipopPlot(maf = laml_full, gene = "NPM1")
dev.off()

# DNMT3A：突变分散在多个功能域，包括错义、剪接位点、移码缺失
# 临床意义：DNA 甲基化酶突变导致表观遗传失调（与 M8 甲基化分析联动）
png("D:/Bio-Informatics Case Study/M6_Mutation/figures/04_lollipop_DNMT3A.png",
    width = 2400, height = 1200, res = 200)
lollipopPlot(maf = laml_full, gene = "DNMT3A")
dev.off()

cat("✅ 图2-4：Lollipop 图已保存\n")

# ── 5. 突变共现 / 互斥矩阵 ───────────────────────────────────
# 用 Fisher 精确检验判断两基因突变是否显著共现（深绿）或互斥（金黄）
# 关键发现：FLT3 + DNMT3A 显著共现（p<0.05），双突变病人预后更差
png("D:/Bio-Informatics Case Study/M6_Mutation/figures/05_somatic_interactions.png",
    width = 2400, height = 2000, res = 200)
somaticInteractions(maf = laml_full, top = 20, pvalue = c(0.05, 0.1))
dev.off()

cat("✅ 图5：共现矩阵已保存\n")

# ── 6. FLT3 突变亚组 vs mRNA 表达量箱线图 ───────────────────
# 整合 M5 生存分析的表达数据，验证 DNA 突变是否对应 RNA 表达变化

# 提取 FLT3 突变样本 ID
flt3_mutated <- subsetMaf(maf = laml_full, genes = "FLT3")
flt3_mut_samples <- getSampleSummary(flt3_mutated)$Tumor_Sample_Barcode
flt3_mut_ids <- substr(as.character(flt3_mut_samples), 1, 12)  # 取前12位病人ID

# 加载 M5 表达 + 临床整合数据框
analysis_df <- readRDS("D:/Bio-Informatics Case Study/M5_Survival/analysis_df.rds")

# 添加 FLT3 突变状态标签
analysis_df$FLT3_mutation <- ifelse(
  analysis_df$patient_id %in% flt3_mut_ids,
  "Mutated", "Wild-type"
)

cat("FLT3 突变组：", sum(analysis_df$FLT3_mutation == "Mutated"), "人\n")
cat("FLT3 野生型：", sum(analysis_df$FLT3_mutation == "Wild-type"), "人\n")

# Wilcoxon 检验
wt <- wilcox.test(FLT3 ~ FLT3_mutation, data = analysis_df)
cat("Wilcoxon p 值：", round(wt$p.value, 3), "\n")

# 画箱线图
png("D:/Bio-Informatics Case Study/M6_Mutation/figures/06_FLT3_mutation_vs_expression.png",
    width = 1600, height = 1600, res = 200)

ggplot(analysis_df, aes(x = FLT3_mutation, y = log2(FLT3 + 1),
                         fill = FLT3_mutation)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("Mutated" = "#E74C3C", "Wild-type" = "#3498DB")) +
  labs(
    title = "FLT3 mRNA Expression by Mutation Status",
    subtitle = "TCGA-LAML (n=130)",
    x = "FLT3 Mutation Status",
    y = "FLT3 Expression (log2 counts + 1)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  annotate("text", x = 1.5, y = max(log2(analysis_df$FLT3 + 1)) * 0.98,
           label = paste0("p = ", round(wt$p.value, 3)), size = 5)

dev.off()
cat("✅ 图6：FLT3 表达量箱线图已保存\n")

cat("\n=== M6 分析全部完成！===\n")
cat("输出目录：D:/Bio-Informatics Case Study/M6_Mutation/figures/\n")
cat("共生成 6 张图：\n")
cat("  01_waterfall_top20.png/pdf  - 瀑布图\n")
cat("  02_lollipop_FLT3.png        - FLT3 突变位点\n")
cat("  03_lollipop_NPM1.png        - NPM1 突变位点\n")
cat("  04_lollipop_DNMT3A.png      - DNMT3A 突变位点\n")
cat("  05_somatic_interactions.png - 突变共现矩阵\n")
cat("  06_FLT3_mutation_vs_expression.png - 表达量箱线图\n")
