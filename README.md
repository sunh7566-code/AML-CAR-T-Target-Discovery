# AML CAR-T Target Discovery — Bulk RNA-seq Differential Expression Analysis

## 项目简介 / Overview

本项目通过分析 GEO 公开数据集 **GSE6891**（537 个样本，包含 AML 患者与正常骨髓），利用 limma 差异表达分析框架，系统鉴定急性髓系白血病（AML）中的潜在 CAR-T 细胞治疗靶点。

This project analyzes the public GEO dataset **GSE6891** (537 samples: AML patients vs. normal bone marrow) using the limma differential expression framework to systematically identify candidate CAR-T therapy targets in Acute Myeloid Leukemia (AML).

---

## 数据来源 / Data Source

| 项目 | 内容 |
|------|------|
| **数据集** | GSE6891 |
| **平台** | Affymetrix Human Genome U133 Plus 2.0 Array (GPL570) |
| **样本数** | 537（AML: 461，Normal: 76） |
| **来源** | NCBI Gene Expression Omnibus (GEO) |
| **数据类型** | 微阵列（Microarray），RMA 归一化后的表达值 |

---

## 分析方法 / Methods

1. **数据下载**：GEOquery 从 GEO 下载 GSE6891 series matrix 文件
2. **样本分组**：根据 `characteristics_ch1` 字段自动识别 AML vs Normal
3. **差异分析**：limma 线性模型 + eBayes 经验贝叶斯收缩
4. **多重检验校正**：Benjamini-Hochberg (BH) 方法，控制 FDR
5. **探针注释**：hgu133plus2.db 将 Affymetrix 探针 ID 转换为基因 symbol
6. **去重**：每个基因保留最显著探针（`slice_min(adj.P.Val)`）
7. **可视化**：EnhancedVolcano 生成专业火山图

**筛选阈值**：|logFC| > 1 且 adj.P.Val < 0.05

---

## 主要结果 / Key Results

- 共检测到显著差异表达基因 **17,000+** 个
- **已知 AML CAR-T 靶点验证**（正向对照）：

| 基因 | logFC | 生物学意义 |
|------|-------|-----------|
| CLEC12A (CLL-1) | +2.78 | AML 高表达，临床 CAR-T 开发中 |
| HAVCR2 (TIM-3) | +2.15 | AML 免疫检查点，CAR-T 靶点候选 |
| IL3RA (CD123) | 上调 | AML/BPDCN 靶点，临床 II 期 |
| FLT3 | 上调 | AML 突变型靶点，临床 CAR-T 中 |

---

## 文件结构 / Repository Structure

```
AML-CAR-T-Target-Discovery/
├── AML Differential Expression Analysis.Rmd   # 主分析脚本（RMarkdown）
├── AML Differential Expression Analysis.html  # 生成的 HTML 报告
├── Bioinformatics_Pipeline_SOP.md             # 标准操作流程（SOP）
├── CLAUDE.md                                  # 项目配置文档
├── M1_bulk_RNAseq/
│   ├── data/          # 原始数据（.gitignore 排除，需自行从 GEO 下载）
│   └── plots/         # 输出图表
└── README.md
```

---

## 环境依赖 / Dependencies

```r
# R 4.6.0
library(GEOquery)          # 下载 GEO 数据
library(limma)             # 差异分析
library(hgu133plus2.db)    # 探针注释
library(EnhancedVolcano)   # 火山图
library(tidyverse)         # 数据处理
library(kableExtra)        # 表格美化
```

---

## 如何复现 / How to Reproduce

1. 从 GEO 手动下载 `GSE6891_series_matrix.txt.gz`，放入 `M1_bulk_RNAseq/data/`
2. 用 RStudio 打开 `AML Differential Expression Analysis.Rmd`
3. 点击 **Knit** 按钮，生成 HTML 报告

---

## 里程碑进度 / Milestone Progress

- [x] **M1**：Bulk RNA-seq 差异分析（✅ 完成）
- [ ] **M2**：单细胞 RNA-seq 靶点定位（Scanpy）
- [ ] **M3**：蛋白质结构预测 + 分子对接（AlphaFold + HDOCK）
- [ ] **M4**：临床数据库 + 因果推断（选做）

---

## 作者 / Author

**HAO** | Pharmacy/Medical Background | Bioinformatics Learner  
GitHub: [@sunh7566-code](https://github.com/sunh7566-code)
