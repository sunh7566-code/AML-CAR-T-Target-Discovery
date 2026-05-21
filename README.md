# AML CAR-T Target Discovery — Full Bioinformatics Portfolio

> **主题**：针对急性髓系白血病（AML）的新型 CAR-T 靶点发现、结构设计、临床验证与因果推断  
> **作者**：Hao | 药学/医学背景 | 系统学习生信 → pharma/biotech 求职  
> **GitHub**：[@sunh7566-code](https://github.com/sunh7566-code)

---

## 项目愿景

本项目从零构建一条完整的**转化生物信息学研究链**：

```
公共多组学数据（GEO / TCGA）
        ↓
M1  Bulk RNA-seq 差异分析          → 哪些基因在 AML 中异常表达？
        ↓
M2  单细胞 RNA-seq 靶点定位        → 哪些细胞类型在表达这些基因？
        ↓
M3  蛋白质结构预测 + 分子对接      → 靶点结构如何？抗体能结合吗？
        ↓
M4  临床数据库 API + 因果推断      → 这些靶点有多少临床试验？有因果证据吗？
        ↓
M5  生存分析 + 预后模型            → 高表达靶点的患者预后是否更差？
        ↓
M6  基因组变异分析（WES/MAF）      → AML 突变图谱如何影响靶点？
        ↓
M7  免疫微环境分析（TME）          → CAR-T 能在骨髓微环境里存活吗？
        ↓
M8  DNA 甲基化 / 表观基因组        → 靶点表达受表观遗传调控吗？
        ↓
M9  空间转录组                     → 靶点在骨髓组织中的空间分布？
        ↓
M10 多组学整合                     → 所有证据汇聚，最终靶点排名
        ↓
自动化 Pipeline（终极目标）        → 输入数据集 ID → 全流程自动输出报告
```

**每个里程碑独立可交付，完成 M1 即可写进简历。**

---

## 里程碑进度

| # | 主题 | 技术栈 | 状态 | 交付物 |
|---|------|--------|------|--------|
| **M1** | Bulk RNA-seq 差异分析 | R / limma | ✅ 完成 | 火山图 + HTML 报告 |
| **M2** | 单细胞 RNA-seq 靶点定位 | Python / Scanpy | ✅ 完成 | UMAP + 候选靶点清单 |
| **M3** | 蛋白质结构 + 分子对接 | AlphaFold3 / ClusPro | ✅ 完成 | 3D 结构图 + 结合能 |
| **M4** | 临床数据库 API + MR 因果推断 | Python / R TwoSampleMR | ✅ 完成 | 临床分析报告 + MR 图 |
| **M5** | 生存分析 + 预后模型 | R / survival / TCGAbiolinks | ✅ 完成 | K-M 曲线 + Cox 报告 |
| **M6** | 基因组变异分析（WES/MAF） | R / maftools | ✅ 完成 | 瀑布图 + oncoprint |
| **M7** | 免疫微环境分析（TME） | R / IOBR 2.2.1 / xCell | ✅ 完成 | 免疫细胞热图 + 相关性矩阵 |
| **M8** | DNA 甲基化 / 表观基因组 | R / minfi / ChAMP | 🔲 待开始 | 甲基化差异热图 |
| **M9** | 空间转录组 | Python / Squidpy | 🔲 待开始 | 空间表达图 |
| **M10** | 多组学整合 | Python / MOFA+ | 🔲 待开始 | 整合因子图 + 最终靶点排名 |

---

## 核心发现（已完成部分）

### M1 + M2 候选靶点

| 靶点 | M1 logFC | M2 AML/Normal 倍数 | 临床试验数（M4） | 综合评级 |
|------|----------|-------------------|----------------|---------|
| FLT3 | 上调 | 10.54x | 7 | ⭐⭐⭐ 高表达低开发，研究价值极高 |
| CD33 | 上调 | 9.92x | 30 | ⭐⭐⭐ 最成熟靶点，正向对照 |
| CD123 (IL3RA) | 上调 | 2.73x | 33 | ⭐⭐⭐ 临床最活跃 |
| CLEC12A | +2.78 | 2.10x | 4 | ⭐⭐⭐ M1 最显著，临床极少，新兴靶点 |
| CD38 | 上调 | 1.41x | 11 | ⭐⭐ 多癌种靶点 |

### M3 结构分析
- AlphaFold3 预测 CLEC12A 胞外域：pTM = 0.75（高可信）
- ClusPro 对接 Cluster 0：Members = 76，Lowest Energy = -749.9
- 结合界面：氢键 + 碳氢键 + 疏水接触 + 盐桥，四种相互作用

### M4 临床 + 因果
- AML CAR-T 临床试验共 90 条（ClinicalTrials.gov，2026-05-20）
- Phase I 占主导（47条），领域仍处早期探索阶段
- MR 分析（模拟演示）：FLT3 高表达 → AML 风险 OR ≈ 1.28

### M5 生存分析（TCGA-LAML，n=130）
- IL3RA：HR=1.711，p=0.00188 ⭐⭐（唯一显著独立预后因子）
- FLT3：HR=1.128，p=0.241（趋势但不显著）
- CD33 / CLEC12A / CD38：均不显著

### M6 基因组变异分析（TCGA-LAML，131样本，3900条突变）
- NPM1（8.4%）：移码插入，C端热点，AML WHO独立亚型
- DNMT3A（6.1%）：多类型突变，表观遗传失调
- FLT3（4.6%）：激酶域错义突变，midostaurin 靶点
- FLT3 + DNMT3A 显著共现（p<0.05），双突变预后差
- FLT3 突变不影响 mRNA 表达量（p=0.255），突变影响蛋白功能而非转录

### M7 免疫微环境分析（TCGA-LAML，xCell，n=151）
- AML 是典型"冷肿瘤"（Cold Tumor）：CD8+ T 细胞浸润极低，这是血液肿瘤的固有特征，不影响 CAR-T 策略
- **CD38** 与 CD8+ Tem（效应记忆 T 细胞）显著**负相关**（p<0.001）→ CD38 高表达 = 免疫抑制微环境，支持 Daratumumab + CAR-T 联合策略
- **CLEC12A / CD33** 与髓系细胞（中性粒细胞 \* / 单核细胞 \*\*）正相关 → 符合髓系抗原的生物学特性，验证靶点特异性
- **FLT3 / IL3RA**：与免疫细胞浸润无显著相关，靶点独立于微环境

---

## 目录结构

```
AML-CAR-T-Target-Discovery/
├── README.md                          # 项目主页（本文件）
├── CLAUDE.md                          # AI 协作上下文文档
├── Bioinformatics_Pipeline_SOP.md     # 全流程标准操作手册（技术细节）
│
├── M1_bulk_RNAseq/                    # ✅ Bulk RNA-seq 差异分析
│   ├── AML Differential Expression Analysis.Rmd
│   ├── AML-Differential-Expression-Analysis.html
│   ├── data/                          # 原始数据（.gitignore，需自行下载）
│   └── plots/
│
├── M2_scRNAseq/                       # ✅ 单细胞 RNA-seq
│   ├── M2_AML_scRNAseq_Analysis.ipynb
│   └── figures/
│
├── M3_Structure/                      # ✅ 蛋白结构 + 分子对接
│   ├── M3_Report.Rmd
│   ├── M3_Report.docx
│   └── figures/
│
├── M4_Clinical_MR/                    # ✅ 临床数据库 + 孟德尔随机化
│   ├── 01_clinicaltrials_api.py
│   ├── 02_analysis_viz.py
│   ├── M4_Part1_ClinicalTrials.ipynb
│   ├── M4_Part2_MendelianRandomization.Rmd
│   └── figures/
│
├── M5_Survival/                       # ✅ 生存分析
│   ├── 01_download_data.R
│   ├── 02_full_analysis.R
│   └── figures/                       # 5张K-M曲线 + Cox森林图
│
├── M6_Mutation/                       # ✅ 基因组变异分析
│   ├── 01_M6_analysis.R
│   ├── M6_Report.Rmd
│   └── figures/                       # 瀑布图 + lollipop图 + 共现矩阵 + 箱线图
├── M7_TME/                            # ✅ 免疫微环境
│   ├── 02_M7_deconvo_and_plots.R      # 正式分析脚本（xCell + 4张图）
│   ├── tpm_matrix.rds                 # TPM 矩阵缓存（复用 M5 数据）
│   ├── immune_xcell.rds               # xCell 反卷积结果缓存
│   └── figures/                       # 4张图（热图+散点图+堆叠图+相关性矩阵）
├── M8_Methylation/                    # 🔲 DNA 甲基化
├── M9_SpatialTranscriptomics/         # 🔲 空间转录组
└── M10_MultiOmics/                    # 🔲 多组学整合
```

---

## 数据来源

| 里程碑 | 数据集 | 来源 | 说明 |
|--------|--------|------|------|
| M1 | GSE6891 | GEO | AML vs 正常骨髓，537样本，微阵列 |
| M2 | GSE116256 | GEO | AML 单细胞图谱，1244个细胞 |
| M3 | PDB 8W9J / UniProt P20138 | PDB / UniProt | CLEC12A 晶体结构 |
| M4 | ClinicalTrials.gov / IEU Open GWAS | API | 临床试验 + GWAS 汇总 |
| M5 | TCGA-LAML | TCGAbiolinks | 200 例 AML，含生存数据 |
| M6 | TCGA-LAML MAF | TCGAbiolinks | AML 体细胞突变图谱 |
| M7 | TCGA-LAML（复用 M5 数据） | GDC API | xCell 免疫细胞反卷积（151样本 × 68指标） |
| M8 | GEO EPIC array | GEO | AML DNA 甲基化数据 |
| M9 | 10x Visium 骨髓 | GEO / 10x | 空间转录组（待确定数据集） |
| M10 | 上述所有 | 多源 | MOFA+ 多组学整合 |

---

## 环境依赖

```
Python:  conda env scrna（Scanpy 1.11.5 + Squidpy + MOFA+）
R:       4.6.0（limma / DESeq2 / survival / maftools / IOBR / minfi）
其他:    AlphaFold Server（网页）/ ClusPro（网页）
```

详见 `Bioinformatics_Pipeline_SOP.md`。

---

## 自动化 Pipeline 愿景

本项目的终极目标是将 M1–M10 整合为一个**半自动化分析系统**：

```
输入: 癌症类型 + GEO/TCGA 数据集 ID
  ↓
自动执行: 数据下载 → QC → 差异分析 → 单细胞 → 结构预测接口 → 临床查询 → 生存分析 → 变异 → TME → 甲基化 → 空间 → 多组学整合
  ↓
输出: 标准化 HTML 报告 + 图表包 + 候选靶点排名表
```

详细架构设计见 `Bioinformatics_Pipeline_SOP.md` Part 11。
