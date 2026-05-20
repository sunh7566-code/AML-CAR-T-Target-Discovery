# 🧬 生物信息学项目 — 会话初始化文档
# CLAUDE.md — 每次新对话必须第一时间读取此文件

> **⚠️ 强制规则**：每次开启新对话，Claude 必须先读取此文件，再回应任何请求。

---

## 👤 用户基本信息

| 项目 | 内容 |
|------|------|
| **姓名** | Hao |
| **母语** | 中文简体（优先用中文交流） |
| **英语** | 第二语言 |
| **编程背景** | Python 有基础，R 零基础 |
| **专业背景** | 药学 / 医学 |
| **学习目标** | 从零开始系统学习生信，建立 portfolio，面向 pharma/biotech 求职 |

---

## 💻 用户硬件环境（已确认）

```
CPU:    Intel Core Ultra 7 155H (16 核)
RAM:    32 GB
GPU:    Intel Arc 集成显卡（⚠️ 无 NVIDIA，不能本地跑 AlphaFold）
存储:    954 GB 总，约 500 GB 可用
系统:    Windows 11 家庭中文版
```

**关键结论：全程本地运行，不需要购买任何云服务。AlphaFold 用网页版 Server。**

---

## 📁 工作目录

- **项目根目录**：`D:\Bio-Informatics Case Study\`
- **SOP 文件**：`D:\Bio-Informatics Case Study\Bioinformatics_Pipeline_SOP.md`
- **本文件**：`D:\Bio-Informatics Case Study\CLAUDE.md`

---

## 🗺️ 项目总览

**主题**：针对急性髓系白血病（AML）的新型 CAR-T 靶点发现、结构设计与药效预测

**逻辑链**：
```
公共多组学数据 (GEO / TCGA)
        ↓
Bulk RNA-seq 差异分析 (DESeq2)
        ↓
单细胞 RNA-seq 精准定位 (Scanpy)
        ↓
候选靶点 shortlist
        ↓
AlphaFold 结构预测
        ↓
分子对接 / 抗原表位分析 (HDOCK / AutoDock Vina)
        ↓
临床转化与因果验证（选做）
```

---

## 🏁 四个里程碑

| 里程碑 | 主题 | 主要技术栈 | 交付物 | 预计时长 |
|--------|------|------------|--------|----------|
| **M1** | Bulk RNA-seq 差异分析 | R + DESeq2 | 火山图 + RMarkdown 报告 | 2–3 周 |
| **M2** | 单细胞 RNA-seq 靶点定位 | Python + Scanpy | UMAP 图 + 候选靶点清单 | 3–4 周 |
| **M3** | 蛋白质结构 + 分子对接 | AlphaFold + PyMOL + AutoDock Vina | 3D 结构图 + 结合能数据 | 2–3 周 |
| **M4** | 临床数据库 API + 因果推断（选做） | Python (requests) + R (TwoSampleMR) | 临床试验分析报告 + MR 因果图 | 3–4 周 |

**每个里程碑独立可交付，完成 M1 即可写进简历。**

---

## 🛠️ 已安装 / 待安装环境（按需安装策略）

| 环境 | 状态 | 用途 |
|------|------|------|
| Git for Windows（v2.54.0，装在 D:\Git，用 Git Bash 而非 PowerShell） | ✅ 已安装 | 版本管理 |
| VS Code（v1.120.0）+ Python/Jupyter/R/Markdown 扩展 | ✅ 已安装 | 编辑器 |
| Miniconda（装在 `D:\miniconda3`，勾选 PATH） | ✅ 已安装（conda 26.3.2） | Python 环境管理 |
| R 4.6.0（装在 `D:\R\R-4.6.0`，PATH 已写入 ~/.bashrc） | ✅ 已安装 | DESeq2 差异分析 |
| RStudio Desktop（免费开源版） | ✅ 已安装 | R 的 IDE |
| Rtools45（装在 `D:\rtools45`） | ✅ 已安装 | 编译 R 包（DESeq2 等需要） |
| R 包：DESeq2/GEOquery/EnhancedVolcano/tidyverse 等 M1 全套 | ✅ 已安装 | M1 差异分析 |
| conda env `scrna`（Scanpy 1.11.5 + leiden + CellTypist + Jupyter，装在 `D:\miniconda3\envs\scrna`） | ✅ 已安装（2026-05-18） | 单细胞分析 |
| conda env `strucbio`（PyMOL 开源版） | 待确认（M3 开始前） | 结构可视化 |
| AutoDock Vina（`C:\Tools\vina`，加入 PATH） | 待确认（M3 开始前） | 分子对接 |

> **更新提示**：每当某环境成功安装后，请告知 Claude，Claude 将更新此表格中的"状态"栏。

---

## 📊 推荐数据集

| 里程碑 | 数据集 | 说明 |
|--------|--------|------|
| M1 | GSE6891 | AML vs 正常骨髓，分组清晰，样本量大，适合教学 |
| M1 备选 | GSE37642 | AML 多亚型 RNA-seq，含正常对照 |
| M2 练手 | PBMC 3k（Scanpy 内置） | `sc.datasets.pbmc3k()` 一行下载 |
| M2 真实 | GSE116256 | AML 单细胞图谱，包含 AML 患者与正常骨髓细胞 |
| M3 靶点 | CD33（P20138）| AML 最成熟 CAR-T 靶点，UniProt 下载序列 |

---

## 🔑 DESeq2 核心筛选阈值（M1）

```
|log2FoldChange| > 1   （即 fold change > 2 倍）
padj < 0.05            （看 padj，不看 p-value）
低表达过滤：至少 3 个样本里有 ≥10 reads
```

---

## 🔑 CAR-T 靶点筛选标准（M2）

理想靶点必须满足：
1. ✅ 细胞表面表达（膜蛋白）
2. ✅ 肿瘤细胞高表达
3. ✅ 正常组织低/无表达（safety）
4. ✅ 抗原稳定（不容易丢失）

**正向对照靶点**（跑完看列表里有没有这些）：
- CD33（AML，临床试验活跃，Gemtuzumab 已批）
- FLT3（AML，突变型，临床 CAR-T 中）
- CD123（IL3RA，AML/BPDCN，临床 II 期）
- CLL-1（CLEC12A，AML，临床中）
- CD38（多发性骨髓瘤，Daratumumab 已批）

---

## 🧪 分子对接工具选择（M3）

| 对接类型 | 正确工具 |
|----------|----------|
| 蛋白-蛋白对接（抗体-抗原） | HDOCK（网页）或 ClusPro |
| 小分子-蛋白对接 | AutoDock Vina |
| 结构可视化 | PyMOL 开源版（conda 安装） |
| 高级可视化备选 | ChimeraX（免费学术版） |

> ⚠️ **CAR-T 主题应做蛋白-蛋白对接**，不是小分子对接！

---

## 🌐 网络优化（中国大陆）

```r
# R / Bioconductor 镜像（写入 ~/.Rprofile）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
```

```bash
# Conda 换源（如果 conda install 很慢）
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
conda config --set show_channel_urls yes
```

---

## ⚠️ 常见坑（避免重蹈覆辙）

1. **Miniconda 安装路径**：不要带中文和空格（正确：`D:\miniconda3`）
2. **DESeq2 需要 raw count**，不是 TPM/FPKM。GEO 下载后先看 series_matrix 说明
3. **PyMOL 有两个版本**：商业版（Schrödinger）和开源版，我们用 conda 装的**开源版**
4. **AlphaFold Server 配额**：每日约 30 个 job，够用，不需要本地部署
5. **分组变量识别**：每个 GEO 数据集 pheno 列名不同，需要人工看 `pData()` 输出
6. **Leiden 聚类**：`pip install scanpy[leiden]`，中括号不能省

---

## 📋 M1 完成检查清单

```
☑ R + RStudio 装好，所有包能正常加载（✅ 2026-05-17 完成）
☑ 从 GEO 成功下载 AML 数据集（GSE6891，537样本，✅ 2026-05-18 完成）
☑ 跑通 limma 差异分析，得到 AML vs 正常骨髓差异基因表（✅ 2026-05-18 完成）
☑ 画出专业火山图 PNG + PDF（✅ 2026-05-18 完成）
☑ 写完 RMarkdown 报告，Knit 成 HTML（✅ 2026-05-18 完成）
☑ 把所有代码 push 到 GitHub（✅ 2026-05-18 完成）
☑ 仓库 README 说明数据来源、方法、主要结果（✅ 2026-05-18 完成）
```

**M1 关键结果**：CLEC12A (logFC=2.78)、HAVCR2 (logFC=2.15) 为最显著上调靶点
**注意**：GSE6891 是微阵列数据（Affymetrix GPL570），用 limma 而非 DESeq2

---

## 📋 M2 完成检查清单

```
☑ Scanpy 环境装好（scrna env，D:\miniconda3\envs\scrna，✅ 2026-05-18 完成）
☑ 跑通 PBMC 3k 标准流程练手（✅ 2026-05-18 完成）
☑ 处理 GSE116256 真实 AML scRNA-seq 数据（1244个细胞，✅ 2026-05-18 完成）
☑ 出 UMAP 图（按聚类 + 按AML/Normal）（✅ 2026-05-18 完成）
☑ 输出 5 个 CAR-T 候选靶点（✅ 2026-05-18 完成）
☑ 整理成 Jupyter Notebook 推到 GitHub（✅ 2026-05-18 完成）
```

**M2 关键结果**（AML vs Normal 表达倍数）：
- FLT3：10.54x ⭐⭐⭐
- CD33：9.92x ⭐⭐⭐
- IL3RA：2.73x ⭐⭐⭐
- CLEC12A：2.10x ⭐⭐⭐
- CD38：1.41x ⭐⭐

**下一步：M4** — 临床数据库 API + 孟德尔随机化因果推断（选做）

---

## 📋 M3 完成检查清单

```
☑ AlphaFold3 预测 CLEC12A 胞外域结构（pTM=0.75，✅ 2026-05-19 完成）
☑ 下载 PDB 8W9J 实验晶体结构（CLEC12A + 50C1 Fab，3.50 Å）
☑ 提取 Chain C（CLEC12A）和 Chain L+H（50C1 Fab）
☑ ClusPro 蛋白-蛋白对接（Job 1436585，✅ 2026-05-19 完成）
☑ Discovery Studio 可视化 + 分子间相互作用分析（5张图）
☑ 撰写 M3_Report.Rmd + M3_Report.docx（✅ 2026-05-20 完成）
☑ 推送到 GitHub（✅ 2026-05-20 完成）
```

**M3 关键结果**：
- AlphaFold3 pTM = 0.75（高可信，与 8W9J 晶体结构高度吻合）
- ClusPro Cluster 0：Members = 76（最大聚类），Lowest Energy = -749.9
- 结合界面：CLEC12A CRD 顶面，含氢键、碳氢键、疏水接触、盐桥四种相互作用
- 工具链：AlphaFold Server → PDB 8W9J → ClusPro 2.1 → Discovery Studio 2021

---

## 📋 M4 完成检查清单

```
☑ 编写 ClinicalTrials.gov v2 API 抓取脚本（01_clinicaltrials_api.py，✅ 2026-05-20 完成）
☑ 编写数据清洗与可视化脚本（02_analysis_viz.py，5张图，✅ 2026-05-20 完成）
☑ 整理完整 Jupyter Notebook（M4_Part1_ClinicalTrials.ipynb，✅ 2026-05-20 完成）
☑ 编写孟德尔随机化 RMarkdown（M4_Part2_MendelianRandomization.Rmd，✅ 2026-05-20 完成）
○ 本地运行 01_clinicaltrials_api.py 生成真实数据（待执行）
○ Knit MR Rmd 生成 HTML 报告（待执行，需安装 TwoSampleMR）
○ 推送到 GitHub（待执行）
```

**M4 文件结构**（`D:\Bio-Informatics Case Study\M4_Clinical_MR\`）：
- `01_clinicaltrials_api.py` — API 抓取脚本（本地运行）
- `02_analysis_viz.py` — 5张可视化图
- `M4_Part1_ClinicalTrials.ipynb` — Jupyter Notebook
- `M4_Part2_MendelianRandomization.Rmd` — MR 因果分析
- `README_M4_运行指南.md` — 运行步骤说明

**M4 核心工具**：
- Part1：Python requests + pandas + matplotlib（ClinicalTrials.gov v2 API）
- Part2：R TwoSampleMR 包 + IEU Open GWAS（含离线模拟演示模式）

**下一步**：本地运行 `01_clinicaltrials_api.py` 获取真实数据 → Knit M4 Part2 Rmd → push GitHub

---

## 📚 核心学习资源

- DESeq2 官方 vignette：https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- Scanpy 教程：https://scanpy.readthedocs.io/en/stable/tutorials.html
- OSCA 单细胞分析书：https://bioconductor.org/books/release/OSCA/
- 中文社区：生信技能树（biotrainee.com）、生信菜鸟团（公众号）

---

## 🔄 Claude 工作规范（给自己的提示）

1. **语言**：默认用中文简体回复，代码注释可中英双语
2. **解释深度**：Hao 有药学/医学背景和 Python 基础，但 R 零基础。R 相关内容要解释得更细
3. **遇到报错**：引导 Hao 先看错误信息的关键部分，再给具体修复命令
4. **代码风格**：每段代码都要有中文注释，解释"为什么这样做"不只是"做了什么"
5. **里程碑进度**：主动跟踪当前在哪个里程碑，哪些步骤已完成
6. **文件保存**：所有输出文件保存到 `D:\Bio-Informatics Case Study\` 下对应子目录

---

*文件版本：v1.6 | 基于 Bioinformatics_Pipeline_SOP.md v1.1 生成*
*更新时间：2026-05-20（M1 ✅ + M2 ✅ + M3 ✅ + M4 脚本全部编写完成；待本地运行 API 脚本并 push GitHub）*
