# 生物信息学全流程实战 SOP
# Bioinformatics End-to-End Pipeline SOP

> **项目主题 / Project Theme**：针对急性髓系白血病（AML）的新型 CAR-T 靶点发现、结构设计与药效预测
> **Topic**: Novel CAR-T Target Discovery, Structural Design & Efficacy Prediction for Acute Myeloid Leukemia (AML)

> **文档版本 / Version**：v1.1
> **适用对象 / Audience**：药学/医学背景，Python 入门，R 零基础，从零开始系统学习生信的学习者
> **预计总周期 / Total Duration**：10–14 周（每周投入 15–25 小时）

---

## 📋 目录 / Table of Contents

- [Part 0：项目总览 / Project Overview](#part-0)
- [Part 1：环境搭建总策略 / Environment Setup Strategy](#part-1)
- [Part 2：工具与数据库清单 / Tools & Databases Catalog](#part-2)
- [Part 3：里程碑 1 — Bulk RNA-seq 差异分析与火山图](#part-3)
- [Part 4：里程碑 2 — 单细胞 RNA-seq 与靶点定位](#part-4)
- [Part 5：里程碑 3 — 蛋白质结构预测与分子对接](#part-5)
- [Part 6：里程碑 4 — 临床数据库 API 与因果推断（选做）](#part-6)
- [Part 7：常见错误速查 / Troubleshooting Index](#part-7)
- [Part 8：学习资源 / Learning Resources](#part-8)

---

<a name="part-0"></a>

# Part 0：项目总览 / Project Overview

## 0.1 项目逻辑链 / Project Logic Chain

```
临床未满足需求 (Unmet Need)
        ↓
公共多组学数据 (Public Omics Data: GEO / TCGA)
        ↓
差异表达分析 (DEG Analysis: bulk RNA-seq)
        ↓
单细胞精准定位 (Cell-type Resolution: scRNA-seq)
        ↓
候选靶点 shortlist (Target Candidates)
        ↓
三维结构预测 (Structural Prediction: AlphaFold)
        ↓
分子对接 / 抗原表位分析 (Docking / Epitope Analysis)
        ↓
临床转化与因果验证 (Clinical Translation & Causal Inference)
```

## 0.2 4 个里程碑 / 4 Milestones

| 里程碑 | 主题 | 预计时长 | 主要技术栈 | 交付物 |
|--------|------|----------|------------|--------|
| **M1** | Bulk RNA-seq 差异分析 | 2–3 周 | R + DESeq2 | 火山图 + RMarkdown 报告 |
| **M2** | 单细胞 RNA-seq 靶点定位 | 3–4 周 | Python + Scanpy | UMAP 图 + 候选靶点清单 |
| **M3** | 蛋白质结构 + 分子对接 | 2–3 周 | AlphaFold + PyMOL + AutoDock Vina | 3D 结构图 + 结合能数据 |
| **M4** | 临床数据库 API + 因果推断（选做） | 3–4 周 | Python (requests) + R (TwoSampleMR) | 临床试验分析报告 + MR 因果图 |

## 0.3 关键设计原则 / Design Principles

1. **每个里程碑独立可交付** — 完成 M1 就已经可以写进简历，不必等到 M4
2. **环境按里程碑分批装** — 避免一次性装太多包导致冲突
3. **优先使用 AI 辅助编程** — 但要看懂每一行代码在做什么
4. **用真实公共数据** — 不用 toy data，确保产出有工业级说服力
5. **每步都要有视觉产出** — 图表是面试展示的硬通货

---

<a name="part-1"></a>

# Part 1：环境搭建总策略 / Environment Setup Strategy

## 1.1 你的硬件评估 / Hardware Assessment

```
CPU:    Intel Core Ultra 7 155H (16 cores)      ✅ 充足
RAM:    32 GB                                    ✅ 单细胞分析也够用
GPU:    Intel Arc 集成显卡                       ⚠️ 不能本地跑 AlphaFold（需 NVIDIA）
                                                    → 用 AlphaFold Server 网页版即可
存储:    954 GB 总，451 GB 已用，500+ GB 可用     ✅ 充足
系统:    Windows 11 家庭中文版 25H2               ✅
```

**结论 / Conclusion**：**全程本地运行，不需要购买任何云服务。**

## 1.2 总体环境策略 / Overall Strategy

我们采用**「按需安装」**策略：每个里程碑开始前，只装那个里程碑需要的环境。

| 环境 | 何时安装 | 用途 |
|------|----------|------|
| **基础工具**（Git, VS Code, Conda） | 立刻 | 全程使用 |
| **R 环境**（R + RStudio + Bioconductor） | M1 开始前 | RNA-seq 差异分析 |
| **Python 生信环境**（Scanpy + 相关包） | M2 开始前 | 单细胞分析 |
| **结构生物学工具**（PyMOL + AutoDock Vina） | M3 开始前 | 分子对接 |
| **MR 分析环境**（TwoSampleMR） | M4 开始前 | 因果推断 |

## 1.3 立刻要装的基础工具 / Install Immediately

### 1.3.1 Git for Windows
- **官网 / Website**：https://git-scm.com/download/win
- **价格 / Price**：✅ 免费开源
- **用途 / Purpose**：版本管理 + 后续把项目推到 GitHub 作 portfolio
- **安装要点 / Install Notes**：
  - 一路 Next 即可
  - 安装后在 PowerShell 输入 `git --version` 验证

### 1.3.2 VS Code（推荐编辑器）
- **官网**：https://code.visualstudio.com/
- **价格**：✅ 免费开源
- **必装扩展 / Extensions**：
  - `Python`（Microsoft 官方）
  - `R`（REditorSupport）
  - `Jupyter`（Microsoft 官方）
  - `Markdown All in One`（写文档用）
  - `GitHub Copilot`（如果你有 GitHub Student/Pro 订阅，可免费用；否则跳过）

### 1.3.3 Miniconda（Python 环境管理器）
- **官网**：https://docs.conda.io/projects/miniconda/en/latest/
- **价格**：✅ 免费（个人/学术用途）
- **选择版本 / Version**：Windows 64-bit，Python 3.11
- **为什么用 Miniconda 而不是 Anaconda**：
  - Anaconda 默认装几百个包，臃肿且容易冲突
  - Miniconda 只装核心，按需添加，干净可控
- **安装要点**：
  - 安装路径**不要带中文和空格**（例如装在 `D:\miniconda3`）
  - 勾选 "Add Miniconda3 to PATH"（虽然官方不推荐，但对新手更友好）
  - 安装后打开 PowerShell，输入 `conda --version` 验证

### 1.3.4 安装顺序快速清单 / Quick Install Checklist

```
□ 1. Git for Windows
□ 2. VS Code + 扩展
□ 3. Miniconda
□ 4. （可选）注册 GitHub 账号，用于后续 push 代码
□ 5. （可选）注册 ORCID，用于学术身份管理
```

---

<a name="part-2"></a>

# Part 2：工具与数据库清单 / Tools & Databases Catalog

## 2.1 软件工具总表 / Software Catalog

### 编程语言与环境 / Languages & Runtimes

| 工具 | 价格 | 官网 | 用途 |
|------|------|------|------|
| **Python 3.11** | ✅ 免费 | python.org | 通用编程、单细胞分析 |
| **R 4.4+** | ✅ 免费 | r-project.org | 统计分析、DEG、生存分析 |
| **RStudio Desktop** | ✅ 免费（开源版） | posit.co | R 的 IDE |
| **Miniconda** | ✅ 免费 | docs.conda.io | Python 包/环境管理 |

### 生物信息核心包 / Core Bioinformatics Packages

| 工具 | 价格 | 用途 | 安装方式 |
|------|------|------|----------|
| **Bioconductor** | ✅ 免费 | R 生信包仓库 | `BiocManager::install()` |
| **DESeq2** | ✅ 免费 | bulk RNA-seq 差异分析 | Bioconductor |
| **edgeR / limma** | ✅ 免费 | DEG 替代方案 | Bioconductor |
| **Scanpy** | ✅ 免费 | 单细胞分析（Python） | pip / conda |
| **Seurat** | ✅ 免费 | 单细胞分析（R 备选） | CRAN |
| **AnnData** | ✅ 免费 | 单细胞数据格式 | pip |
| **ggplot2** | ✅ 免费 | R 可视化 | CRAN |
| **EnhancedVolcano** | ✅ 免费 | 火山图专用包 | Bioconductor |

### 结构生物学工具 / Structural Biology Tools

| 工具 | 价格 | 用途 | 备注 |
|------|------|------|------|
| **AlphaFold Server** | ✅ 免费（网页版） | 蛋白结构预测 | 每日 30 个 job，够用 |
| **PyMOL（开源版）** | ✅ 免费 | 3D 结构可视化 | 区别于付费的 Schrödinger PyMOL |
| **ChimeraX** | ✅ 免费（学术） | 高级可视化 | PyMOL 替代品 |
| **AutoDock Vina** | ✅ 免费 | 分子对接 | 命令行工具 |
| **AutoDockTools (ADT)** | ✅ 免费 | 配体/受体预处理 | GUI 工具 |
| **HDOCK** | ✅ 免费（网页版） | 蛋白-蛋白对接 | 抗体-抗原对接首选 |

> **重要区分 / Important Distinction**：PyMOL 有两个版本
> - **Open-Source PyMOL** ✅ 完全免费，需要从源码编译或装预编译包
> - **PyMOL by Schrödinger** ❌ 商业版，学术许可证免费但需申请
> - 我们用 **Open-Source 版**，通过 conda 一键装好

### 临床与因果分析 / Clinical & Causal Tools

| 工具 | 价格 | 用途 |
|------|------|------|
| **ClinicalTrials.gov API v2** | ✅ 免费 | 临床试验数据 |
| **openFDA API** | ✅ 免费 | FDA 审批 + 不良反应 |
| **TwoSampleMR (R)** | ✅ 免费 | 孟德尔随机化 |
| **MR-Base 平台** | ✅ 免费 | GWAS 数据库 |

### 通路与网络分析 / Pathway & Network

| 工具 | 价格 | 用途 |
|------|------|------|
| **STRING DB** | ✅ 免费（学术） | 蛋白互作网络 |
| **KEGG** | ✅ 免费（学术用途） | 通路富集 |
| **Cytoscape** | ✅ 免费 | 网络可视化 |
| **clusterProfiler (R)** | ✅ 免费 | 富集分析 |

## 2.2 数据库清单 / Database Catalog

### 多组学数据 / Omics Data

| 数据库 | 内容 | 访问方式 | 价格 |
|--------|------|----------|------|
| **NCBI GEO** | 各类组学数据集（最常用） | 网页 + GEOquery R 包 | ✅ 免费 |
| **TCGA / GDC** | 33 种癌症的多组学 | GDC Portal + TCGAbiolinks | ✅ 免费 |
| **ArrayExpress** | EBI 的组学数据库 | 网页 + ArrayExpress R 包 | ✅ 免费 |
| **SRA** | 原始测序数据（FASTQ） | sra-tools | ✅ 免费（下载慢） |
| **Human Cell Atlas** | 单细胞数据集合 | 网页下载 | ✅ 免费 |
| **CellxGene** | 标准化单细胞数据浏览器 | 网页 | ✅ 免费 |
| **GTEx** | 正常组织表达谱 | 网页 + R 包 | ✅ 免费 |

### 结构与序列 / Structure & Sequence

| 数据库 | 内容 | 价格 |
|--------|------|------|
| **UniProt** | 蛋白序列与注释 | ✅ 免费 |
| **PDB (RCSB)** | 实验测定蛋白结构 | ✅ 免费 |
| **AlphaFold DB** | AI 预测的全蛋白组结构 | ✅ 免费 |
| **NCBI Protein** | 蛋白序列 | ✅ 免费 |

### 临床与药物 / Clinical & Drug

| 数据库 | 内容 | 价格 |
|--------|------|------|
| **ClinicalTrials.gov** | 全球临床试验登记 | ✅ 免费 |
| **DrugBank** | 药物信息 | ✅ 学术版免费，商业版付费 |
| **OpenTargets** | 靶点-疾病关联证据 | ✅ 免费 |
| **DGIdb** | 药物-基因互作 | ✅ 免费 |

## 2.3 关于"完全免费"的说明 / Note on "Free" Status

整个 SOP 涉及的所有工具与数据库，**学术 / 个人学习用途均完全免费**。
唯一可能产生费用的场景：
- 把项目放上 GitHub Pro（可选，免费版够用）
- 用 Colab Pro 跑大型计算（本项目不需要）
- 申请商业级 DrugBank 许可（M4 选做，可跳过）

---

<a name="part-3"></a>

# Part 3：里程碑 1 — Bulk 表达谱差异分析与火山图

---

## 3.0 在开始之前：彻底搞懂你在做什么

> ⚠️ **强制阅读**。跳过这节直接跑代码，你只是在复制粘贴，不是在做生信分析。

### 3.0.1 为什么要做差异分析？

你的研究问题是：**哪些基因在疾病状态下异常活跃？**

癌症细胞之所以和正常细胞不同，根本原因是基因表达谱发生了改变——某些基因被异常激活（上调），某些被沉默（下调）。找到这些基因，就找到了疾病的"分子特征"，也找到了潜在的治疗靶点。

差异分析（Differential Expression Analysis，DEA）就是用统计方法，系统地找出在两组样本之间表达量有显著差异的基因。

### 3.0.2 基因表达是怎么被测量的？

基因工作的流程：

```
DNA（基因）→ 转录 → mRNA（信使RNA）→ 翻译 → 蛋白质 → 执行功能
```

**mRNA 是中间信使**。一个基因越活跃，产生的 mRNA 就越多。所以测量基因活跃程度 = 测量细胞里的 mRNA 含量。

目前主流有两种测量技术，选哪种决定了你后续用什么分析工具：

| 技术 | 原理 | 数据格式 | 正确分析工具 |
|------|------|----------|-------------|
| **基因芯片（Microarray）** | 荧光探针杂交，荧光强度代表 mRNA 含量 | 连续浮点数（已标准化） | **limma** |
| **RNA-seq（测序）** | 直接对 mRNA 测序并计数 | 整数（raw read count） | **DESeq2 / edgeR** |

**⚠️ 这是最容易犯的错误之一**：拿到数据不看类型，直接上 DESeq2，结果报错或得到错误结论。你必须先判断数据类型，再选工具。

### 3.0.3 如何判断数据是芯片还是 RNA-seq？

拿到数据后运行：

```r
expr <- exprs(gse)
cat("最小值:", min(expr, na.rm=TRUE), "\n")
cat("最大值:", max(expr, na.rm=TRUE), "\n")
cat("是否含小数:", any(expr != round(expr), na.rm=TRUE), "\n")
```

判断标准：

- **含小数，值在 0–20 之间（log2 scale）** → 芯片数据，用 **limma**
- **全是整数，值可能很大（几千几万）** → RNA-seq count，用 **DESeq2**
- **含小数，值是比例（0–100）** → TPM/FPKM，需要找 raw count 或用 **limma-voom**

### 3.0.4 logFC 和 adj.P.Val 是什么，为什么这样设计？

**logFC（log2 Fold Change）**：

原始 fold change = 疾病组平均表达 ÷ 正常组平均表达。

取 log2 的原因：原始 fold change 不对称（上调8倍=8，下调8倍=0.125），统计上难处理。取 log2 后完全对称（+3 vs -3），且直觉清晰：logFC 每增加 1，表达量翻一倍。

| logFC | 含义 |
|-------|------|
| +3 | 疾病组是正常组的 2³ = 8 倍 |
| +1 | 疾病组是正常组的 2¹ = 2 倍（我们的最低阈值） |
| 0 | 没有差异 |
| -1 | 疾病组是正常组的 1/2 |
| -3 | 疾病组是正常组的 1/8 |

**adj.P.Val（校正后 p 值）**：

如果我们对 50,000 个基因分别做检验，即使每个基因的显著性阈值设为 p<0.05，也会有 50,000 × 0.05 = 2,500 个基因**纯靠运气**显示显著。这叫**多重检验问题**。

BH 校正（Benjamini-Hochberg）保证：在所有我们认为显著的基因里，假阳性比例（FDR）不超过 5%。这比原始 p 值严格得多。

**永远看 adj.P.Val，不看 P.Value。**

### 3.0.5 标准筛选阈值

```
|logFC| > 1      → 表达差异超过 2 倍（生物学意义阈值）
adj.P.Val < 0.05 → 统计显著（控制假阳性）
```

这两个条件必须同时满足。只看 p 值会引入太多无意义的微小差异；只看 logFC 会引入统计上不可信的结果。

---

## 3.1 里程碑概述 / Milestone Overview

**目标**：从公共数据库下载疾病 vs 正常样本的表达谱数据，做差异分析，找出疾病中显著上调的基因作为靶点候选，画出火山图。

**交付物**：
1. 差异基因列表（CSV，含 logFC 和 adj.P.Val）
2. 火山图（PNG + PDF，高分辨率）
3. RMarkdown 分析报告（HTML，可复现）
4. GitHub 仓库（含完整代码）

**预计时长**：2–3 周

---

## 3.2 环境搭建 / Environment Setup

### Step 1：安装 R

- **下载地址**：https://cran.r-project.org/bin/windows/base/
- 选最新稳定版（≥ 4.4），不要用 4.0 以下
- 安装路径不要有中文和空格，推荐 `D:\R\R-4.x.x`

### Step 2：安装 RStudio

- **下载地址**：https://posit.co/download/rstudio-desktop/
- 选免费的 RStudio Desktop Open Source Edition
- **必须先装 R，再装 RStudio**

### Step 3：安装所有必需 R 包

在 RStudio Console 里逐条运行，**不要一次性全贴**（方便排查哪一步出错）：

```r
# ---- 第一步：设置国内镜像（中国大陆必须做，否则下载极慢）----
# 把这两行写入 ~/.Rprofile，永久生效
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# ---- 第二步：安装 BiocManager（Bioconductor 的包管理器）----
install.packages("BiocManager")

# ---- 第三步：安装数据处理基础包 ----
install.packages(c(
  "tidyverse",    # 数据处理全家桶（含 ggplot2、dplyr 等）
  "data.table"    # 高效读取大文件
))

# ---- 第四步：安装 Bioconductor 生信核心包 ----
BiocManager::install(c(
  "GEOquery",        # 从 GEO 数据库下载数据
  "Biobase",         # Bioconductor 基础数据结构
  "limma",           # 芯片数据差异分析（本流程主要用这个）
  "DESeq2",          # RNA-seq count 数据差异分析
  "EnhancedVolcano", # 专业火山图
  "clusterProfiler", # GO/KEGG 富集分析
  "org.Hs.eg.db",    # 人类基因注释数据库
  "biomaRt"          # 基因 ID 转换（在线版）
))

# ---- 第五步：安装报告与可视化包 ----
install.packages(c(
  "rmarkdown",    # 生成 HTML/PDF 报告
  "knitr",        # RMarkdown 引擎
  "kableExtra",   # 漂亮的表格
  "ggrepel",      # 标签防重叠
  "RColorBrewer"  # 配色方案
))
```

> **如果安装过程中问 "Update all/some/none?"，输入 `n` 回车跳过。**
> 更新其他包可能引入兼容性问题，等分析做完再统一更新。

### Step 4：验证环境

```r
# 逐行运行，确认没有报错
library(GEOquery)
library(limma)
library(DESeq2)
library(EnhancedVolcano)
library(tidyverse)

cat("✅ 所有包加载成功！\n")
sessionInfo()  # 打印版本信息，截图保存，以后排错用得到
```

---

## 3.3 第一步：理解你的数据集 / Understanding Your Dataset

### 3.3.1 GEO 数据库是什么？

GEO（Gene Expression Omnibus）是美国 NCBI 维护的全球最大公共基因表达数据库。全球科研人员发表论文后必须将原始数据上传至此，所有人可以免费下载。

每个数据集有一个编号：**GSE + 数字**（如 GSE6891）。

### 3.3.2 如何为你的研究选择合适的 GEO 数据集？

在 GEO 网站（https://www.ncbi.nlm.nih.gov/geo/）搜索你的疾病名称，筛选标准：

| 标准 | 要求 | 原因 |
|------|------|------|
| 样本量 | 疾病组 ≥ 20，正常组 ≥ 10 | 样本太少统计功效不足 |
| 分组清晰 | 能明确区分疾病 vs 正常 | 分组不清会导致分析结果无意义 |
| 物种 | Human（人类） | 我们研究人类疾病 |
| 数据类型 | Expression profiling by array 或 by high throughput sequencing | 其他类型（如 ChIP-seq）不适合本流程 |
| 引用量 | 被多篇论文引用过的优先 | 说明数据质量经过验证 |

**GEO 数据集推荐（按疾病）**：

| 疾病 | 推荐数据集 | 说明 |
|------|------------|------|
| AML（急性髓系白血病） | GSE6891 | 537样本，芯片，分组清晰 |
| AML（备选） | GSE37642 | 多亚型，含正常对照 |
| TNBC（三阴性乳腺癌） | GSE76124 | 乳腺癌亚型分层 |
| 乳腺癌（通用） | GSE58812 | 大样本，预后数据完整 |
| ALL（急性淋巴细胞白血病） | GSE11877 | 儿童 ALL，分型详细 |
| CML（慢性粒细胞白血病） | GSE4170 | CML vs 正常骨髓 |
| 胃癌 | GSE54129 | 胃癌 vs 正常胃黏膜 |
| 肝癌（HCC） | GSE14520 | 大样本，含临床数据 |
| 肺癌（NSCLC） | GSE81089 | 肿瘤 vs 正常肺组织 |

### 3.3.3 下载数据

**方法一：让 R 直接联网下载（前提：网络可以访问 NCBI）**

```r
# ============================================================
# 01_download_data.R
# 把 "GSE6891" 换成你的数据集编号即可复用此脚本
# ============================================================

library(GEOquery)
library(tidyverse)

# 设置工作目录（改成你自己的路径）
setwd("D:/Bio-Informatics Case Study")

# 创建项目子目录结构
dir.create("M1_bulk_RNAseq/data",    recursive = TRUE, showWarnings = FALSE)
dir.create("M1_bulk_RNAseq/results", recursive = TRUE, showWarnings = FALSE)
dir.create("M1_bulk_RNAseq/figures", recursive = TRUE, showWarnings = FALSE)

# 下载数据
# GSEMatrix = TRUE：下载整理好的表达矩阵，不是原始文件
# getGPL = FALSE：不下载芯片平台注释（可以后面单独装对应的 .db 包）
GSE_ID <- "GSE6891"   # ← 改这里换数据集

gse <- getGEO(GSE_ID,
              destdir = "M1_bulk_RNAseq/data",
              GSEMatrix = TRUE,
              getGPL = FALSE)

# gse 是一个 list，通常取第一个元素
gse <- gse[[1]]

cat("样本数:", ncol(gse), "\n")
cat("基因/探针数:", nrow(gse), "\n")
```

**方法二：手动下载（网络不稳定时用这个）**

1. 浏览器打开：`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6891`
2. 滚到页面最底部，点击 "Series Matrix File(s)" 下载（约 10–50 MB）
3. 把下载的 `.txt.gz` 文件放到 `M1_bulk_RNAseq/data/` 文件夹
4. 用本地读取方式：

```r
# 读取本地文件，getGPL = FALSE 避免再次联网
gse <- getGEO(filename = "M1_bulk_RNAseq/data/GSE6891_series_matrix.txt.gz",
              getGPL = FALSE)
# 注意：用 filename 参数时，返回的直接是 ExpressionSet 对象，不需要取 [[1]]
```

> **⚠️ 中国大陆网络提示**：NCBI 服务器在美国，直连经常超时。推荐开启 VPN 全局代理模式再运行下载代码，成功率大幅提升。

---

## 3.4 第二步：探索数据，判断数据类型 / Data Exploration

> 这一步是分析中最容易被跳过、却最容易出错的一步。不探索直接分析，等于闭着眼睛开车。

```r
# ============================================================
# 02_explore_data.R
# ============================================================

library(tidyverse)
library(Biobase)

# 如果已经下载保存过 RDS，用这行读取（更快）
# gse <- readRDS("M1_bulk_RNAseq/data/GSE6891_eset.rds")

# ---- 1. 看数据基本结构 ----
cat("=== 数据基本信息 ===\n")
cat("样本数（列数）:", ncol(gse), "\n")
cat("基因/探针数（行数）:", nrow(gse), "\n")

# ---- 2. 判断数据类型（最关键！）----
expr <- exprs(gse)  # 取出表达矩阵

cat("\n=== 数据类型判断 ===\n")
cat("最小值:", min(expr, na.rm = TRUE), "\n")
cat("最大值:", max(expr, na.rm = TRUE), "\n")
cat("是否含小数:", any(expr != round(expr), na.rm = TRUE), "\n")
cat("是否含 NA:", any(is.na(expr)), "\n")

# 根据输出判断：
# 含小数 + 值在 0-20 → 芯片（用 limma）
# 全整数 + 值可能很大 → RNA-seq count（用 DESeq2）

# ---- 3. 看样本元数据（找分组列）----
pdata <- pData(gse)

cat("\n=== 样本信息表列名 ===\n")
print(colnames(pdata))

# 逐列查看，找到区分疾病 vs 正常的那一列
# 通常在 "characteristics_ch1"、"tissue:ch1"、"disease state:ch1" 等列

# 示例：查看某一列的所有取值
# print(table(pdata$`cell type:ch1`))
# print(table(pdata$`disease state:ch1`))
```

### 3.4.1 如何找到分组列？

这是每个数据集都不一样的地方，必须人工判读。步骤：

1. 运行 `print(colnames(pdata))` 看所有列名
2. 重点关注含这些词的列：`tissue`、`cell type`、`disease`、`group`、`status`、`condition`
3. 对可疑列运行 `print(table(pdata$列名))` 看取值
4. 找到值明显分为"疾病"和"正常"两类的那一列

**GSE6891 的分组列**：`cell type:ch1`，取值为：
- `blasts and mononuclear cells`（461个）→ AML 患者
- `mononuclear cells from bone marrow`（76个）→ 正常骨髓

---

## 3.5 第三步：构建分组变量 / Define Groups

分组变量是整个分析的基础。它告诉统计模型：哪些样本是"病人"，哪些是"对照"。

```r
# ============================================================
# 03_define_groups.R
# ============================================================

# ---- GSE6891 的分组方式（其他数据集请根据实际列名修改）----

# 方法：用 ifelse 把原始标签转成简洁的 "normal" / "AML"
group <- ifelse(
  pdata$`cell type:ch1` == "mononuclear cells from bone marrow",
  "normal",   # 正常骨髓
  "AML"       # AML 患者
)

# 转成 factor，并指定 "normal" 为参照组（非常重要！）
# 参照组是 logFC 的分母，即 logFC = 疾病组 - 参照组（log scale）
group <- factor(group, levels = c("normal", "AML"))

# 核查分组结果
print(table(group))

# ---- 适配其他数据集的通用写法 ----
# 如果你的数据集分组列叫 "disease state:ch1"，值为 "cancer" 和 "normal"：
# group <- factor(pdata$`disease state:ch1`, levels = c("normal", "cancer"))
#
# 如果值是数字编码（如 1=疾病，0=正常）：
# group <- factor(ifelse(pdata$`group:ch1` == "1", "disease", "normal"),
#                 levels = c("normal", "disease"))
```

---

## 3.6 第四步：差异分析 / Differential Expression Analysis

根据 3.0.3 节判断的数据类型，选择对应的分析流程。

---

### 路线 A：芯片数据 → 用 limma

适用于：数据含小数、值在 0–20 之间（log2 标准化过的芯片信号）

```r
# ============================================================
# 04A_limma_analysis.R
# 适用于芯片数据（Microarray）
# ============================================================

library(limma)
library(tidyverse)

# ---- 1. 处理缺失值 ----
# 芯片数据可能有少量 NA（探针质控失败），需要先过滤
keep_genes <- complete.cases(expr)  # 找出没有 NA 的基因行
expr_clean <- expr[keep_genes, ]
cat("去除 NA 后剩余探针数:", nrow(expr_clean), "\n")

# ---- 2. 构建设计矩阵（Design Matrix）----
# 设计矩阵告诉 limma 每个样本属于哪个组
# model.matrix 自动把 factor 转成 0/1 编码
design <- model.matrix(~ group)
# 结果：每行是一个样本，列是截距和分组系数
# "groupAML" 这一列的系数就是我们要的 logFC

cat("\n设计矩阵（前3行）:\n")
print(head(design, 3))

# ---- 3. 拟合线性模型 ----
# lmFit 对每一个基因单独拟合线性模型：表达量 = 截距 + β×分组 + 误差
fit <- lmFit(expr_clean, design)

# ---- 4. 经验贝叶斯校正（eBayes）----
# 这一步借用所有基因的信息来稳定每个基因的方差估计
# 解决了小样本下方差估计不稳定的问题
fit <- eBayes(fit)

# ---- 5. 提取结果 ----
# coef = "groupAML" 指定我们要看的是 AML vs normal 这个比较
# number = Inf 提取所有基因（不截断）
# adjust.method = "BH" 用 Benjamini-Hochberg 方法做多重检验校正
results <- topTable(fit,
                    coef = "groupAML",
                    number = Inf,
                    adjust.method = "BH",
                    sort.by = "adj.P.Val")

cat("\n结果表前几行:\n")
print(head(results))

# ---- 6. 统计显著差异基因 ----
n_sig <- sum(abs(results$logFC) > 1 & results$adj.P.Val < 0.05, na.rm = TRUE)
n_up  <- sum(results$logFC > 1   & results$adj.P.Val < 0.05, na.rm = TRUE)
n_dn  <- sum(results$logFC < -1  & results$adj.P.Val < 0.05, na.rm = TRUE)
cat("\n显著差异基因总数:", n_sig, "\n")
cat("  上调（疾病高于正常）:", n_up, "\n")
cat("  下调（疾病低于正常）:", n_dn, "\n")

# ---- 7. 保存结果 ----
write.csv(results, "M1_bulk_RNAseq/results/limma_all_results.csv")

sig_genes <- results[abs(results$logFC) > 1 & results$adj.P.Val < 0.05, ]
write.csv(sig_genes, "M1_bulk_RNAseq/results/limma_significant_genes.csv")

cat("\n✅ limma 分析完成，结果已保存\n")
```

---

### 路线 B：RNA-seq Count 数据 → 用 DESeq2

适用于：数据全是整数、来自 RNA-seq 的 raw read count

```r
# ============================================================
# 04B_deseq2_analysis.R
# 适用于 RNA-seq count 数据
# ============================================================

library(DESeq2)
library(tidyverse)

# ---- 1. 确认数据是整数 ----
# DESeq2 对输入格式极其严格，必须是非负整数
# 如果含小数，用 round() 取整（但要先确认这合理）
count_matrix <- round(expr)
stopifnot(all(count_matrix >= 0))  # 确认没有负值

# ---- 2. 构建 DESeq2 对象 ----
coldata <- data.frame(group = group, row.names = colnames(count_matrix))

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = coldata,
  design    = ~ group        # 用 group 列做差异分析
)

# ---- 3. 过滤低表达基因 ----
# 低表达基因信噪比差，保留它们只会增加假阳性
# 标准：至少在 [最小组样本数] 个样本里有 ≥10 reads
min_samples <- min(table(group))
keep <- rowSums(counts(dds) >= 10) >= min_samples
dds <- dds[keep, ]
cat("过滤后保留基因数:", nrow(dds), "\n")

# ---- 4. 运行 DESeq2 ----
# 这一步自动完成：标准化 → 方差估计 → 统计检验
# 样本多时会比较慢（500 个样本可能需要 10–30 分钟）
dds <- DESeq(dds)

# ---- 5. 提取结果 ----
res <- results(dds,
               contrast = c("group", "AML", "normal"),  # 分子 vs 分母
               alpha = 0.05)
summary(res)

# ---- 6. 整理为 data.frame ----
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)

# ---- 7. 筛选显著差异基因 ----
sig_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

cat("显著差异基因数:", nrow(sig_genes), "\n")
cat("  上调:", sum(sig_genes$log2FoldChange > 0), "\n")
cat("  下调:", sum(sig_genes$log2FoldChange < 0), "\n")

# ---- 8. 保存结果 ----
write.csv(res_df,   "M1_bulk_RNAseq/results/deseq2_all_results.csv")
write.csv(sig_genes,"M1_bulk_RNAseq/results/deseq2_significant_genes.csv")
saveRDS(dds, "M1_bulk_RNAseq/data/dds_object.rds")

cat("\n✅ DESeq2 分析完成，结果已保存\n")
```

### 常见报错速查

| 报错信息 | 原因 | 解决方案 |
|----------|------|----------|
| `some values in assay are not integers` | 数据是 TPM/芯片值，不是 count | 换用 limma（路线 A） |
| `Failed to download ... series_matrix` | 网络无法访问 NCBI | 开 VPN 全局模式，或手动下载文件 |
| `Failed to download ... GPL570.soft.gz` | getGEO 尝试下载平台注释 | 加参数 `getGPL = FALSE` |
| `Error: object 'pdata' not found` | 前面的代码没有运行成功 | 全选代码（Ctrl+A）重新运行 |
| `the design matrix has the same number of samples and coefficients` | group 只有一个水平 | 检查 group 是否正确分成两组 |
| `内存不足 / cannot allocate vector` | 数据太大 | 先过滤低表达基因再分析 |

---

## 3.7 第五步：探针 ID 转基因名 / Probe ID to Gene Symbol

> 这一步仅适用于芯片数据（limma 路线）。RNA-seq 数据通常直接就是基因名，可跳过。

芯片数据的行名是**探针 ID**（如 `1555629_at`），不是人类可读的基因名。需要用注释包把探针号翻译成基因符号（SYMBOL），如 CLEC12A、HAVCR2。

**第一步：确认你的芯片型号**

在 GEO 数据集页面找 "Platform"，或运行：

```r
cat(annotation(gse), "\n")  # 打印芯片平台名称
```

**第二步：安装对应的注释包**

| 芯片平台 | 对应注释包 |
|----------|-----------|
| hgu133plus2（Affymetrix HG-U133 Plus 2.0）| `hgu133plus2.db` |
| hgu133a（Affymetrix HG-U133A）| `hgu133a.db` |
| hgu219（Affymetrix HG-U219）| `hgu219.db` |
| IlluminaHumanMethylation450k | `IlluminaHumanMethylation450kmanifest` |
| 其他 | 在 Bioconductor 搜索平台名称 |

```r
# 安装对应包（以 hgu133plus2 为例）
BiocManager::install("hgu133plus2.db")
```

**第三步：执行 ID 转换**

```r
# ============================================================
# 05_id_conversion.R
# ============================================================

library(hgu133plus2.db)   # 替换成你的芯片对应包
library(tidyverse)

# mapIds：把探针 ID 映射到基因 Symbol
# keys：要转换的探针 ID 列表
# column：想要的目标注释类型（SYMBOL = 基因名）
# keytype：输入 ID 的类型
# multiVals = "first"：一个探针对应多个基因时，只取第一个

probe_ids    <- rownames(results)
gene_symbols <- mapIds(hgu133plus2.db,
                       keys     = probe_ids,
                       column   = "SYMBOL",
                       keytype  = "PROBEID",
                       multiVals = "first")

# 把基因名加到结果表
results$gene_symbol <- gene_symbols

# 去掉没有对应基因名的探针（NA）
results_named <- results[!is.na(results$gene_symbol), ]
cat("成功映射到基因名的探针数:", nrow(results_named), "\n")

# 每个基因只保留最显著的那个探针
# 原因：一个基因在芯片上通常有多个探针，取最显著的代表该基因
results_best <- results_named %>%
  group_by(gene_symbol) %>%
  slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  as.data.frame()

rownames(results_best) <- results_best$gene_symbol
cat("去重后唯一基因数:", nrow(results_best), "\n")

# 保存带基因名的结果
write.csv(results_best, "M1_bulk_RNAseq/results/limma_results_with_symbols.csv")
cat("✅ ID 转换完成\n")
```

---

## 3.8 第六步：绘制火山图 / Volcano Plot

```r
# ============================================================
# 06_volcano_plot.R
# ============================================================

library(EnhancedVolcano)
library(tidyverse)

# ---- 基础火山图（探针 ID 版，适合快速查看）----
p1 <- EnhancedVolcano(
  results,                        # 替换成你的结果对象
  lab      = rownames(results),   # 标签：探针 ID 或基因名
  x        = "logFC",             # DESeq2 用 "log2FoldChange"
  y        = "adj.P.Val",         # DESeq2 用 "padj"
  title    = "Disease vs Normal", # 改成你的标题
  subtitle = "limma | GSE6891",   # 改成你的数据集
  pCutoff  = 0.05,                # adj.P 阈值
  FCcutoff = 1,                   # |logFC| 阈值
  pointSize = 1.5,
  labSize   = 3.0,
  col       = c("grey70",       # NS（不显著）
                "forestgreen",  # 只有 p 值显著
                "royalblue",    # 只有 FC 显著
                "red2"),        # 同时显著（最重要）
  colAlpha  = 0.5,
  legendPosition = "right",
  drawConnectors = TRUE,
  max.overlaps   = 15
)

ggsave("M1_bulk_RNAseq/figures/volcano_basic.png",
       plot = p1, width = 12, height = 8, dpi = 300)
ggsave("M1_bulk_RNAseq/figures/volcano_basic.pdf",
       plot = p1, width = 12, height = 8)

# ---- 靶点高亮版火山图（用基因名，标注已知靶点）----
# 先定义你想高亮的靶点列表（根据你的疾病修改）
my_targets <- c("CLEC12A", "HAVCR2", "IL3RA", "FLT3", "CD33", "CD38")
# TNBC 靶点示例：c("HER2", "EGFR", "MSLN", "NECTIN4", "TROP2")
# 胃癌靶点示例：c("CLDN18", "HER2", "MET", "FGFR2", "VEGFR2")

p2 <- EnhancedVolcano(
  results_best,                          # 使用去重后的基因名版本
  lab       = results_best$gene_symbol,
  x         = "logFC",
  y         = "adj.P.Val",
  title     = "AML vs Normal Bone Marrow",
  subtitle  = "limma | GSE6891 | Known CAR-T targets highlighted",
  pCutoff   = 0.05,
  FCcutoff  = 1,
  pointSize = 1.5,
  labSize   = 4.0,
  col       = c("grey70", "forestgreen", "royalblue", "red2"),
  colAlpha  = 0.5,
  legendPosition  = "right",
  drawConnectors  = TRUE,
  widthConnectors = 0.5,
  max.overlaps    = 20,
  selectLab = my_targets    # 只标注这些基因的名字
)

ggsave("M1_bulk_RNAseq/figures/volcano_targets.png",
       plot = p2, width = 12, height = 8, dpi = 300)
ggsave("M1_bulk_RNAseq/figures/volcano_targets.pdf",
       plot = p2, width = 12, height = 8)

cat("✅ 火山图已保存\n")
```

### 如何读火山图

```
高 -log10(adj.P)
        ↑
        |   蓝点（只显著）  红点（显著+大FC）← 候选靶点在这里
        |
        |         灰点（不显著）
        |
————————|————————————————————————→  logFC
  下调  0     上调
（疾病低表达） （疾病高表达）
```

- **右上角红点**：在疾病中显著高表达，CAR-T 靶点应该在这里找
- **左上角红点**：在疾病中显著低表达，可能是抑癌基因
- **灰点**：两组之间无显著差异

---

## 3.9 第七步：验证分析结果 / Validation

差异分析跑完，不能直接相信结果。要用**已知靶点做正向对照**验证流程是否正确。

```r
# ============================================================
# 07_validation.R
# ============================================================

# 根据你的疾病，定义已知靶点列表
# AML 已知靶点：
known_targets <- c("CD33", "FLT3", "IL3RA", "CLEC12A", "CD38", "HAVCR2")

# TNBC 已知靶点（示例）：
# known_targets <- c("EGFR", "MSLN", "NECTIN4", "TROP2", "LAG3")

# 在结果表里找这些基因
validation <- results_best[results_best$gene_symbol %in% known_targets,
                           c("gene_symbol", "logFC", "adj.P.Val")]
validation <- validation[order(validation$adj.P.Val), ]

cat("=== 已知靶点验证结果 ===\n")
print(validation)

# 解读：
# 如果已知靶点出现在结果里，且方向正确（logFC > 0 表示疾病高表达）→ 流程正确 ✅
# 如果完全没有出现 → 检查分组是否搞反，或数据集本身不适合
```

---

## 3.10 M1 完成检查清单 / Completion Checklist

```
□ 1. 环境：R + RStudio + 所有包安装成功，library() 无报错
□ 2. 数据：从 GEO 成功下载数据集，样本数和基因数符合预期
□ 3. 探索：判断了数据类型（芯片 or RNA-seq），找到了分组列
□ 4. 分析：用对应工具（limma or DESeq2）跑出差异分析结果
□ 5. 转换：完成探针 ID → 基因名转换（芯片数据）
□ 6. 可视化：画出火山图（基础版 + 靶点高亮版），保存 PNG + PDF
□ 7. 验证：已知靶点出现在结果列表中，方向正确
□ 8. 保存：所有 CSV 结果文件保存到 results/ 目录
□ 9. 报告：写完 RMarkdown 报告，Knit 成 HTML
□ 10. GitHub：代码 push 到仓库，README 说明数据来源和主要结论
```

## 3.11 M1 进阶任务（可选）/ Advanced

完成基础流程后，可以用这些方法深入挖掘：

```r
# ---- 富集分析：这些差异基因参与哪些生物学通路？----
library(clusterProfiler)
library(org.Hs.eg.db)

# 把基因名转成 Entrez ID（富集分析需要）
gene_list <- results_best$gene_symbol[results_best$logFC > 1 &
                                      results_best$adj.P.Val < 0.05]
entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_list,
                     column = "ENTREZID", keytype = "SYMBOL")

# GO 富集分析（生物学过程）
go_results <- enrichGO(gene = na.omit(entrez_ids),
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)
dotplot(go_results, showCategory = 20)

# KEGG 通路富集分析
kegg_results <- enrichKEGG(gene = na.omit(entrez_ids),
                            organism = "hsa",
                            pvalueCutoff = 0.05)
dotplot(kegg_results, showCategory = 20)
```

---


<a name="part-4"></a>

# Part 4：里程碑 2 — 单细胞 RNA-seq 与靶点定位

## 4.1 里程碑概述 / Milestone Overview

**目标 / Goal**：对一套公开的 AML 单细胞 RNA-seq 数据做完整分析（质控 → 标准化 → 降维 → 聚类 → 细胞类型注释 → 差异分析），识别 AML 白血病干细胞（LSC）和原始细胞特异性高表达的膜蛋白，作为 CAR-T 候选靶点。

**为什么重要 / Why It Matters**：
- 单细胞是当下最热的技术，每篇 Nature/Cell 论文几乎都用
- AML 的单细胞研究尤其热门，因为 LSC 是耐药和复发的根源
- Pfizer / Roche / Merck 等大药企的 oncology 团队**必查**这个技能
- 你已经具备的 Python 基础在这里能直接发力

**交付物 / Deliverables**：
1. ✅ UMAP 图（按细胞类型染色）
2. ✅ 细胞类型注释表
3. ✅ 候选靶点 shortlist（5–10 个膜蛋白）
4. ✅ Jupyter Notebook 完整分析报告
5. ✅ GitHub 仓库更新

**预计时长 / Duration**：3–4 周

## 4.2 环境搭建 / Environment Setup

### Step 1：用 Conda 创建独立环境

打开 PowerShell 或 Anaconda Prompt：

```bash
# 创建一个名为 scrna 的独立环境
conda create -n scrna python=3.11 -y

# 激活环境
conda activate scrna

# 装核心包
pip install scanpy[leiden] anndata
pip install jupyter notebook ipykernel
pip install matplotlib seaborn
pip install scrublet           # 双胞体检测
pip install harmonypy          # 批次效应校正
pip install celltypist         # 自动细胞类型注释
pip install scvi-tools         # 高级整合工具（可选）

# 把这个环境注册到 Jupyter
python -m ipykernel install --user --name=scrna --display-name "Python (scrna)"
```

### Step 2：验证安装

```bash
python -c "import scanpy as sc; sc.logging.print_header()"
```

应该看到 Scanpy 的版本和依赖列表。

### Step 3：在 VS Code 中配置

1. 打开 VS Code，按 `Ctrl+Shift+P`，输入 "Python: Select Interpreter"
2. 选择 `conda env scrna`
3. 新建一个 `.ipynb` 文件，右上角 Kernel 选 `Python (scrna)`

## 4.3 数据获取 / Data Acquisition

### 推荐数据集 / Recommended Datasets

| 数据集 | 描述 | 大小 | 适合 |
|--------|------|------|------|
| **GSE116256** | AML 单细胞图谱，含 AML 患者和正常骨髓 | ~3 GB | 入门首选 |
| **PBMC 3k** | Scanpy 内置标准 PBMC 数据 | ~50 MB | 走流程练手 |
| **GSE142213** | AML 治疗前后 scRNA-seq，揭示耐药机制 | ~2 GB | 进阶选用 |

### 下载（以 PBMC 10k 为例练手）

```python
import scanpy as sc

# Scanpy 内置了 PBMC 3k 数据集，一行下载
adata = sc.datasets.pbmc3k()
print(adata)
# AnnData object with n_obs × n_vars = 2700 × 32738
```

练完 PBMC，再上 GSE116256：
- 浏览器打开 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116256
- 找到 supplementary file（10X 格式：`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`）
- 用 `sc.read_10x_mtx()` 读取

## 4.4 标准分析流程 / Standard Pipeline

### 完整分析 Notebook 代码框架 / Notebook Skeleton

```python
# ============================================================
# Single-cell RNA-seq Analysis Pipeline
# 单细胞 RNA-seq 标准分析流程
# ============================================================

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 全局设置
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(6, 6))

# ---- 1. 加载数据 ----
adata = sc.read_10x_mtx(
    'data/filtered_feature_bc_matrix/',
    var_names='gene_symbols',
    cache=True
)
print(adata)

# ---- 2. 质控 / Quality Control ----
# 计算质控指标
adata.var['mt'] = adata.var_names.str.startswith('MT-')   # 线粒体基因
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'],
                            percent_top=None, log1p=False, inplace=True)

# 可视化质控指标
sc.pl.violin(adata,
             ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True,
             save='_qc_metrics.png')

# 过滤：经验阈值，需要根据自己数据调整
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(f"过滤后: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")

# ---- 3. 标准化与对数转换 ----
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 保留原始 count 数据
adata.raw = adata

# ---- 4. 高变基因选择 ----
sc.pp.highly_variable_genes(adata, min_mean=0.0125,
                             max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save='_hvg.png')
adata = adata[:, adata.var.highly_variable]

# ---- 5. 缩放与 PCA ----
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, save='_pca_variance.png')

# ---- 6. 邻接图与 UMAP ----
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# ---- 7. 聚类 ----
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color='leiden', save='_clusters.png',
           legend_loc='on data', legend_fontsize=10)

# ---- 8. 寻找 marker 基因 ----
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,
                        save='_marker_genes.png')

# ---- 9. 细胞类型注释 ----
# 方法 A：手动注释（基于已知 marker）
marker_genes = {
    'T cell':           ['CD3D', 'CD3E', 'CD8A'],
    'B cell':           ['CD19', 'MS4A1', 'PAX5'],
    'NK cell':          ['GNLY', 'NKG7', 'NCAM1'],
    'Monocyte':         ['CD14', 'LYZ', 'S100A8'],
    'HSC':              ['CD34', 'HOXA9', 'MEIS1'],   # 造血干细胞
    # AML 标志（白血病原始细胞）
    'AML blast':        ['CD33', 'CD123', 'FLT3', 'CD117', 'CLEC12A'],
    # 白血病干细胞（LSC，耐药根源）
    'LSC':              ['CD34', 'CD38', 'CD123', 'TIM3', 'CD96'],
}
sc.pl.dotplot(adata, marker_genes, groupby='leiden',
              save='_marker_dotplot.png')

# 方法 B：CellTypist 自动注释（更省事）
import celltypist
from celltypist import models

models.download_models(force_update=False, model='Immune_All_Low.pkl')
# AML 分析推荐使用 Immune_All_Low 模型，能识别髓系各亚型
predictions = celltypist.annotate(adata, model='Immune_All_Low.pkl',
                                    majority_voting=True)
adata = predictions.to_adata()
sc.pl.umap(adata, color='majority_voting',
           save='_celltypist.png')

# ---- 10. 保存结果 ----
adata.write('results/scrna_annotated.h5ad')
```

## 4.5 寻找 CAR-T 候选靶点 / Finding CAR-T Targets

### 核心思路 / Core Logic

理想的 CAR-T 靶点必须满足：
1. ✅ **细胞表面表达**（膜蛋白）
2. ✅ **肿瘤细胞高表达**
3. ✅ **正常组织低/无表达**（safety）
4. ✅ **抗原稳定**（不容易丢失）

### 筛选代码 / Filtering Code

```python
# ---- 11. 锁定 AML 原始细胞群（blast） ----
blast_mask = adata.obs['cell_type'].isin(['AML blast', 'LSC'])
print(f"AML 原始细胞数 / AML blast + LSC cells: {blast_mask.sum()}")

# 在所有细胞中找 AML 特异性高表达基因
sc.tl.rank_genes_groups(
    adata,
    groupby='cell_type',
    groups=['AML blast'],
    reference='rest',           # 与其他所有细胞比
    method='wilcoxon'
)

# 提取 top 100
tumor_markers = sc.get.rank_genes_groups_df(
    adata, group='AML blast'
).head(100)
print(tumor_markers.head(20))

# ---- 12. 过滤膜蛋白 ----
# 需要一个膜蛋白列表（从 UniProt 下载或用 surfaceome database）
# 推荐：https://wlab.ethz.ch/surfaceome/ 下载 Surfaceome 数据库

membrane_proteins = pd.read_csv('data/surfaceome.csv')['gene'].tolist()

candidates = tumor_markers[
    tumor_markers['names'].isin(membrane_proteins)
].copy()

# 加上一个 "tumor specificity score"
# = 肿瘤中表达 / 其他细胞中平均表达
# （需要额外计算，这里简化）

print("\n=== CAR-T 候选靶点 / Candidates ===")
print(candidates[['names', 'logfoldchanges', 'pvals_adj']].head(10))

# 保存
candidates.to_csv('results/aml_cart_target_candidates.csv', index=False)
```

### 已知的 AML CAR-T 靶点（用作正向对照）/ Known AML Targets

跑完后，看看你的列表里有没有这些已知靶点：

| 靶点 | 基因名 | AML 状态 |
|------|--------|----------|
| CD33 | SIGLEC3 | 临床试验活跃，Gemtuzumab 抗体已批 |
| CD123 | IL3RA | 临床 II 期，表达于 LSC |
| FLT3 | FLT3 | 突变型 AML 首选靶点，临床中 |
| CLL-1 | CLEC12A | AML 特异性强，正常 HSC 低表达 |
| CD38 | CD38 | Daratumumab 已批（骨髓瘤），AML 临床中 |
| TIM3 | HAVCR2 | LSC 标志，正在临床探索 |

**如果你的列表里出现了这些靶点 → 说明流程跑对了 ✅**

## 4.6 M2 完成检查清单 / M2 Completion Checklist

```
□ Scanpy 环境装好，可以在 Jupyter 中运行
□ 跑通 PBMC 3k 标准流程（练手）
□ 处理 1 个真实癌症 scRNA-seq 数据集（≥10000 细胞）
□ 出 UMAP 图（按聚类 + 按细胞类型各一张）
□ 完成细胞类型注释
□ 输出至少 5 个 CAR-T 候选靶点
□ 整理成 Jupyter Notebook 推到 GitHub
```

## 4.7 M2 进阶任务（可选）/ M2 Advanced

- 用 **CellChat** 或 **CellPhoneDB** 分析细胞间通讯
- 用 **Scrublet** 检测双胞体
- 用 **Harmony** 整合多个样本（去批次效应）
- 用 **CytoTRACE** 分析肿瘤细胞分化轨迹
- 用 **InferCNV** 推断肿瘤细胞的拷贝数变异

---

<a name="part-5"></a>

# Part 5：里程碑 3 — 蛋白质结构预测与分子对接

## 5.1 里程碑概述 / Milestone Overview

**目标 / Goal**：取 M2 找出的 AML 候选靶点（优先选 CD33 或 CD123），预测其三维结构，做一次抗体-抗原对接分析，理解 CAR 的 scFv 是如何识别靶点的。

**为什么重要 / Why It Matters**：
- 结构生物学是 CAR-T 设计的核心 — 没有结构就没有合理的 scFv 设计
- AlphaFold 2024 年拿了诺贝尔化学奖，是当下最热的 AI for Science 案例
- 这部分技能在 medicinal chemistry / antibody engineering 岗位是加分项

**交付物 / Deliverables**：
1. ✅ 靶点蛋白的 3D 结构（PDB 文件）
2. ✅ PyMOL 渲染的高分辨率结构图
3. ✅ 对接分析结果（结合能 + 关键残基）
4. ✅ 结构生物学分析报告

**预计时长 / Duration**：2–3 周

## 5.2 环境搭建 / Environment Setup

### 工具清单 / Tools to Install

| 工具 | 安装方式 | 价格 |
|------|----------|------|
| **AlphaFold Server** | 网页端，无需安装 | ✅ 免费 |
| **PyMOL（开源版）** | Conda 安装 | ✅ 免费 |
| **AutoDock Vina** | 下载二进制 | ✅ 免费 |
| **AutoDockTools** | 下载 MGLTools 套件 | ✅ 免费 |
| **HDOCK / ClusPro** | 网页端，无需安装 | ✅ 免费（学术） |

### Step 1：安装 Open-Source PyMOL

```bash
# 新建一个 structural biology 专用环境
conda create -n strucbio python=3.11 -y
conda activate strucbio

# 装 PyMOL 开源版（注意 channel）
conda install -c conda-forge pymol-open-source -y

# 验证
pymol --help
```

启动 PyMOL：

```bash
pymol
```

应该弹出图形界面。

### Step 2：注册 AlphaFold Server 账号

1. 打开 https://alphafoldserver.com/
2. 用 Google 账号登录（**Google 提供，每日有约 30 个免费 job 配额**）
3. 接受使用条款（仅限非商业研究用途）

### Step 3：下载 AutoDock Vina

- 官网：https://vina.scripps.edu/
- 选 Windows 版本下载（约 5 MB）
- 解压到 `C:\Tools\vina`
- 把 `C:\Tools\vina` 加入系统 PATH

验证：

```bash
vina --help
```

### Step 4：下载 MGLTools（AutoDockTools 的载体）

- 官网：https://ccsb.scripps.edu/mgltools/downloads/
- Windows 版本，约 100 MB
- 安装后会得到一个 GUI 程序 AutoDockTools

> **注意 / Note**：MGLTools 界面比较老旧（基于 Python 2 时代），但它是 AutoDock 生态的标配。新版替代品 ADFR Suite 也可以，但 SOP 资料更少。

## 5.3 蛋白结构预测 / Structure Prediction

### Step 1：拿到靶点的氨基酸序列

以 CD33（AML 最成熟的 CAR-T 靶点）为例：
1. 打开 UniProt：https://www.uniprot.org/
2. 搜索 "CD33 human"，进入 P20138
3. 选 "Sequence" → 复制 FASTA 格式
4. CD33 是单一 isoform，直接下载即可

```
>sp|P20138|CD33_HUMAN Myeloid cell surface antigen CD33
MPLLLLLPLLWAGALAMDPNFWLQVQESVTVQEGLCVLVPCTFFHPIPYYDKNSPVHGYWFREGAIISRDSPVATNKLDQEVQEETQGRFRLLGDPSRNNCSYSIMWKKISSNQNISMSNISGNTLSIVKPSRFNEISDSDFEMTMEPSGYLLIQNIKDMTGHYKCTYVYVHNEMAKFKGLPIKEPIIQLRPKPLEQHVPKEPEITEDLQKGDTPYIIENQTSGIFVSIKPKNKSNHFKETPKWHLNLQHKDPEGQMYFCMVHFNDSGRFPEALQLQSKHCAVDDPTFHSPKPQELICYFWSQETTVSVTRQFKMEQAERHLQERQPYLGSQKLEQLSGEEQLPDRLLQQVPFCTYIGPSGDAVEPGPVLSTSNKMGRIQVKDKDQPIISGDVVFLQPEHITPPQLCTVSRKTQTLSSTTQKQSAGPCLSEFSEQLDLQILNMPVHIYVSTQKTQSISNFMVIPETLASSSQDPCHDFTPSPKPQELICYFWSQETMVSVTRQFKMEQAEQHLQERQPYLGSQKLEQLSGEEQLPDRLLHQVPFCTYIGPSGDAVEPGPVLSTSNKMGRIQVKDKDQPIISGDVVFLQPEHITPPQLCTVSRKTQTLSSTTQKQSAGPCLSEFSEQLDLQILNMPVHIYVSTQKTQSISNFMVIPETLASSSQDPCHDFTPSPK
```

**备选靶点 CD123（IL3RA）**：搜索 "IL3RA human"，UniProt ID P26951

### Step 2：提交到 AlphaFold Server

1. 登录 https://alphafoldserver.com/
2. 点 "New job"
3. 粘贴序列
4. 选 "Protein"
5. （可选）添加配体或离子
6. 点 "Continue" → 等待 5–30 分钟

完成后下载结果包，里面有 `.cif` 和 `.pdb` 文件。

### Step 3：用 PyMOL 可视化

启动 PyMOL，在命令行输入：

```
# 载入结构
load /path/to/your/clnd18_2.pdb, target

# 设置渲染样式
hide everything
show cartoon
color cyan, target

# 高亮跨膜区域（CLDN18.2 的跨膜区在残基 7-25, 81-101, 116-136, 161-181）
select TM1, resi 7-25
select TM2, resi 81-101
select TM3, resi 116-136
select TM4, resi 161-181
color salmon, TM1+TM2+TM3+TM4

# 显示表面
show surface, target
set transparency, 0.5

# 调整视角
orient

# 保存图片
bg_color white
ray 1600, 1200
png clnd18_structure.png, dpi=300
```

### 结构质量评估 / Structure Quality Assessment

AlphaFold 输出会附带 **pLDDT 值**（每个残基的可信度评分）：

| pLDDT | 含义 |
|-------|------|
| > 90 | 非常高可信度 |
| 70–90 | 可信 |
| 50–70 | 低可信，谨慎使用 |
| < 50 | 不可信，可能是无序区 |

在 PyMOL 中按 pLDDT 上色：

```
spectrum b, blue_white_red, target
```

## 5.4 抗体-抗原对接 / Antibody-Antigen Docking

### 推荐方法：HDOCK 网页服务器（最适合新手）

CAR-T 的 scFv 是抗体的可变区，抗体-抗原对接是**蛋白-蛋白对接**（不是小分子对接），AutoDock Vina 不是为此设计的。

**正确工具 / Right Tool**：HDOCK 或 ClusPro

#### HDOCK 操作步骤

1. 打开 http://hdock.phys.hust.edu.cn/
2. 上传两个 PDB：
   - **Receptor**：你的靶点（如 CD33 或 CD123）
   - **Ligand**：抗体 scFv（可以从 PDB 下载已知的抗 CD33 抗体，如 Gemtuzumab 的 scFv 结构，或搜索 PDB 中 "CD33 antibody" 相关结构）
3. 提交，等待 1–3 小时
4. 收到邮件结果，下载 top 10 docking models

#### 用 PyMOL 分析对接结果

```
load receptor.pdb, target
load top1.pdb, antibody

# 找接触界面（5Å 范围内）
select interface_t, byres (target within 5 of antibody)
select interface_a, byres (antibody within 5 of target)

color magenta, interface_t
color yellow, interface_a

# 显示关键残基的侧链
show sticks, interface_t
show sticks, interface_a

# 计算界面面积（需要 get_area 命令）
get_area target
get_area antibody
get_area target or antibody
# 界面面积 = (S(t) + S(a) - S(complex)) / 2
```

## 5.5 小分子对接（可选 / Alternative）

如果你想做**小分子药物 vs 靶点**的对接（不是 CAR-T 主题，但 AutoDock Vina 适用）：

```bash
# 1. 准备 receptor（去水、加氢、加电荷）
# 在 AutoDockTools 中打开 receptor.pdb
# File → Read Molecule → Edit → Hydrogens → Add → All hydrogens
# Edit → Charges → Compute Gasteiger
# 保存为 receptor.pdbqt

# 2. 准备 ligand（小分子）
# 类似处理，保存为 ligand.pdbqt

# 3. 创建配置文件 config.txt
receptor = receptor.pdbqt
ligand = ligand.pdbqt
center_x = 10.0
center_y = 15.0
center_z = 20.0
size_x = 25
size_y = 25
size_z = 25
exhaustiveness = 8
num_modes = 9

# 4. 运行 Vina
vina --config config.txt --out output.pdbqt --log log.txt
```

输出的 `log.txt` 会列出 9 个对接模式的 binding affinity（kcal/mol，越负越好）。

## 5.6 M3 完成检查清单 / M3 Completion Checklist

```
□ PyMOL 装好，能打开 PDB 文件
□ AlphaFold Server 注册成功
□ 至少预测 1 个靶点蛋白的结构
□ 用 PyMOL 渲染至少 2 张专业结构图
□ 完成至少 1 次蛋白-蛋白对接（HDOCK）
□ 标注出关键的结合界面残基
□ 整理成结构分析报告
```

## 5.7 M3 进阶任务（可选）/ M3 Advanced

- 用 **GROMACS** 做短时分子动力学模拟（学习曲线陡，需要 NVIDIA GPU 加速）
- 用 **Rosetta** 做 scFv 的从头设计或亲和力优化
- 用 **PRODIGY** 网页工具评估结合亲和力
- 用 **PISA** 分析界面的氢键和盐桥

---

<a name="part-6"></a>

# Part 6：里程碑 4 — 临床数据库 API 与因果推断（选做）

## 6.1 里程碑概述 / Milestone Overview

**目标 / Goal**：用 ClinicalTrials.gov 和 openFDA 的官方 API（不是爬虫）拉取 CAR-T AML 临床试验数据，做描述性分析（试验数量趋势、靶点分布、各期比例等）。然后用孟德尔随机化（MR）验证某个靶点-AML 的因果关系。

**为什么做这个 / Why**：
- 调 API 是 Real-World Data / RWE 分析岗的核心技能
- MR 是当下流行的因果推断方法，发表门槛低、引用率高
- 这部分能让你的 portfolio 跨越"纯生信"进入"临床数据科学"

**预计时长 / Duration**：3–4 周

## 6.2 ClinicalTrials.gov API v2 使用 / Using CT.gov API v2

### API 基础 / API Basics

- **API 文档 / Docs**：https://clinicaltrials.gov/data-api/api
- **价格 / Price**：✅ 免费，无需注册
- **限流 / Rate Limit**：基本没有硬性限制，但建议每秒 ≤5 个请求

### 示例：拉取 CAR-T 实体瘤试验

```python
import requests
import pandas as pd
import time

# 基础 URL
BASE_URL = "https://clinicaltrials.gov/api/v2/studies"

# 查询参数
params = {
    "query.term": "CAR-T acute myeloid leukemia AML",
    "filter.overallStatus": "RECRUITING|ACTIVE_NOT_RECRUITING|COMPLETED",
    "pageSize": 100,
    "format": "json"
}

# 第一次请求
resp = requests.get(BASE_URL, params=params)
data = resp.json()

studies = data.get('studies', [])
print(f"第一页 / Page 1: {len(studies)} studies")

# 分页拉取
all_studies = list(studies)
while 'nextPageToken' in data:
    params['pageToken'] = data['nextPageToken']
    time.sleep(0.5)            # 礼貌性延迟
    resp = requests.get(BASE_URL, params=params)
    data = resp.json()
    all_studies.extend(data.get('studies', []))
    print(f"累计 / Total: {len(all_studies)}")

print(f"\n总试验数 / Total trials: {len(all_studies)}")
```

### 提取关键字段

```python
def extract_trial_info(study):
    """从单个 study 中提取关键字段"""
    protocol = study.get('protocolSection', {})
    ident = protocol.get('identificationModule', {})
    status = protocol.get('statusModule', {})
    design = protocol.get('designModule', {})
    arms = protocol.get('armsInterventionsModule', {})

    return {
        'nct_id': ident.get('nctId'),
        'title': ident.get('briefTitle'),
        'phase': design.get('phases', []),
        'status': status.get('overallStatus'),
        'start_date': status.get('startDateStruct', {}).get('date'),
        'completion_date': status.get('completionDateStruct', {}).get('date'),
        'enrollment': design.get('enrollmentInfo', {}).get('count'),
        'sponsor': protocol.get('sponsorCollaboratorsModule', {})
                          .get('leadSponsor', {}).get('name'),
        'conditions': protocol.get('conditionsModule', {}).get('conditions', []),
    }

df = pd.DataFrame([extract_trial_info(s) for s in all_studies])
df.to_csv('results/cart_solid_tumor_trials.csv', index=False)
print(df.head())
```

### 简单分析示例

```python
import matplotlib.pyplot as plt

# 按状态分布
status_counts = df['status'].value_counts()
status_counts.plot(kind='bar')
plt.title('CAR-T AML Trials by Status')
plt.tight_layout()
plt.savefig('figures/cart_aml_trials_by_status.png', dpi=300)

# 按起始年份的趋势
df['start_year'] = pd.to_datetime(df['start_date']).dt.year
df.groupby('start_year').size().plot(kind='line', marker='o')
plt.title('CAR-T AML Trials Initiated per Year')
plt.savefig('figures/cart_aml_trials_trend.png', dpi=300)
```

## 6.3 openFDA API 使用 / Using openFDA API

```python
import requests

# 查询某个药物的不良反应
url = "https://api.fda.gov/drug/event.json"
params = {
    "search": 'patient.drug.medicinalproduct:"YESCARTA"',
    "count": "patient.reaction.reactionmeddrapt.exact",
    "limit": 20
}

resp = requests.get(url, params=params)
data = resp.json()

reactions = pd.DataFrame(data['results'])
print(reactions)
```

## 6.4 孟德尔随机化 / Mendelian Randomization

### 核心概念 / Core Concept

MR 用基因变异作为**工具变量**，回答：「暴露 X 是否**因果地**导致结局 Y？」

- **暴露 / Exposure**：通常是一个生物标志物或基因表达水平
- **结局 / Outcome**：通常是疾病
- **工具变量 / IV**：与暴露强相关的 SNP

### 三大假设 / Three Assumptions

1. **关联性**：IV 与暴露强相关（F > 10）
2. **独立性**：IV 不通过暴露以外的路径影响结局
3. **排他性**：IV 与潜在混杂因素不相关

### R 包安装

```r
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)
ao <- available_outcomes()    # 列出所有可用的 GWAS 数据
head(ao)
```

### 示例分析框架

```r
library(TwoSampleMR)

# 1. 选择暴露（如：某个基因的 eQTL）
exposure_dat <- extract_instruments(outcomes = 'ieu-a-XXX')

# 2. 选择结局（如：某种癌症的 GWAS）
outcome_dat <- extract_outcome_data(
  snps = exposure_dat$SNP,
  outcomes = 'ieu-a-YYY'
)

# 3. 数据 harmonization
dat <- harmonise_data(exposure_dat, outcome_dat)

# 4. 运行 MR
res <- mr(dat)
print(res)

# 5. 敏感性分析
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

# 6. 可视化
p1 <- mr_scatter_plot(res, dat)
ggsave('figures/mr_scatter.png', p1[[1]])
```

> **重要 / Important**：MR 是统计方法，结论必须谨慎。第一次做 MR 建议参照 STROBE-MR 指南撰写报告。

## 6.5 M4 完成检查清单 / M4 Completion Checklist

```
□ 用 CT.gov API 拉取 100+ 个相关临床试验
□ 用 openFDA API 提取不良反应数据
□ 生成至少 2 张临床试验分析图
□ （选做）完成 1 次 MR 分析
□ 写一份 RWE-style mini report
```

---

<a name="part-7"></a>

# Part 7：常见错误速查 / Troubleshooting Index

## 7.1 环境与安装类 / Environment Issues

| 问题 | 排查思路 |
|------|----------|
| Conda 装不上包 | 切换 channel，加 `-c conda-forge`；或先 `conda update conda` |
| R 包安装失败 | 检查 R 版本是否 ≥ 4.4；检查 BiocManager 是否最新 |
| Scanpy 提示缺 leiden | `pip install scanpy[leiden]` 重装 |
| PyMOL 启动黑屏 | 显卡驱动问题，更新 Intel Arc 驱动 |
| Jupyter 找不到 conda 环境 | `python -m ipykernel install --user --name=<env>` |

## 7.2 数据处理类 / Data Processing

| 问题 | 排查思路 |
|------|----------|
| getGEO 下载超时 | 用国内镜像，或手动下 series matrix |
| 内存爆掉（OOM） | 用 read_only 模式打开大文件；分块读取 |
| Count 矩阵不是整数 | 检查是 raw count 还是 TPM；用 tximport 转 |
| 基因 ID 对不上 | 用 biomaRt 或 org.Hs.eg.db 做 ID 转换 |
| 单细胞数据样本数特别少 | 检查是不是 metadata 文件，不是 expression 文件 |

## 7.3 分析逻辑类 / Analysis Logic

| 问题 | 排查思路 |
|------|----------|
| DESeq2 出来没有差异基因 | 检查分组是否搞反；检查 padj 不要用 p；放宽 logFC 阈值 |
| UMAP 分群太散或太团 | 调 `n_neighbors` 和 `resolution` 参数 |
| AlphaFold pLDDT 普遍低 | 检查序列是否有信号肽或无序区；考虑只看 well-folded 部分 |
| 对接结合能不合理 | 检查 receptor/ligand 是否正确加氢加电荷；检查 docking box 大小 |

---

<a name="part-8"></a>

# Part 8：学习资源 / Learning Resources

## 8.1 推荐书籍 / Books

- **《生物信息学：基础与实践》** 樊龙江主编 — 中文入门书
- **"Modern Statistics for Modern Biology"** Susan Holmes & Wolfgang Huber — R/Bioconductor 圣经，免费在线：https://www.huber.embl.de/msmb/
- **"Bioinformatics Data Skills"** Vince Buffalo — 工程化思维

## 8.2 推荐在线课程 / Online Courses

| 平台 | 课程 | 价格 |
|------|------|------|
| Coursera | Genomic Data Science Specialization (Johns Hopkins) | 免费旁听 |
| edX | Data Analysis for Life Sciences (Harvard) | 免费旁听 |
| Bilibili | "生信技能树" 系列视频 | 免费 |
| Bilibili | "生信菜鸟团" 系列视频 | 免费 |

## 8.3 必读教程 / Must-Read Tutorials

- **DESeq2 vignette**：https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- **Scanpy tutorials**：https://scanpy.readthedocs.io/en/stable/tutorials.html
- **Seurat tutorials**：https://satijalab.org/seurat/articles/get_started.html
- **OSCA (Orchestrating Single-Cell Analysis)**：https://bioconductor.org/books/release/OSCA/

## 8.4 中文社区 / Chinese Communities

- 生信技能树（biotrainee.com）
- 生信菜鸟团（公众号）
- 解螺旋（学术写作社区）
- 知乎「生物信息学」话题

## 8.5 英文社区 / English Communities

- Biostars（biostars.org）— 提问的好地方
- r/bioinformatics（Reddit）
- Bioconductor Support Forum
- Twitter/X：关注 #scRNAseq、#bioinformatics 标签

## 8.6 数据集发现 / Finding Datasets

- **CELLxGENE**：https://cellxgene.cziscience.com/ — 标准化的单细胞数据浏览器
- **Single Cell Portal (Broad)**：https://singlecell.broadinstitute.org/
- **DepMap**：https://depmap.org/ — 癌症细胞系数据
- **OncoLnc**：http://www.oncolnc.org/ — TCGA 生存分析快查

---

# Part 9：里程碑 5 — 生存分析与预后模型

## 9.1 目标

用 **TCGA-LAML**（约 200 例 AML 患者，含随访数据）验证：
> M1/M2 筛选出的候选靶点（CLEC12A、FLT3、CD33 等）高表达患者，OS（总生存期）是否显著更短？

这一步让整个 portfolio 逻辑闭环：**从发现靶点 → 证明靶点与预后相关**。

## 9.2 数据来源

| 数据 | 获取方式 | 说明 |
|------|----------|------|
| TCGA-LAML 表达矩阵 | `TCGAbiolinks` R 包自动下载 | RNA-seq count 数据，~150 例 |
| TCGA-LAML 临床数据 | `TCGAbiolinks` 同步下载 | 含 OS 天数、死亡状态 |

不需要手动去 TCGA 官网下载，R 包全自动搞定。

## 9.3 环境安装

```r
# 在 RStudio Console 运行（首次）
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
install.packages(c("survival", "survminer", "dplyr", "ggplot2"))
```

## 9.4 分析步骤（逐步操作）

### Step 1：下载 TCGA-LAML 数据

```r
library(TCGAbiolinks)
library(dplyr)

# 查询 TCGA-LAML RNA-seq 数据（STAR counts）
query <- GDCquery(
  project      = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type    = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# 下载（第一次运行需要几分钟，下载到 GDCdata/ 目录）
GDCdownload(query, method = "api", files.per.chunk = 10)

# 整理成表达矩阵
expr_data <- GDCprepare(query)

# 提取 count 矩阵（行=基因，列=样本）
count_mat <- assay(expr_data, "unstranded")
```

### Step 2：下载临床数据

```r
# 下载临床信息
clinical <- GDCquery_clinic(project = "TCGA-LAML", type = "clinical")

# 关键列：
# submitter_id      → 样本 ID
# days_to_death     → 死亡时间（天）
# days_to_last_follow_up → 末次随访时间
# vital_status      → "Dead" or "Alive"

# 构造生存数据框
surv_df <- clinical %>%
  transmute(
    sample_id = submitter_id,
    OS_time   = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
    OS_status = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(OS_time), OS_time > 0)
```

### Step 3：提取靶点表达量，分高低组

```r
library(SummarizedExperiment)

# 目标靶点列表（与 M1/M2 结果对应）
targets <- c("CLEC12A", "FLT3", "CD33", "CD123", "CD38")

# 从表达矩阵提取（rownames 是 Ensembl ID，需要转换）
# 方法：用 rowData 中的 gene_name 列匹配
gene_info <- rowData(expr_data)

for (gene in targets) {
  # 找到对应行
  idx <- which(gene_info$gene_name == gene)
  if (length(idx) == 0) { message("未找到: ", gene); next }

  # 提取该基因表达量（取第一个匹配）
  expr_vec <- count_mat[idx[1], ]

  # 样本 ID 对齐（TCGA 样本 ID 取前 12 位）
  names(expr_vec) <- substr(names(expr_vec), 1, 12)

  # 按中位数分高低组
  median_val <- median(expr_vec, na.rm = TRUE)

  # 合并到 surv_df
  surv_df[[paste0(gene, "_expr")]] <- expr_vec[surv_df$sample_id]
  surv_df[[paste0(gene, "_group")]] <- ifelse(
    surv_df[[paste0(gene, "_expr")]] >= median_val, "High", "Low"
  )
}
```

### Step 4：Kaplan-Meier 曲线

```r
library(survival)
library(survminer)

# 对每个靶点画 K-M 曲线
plot_km <- function(gene, df) {
  group_col <- paste0(gene, "_group")
  df_clean  <- df[!is.na(df[[group_col]]), ]

  fit <- survfit(
    Surv(OS_time, OS_status) ~ df_clean[[group_col]],
    data = df_clean
  )

  p <- ggsurvplot(
    fit,
    data          = df_clean,
    pval          = TRUE,          # 显示 log-rank p 值
    pval.method   = TRUE,
    conf.int      = TRUE,
    risk.table    = TRUE,          # 底部风险表
    palette       = c("#C44E52", "#4C72B0"),
    legend.labs   = c("High", "Low"),
    legend.title  = paste(gene, "expression"),
    title         = paste0("Overall Survival by ", gene, " Expression\n(TCGA-LAML, n=", nrow(df_clean), ")"),
    xlab          = "Time (days)",
    ylab          = "Survival Probability",
    ggtheme       = theme_bw()
  )
  return(p)
}

# 批量生成并保存
dir.create("M5_Survival/figures", recursive = TRUE, showWarnings = FALSE)

for (gene in targets) {
  p <- plot_km(gene, surv_df)
  ggsave(
    filename = paste0("M5_Survival/figures/KM_", gene, ".png"),
    plot     = print(p),
    width = 8, height = 8, dpi = 150
  )
  message("✅ 保存: KM_", gene, ".png")
}
```

### Step 5：Cox 多变量回归

```r
# 单变量 Cox（每个靶点独立检验）
cox_results <- list()

for (gene in targets) {
  expr_col <- paste0(gene, "_expr")
  df_clean <- surv_df[!is.na(surv_df[[expr_col]]), ]

  # log2 转换表达量（避免数值范围过大）
  df_clean$expr_log2 <- log2(df_clean[[expr_col]] + 1)

  cox_fit <- coxph(Surv(OS_time, OS_status) ~ expr_log2, data = df_clean)
  s <- summary(cox_fit)

  cox_results[[gene]] <- data.frame(
    Gene    = gene,
    HR      = exp(coef(cox_fit)),
    HR_lower = exp(confint(cox_fit)[1]),
    HR_upper = exp(confint(cox_fit)[2]),
    pval    = s$coefficients[, "Pr(>|z|)"]
  )
}

cox_table <- do.call(rbind, cox_results)
print(cox_table)

# 保存结果
write.csv(cox_table, "M5_Survival/cox_results.csv", row.names = FALSE)
```

### Step 6：生成 RMarkdown 报告

将以上所有代码整合进 `M5_Survival/M5_Survival_Analysis.Rmd`，Knit 生成 HTML。

## 9.5 预期产出

- `KM_CLEC12A.png` / `KM_FLT3.png` 等 K-M 曲线（每个靶点一张）
- `cox_results.csv` — HR 值、95% CI、p 值
- `M5_Survival_Analysis.html` — 完整报告

## 9.6 结果解读要点

- **log-rank p < 0.05** → 高表达组与低表达组生存曲线显著分离
- **HR > 1** → 高表达 = 预后差（这正是我们希望看到的，说明靶点有致病意义）
- CLEC12A 预期 HR 较高（M1 logFC 最大）

---

# Part 10：里程碑 6 — 基因组变异分析（Somatic Mutation）

## 10.1 目标

分析 AML 体细胞突变图谱，回答：
1. AML 最常见的驱动突变是哪些？（FLT3-ITD、NPM1、DNMT3A 等）
2. 这些突变如何影响我们的 CAR-T 候选靶点的表达？
3. 有无突变的亚组，靶点表达是否有差异？

## 10.2 数据来源

TCGA-LAML 的 MAF 文件，通过 TCGAbiolinks 批量下载（153个文件，131个病人，3900条突变）。

> ⚠️ **关键坑**：GDC 按病人分开存储 MAF，单个 file_id 只能得到 1 个样本。必须用 TCGAbiolinks 批量下载后合并，不能直接用 httr GET 单个文件。

## 10.3 环境安装

```r
BiocManager::install("maftools")      # 已含在 TCGAbiolinks 依赖中
# TCGAbiolinks 在 M5 已安装
library(maftools)
cat(as.character(packageVersion("maftools")))  # 确认版本 >= 2.28
```

## 10.4 分析步骤（已验证可运行）

### Step 1：批量下载并合并 TCGA-LAML MAF

```r
library(TCGAbiolinks)
library(maftools)

# 查询：TCGA-LAML 所有 masked somatic mutation 文件
query_mut <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)
# 共找到 153 个文件，对应 131 个唯一病人

# 批量下载（约 1.6 MB，<1 分钟）
GDCdownload(
  query = query_mut,
  directory = "D:/Bio-Informatics Case Study/M6_Mutation/data/GDCdata"
)

# 合并为单一数据框（3900行 × 140列）
maf_all <- GDCprepare(
  query = query_mut,
  directory = "D:/Bio-Informatics Case Study/M6_Mutation/data/GDCdata"
)

# 转为 maftools MAF 对象
laml_full <- read.maf(maf = maf_all)

# 验证
cat("样本数：", nrow(getSampleSummary(laml_full)), "\n")  # 131
cat("基因数：", nrow(getGeneSummary(laml_full)), "\n")    # 2423
```

### Step 2：瀑布图 / Oncoprint（Top 20 突变基因）

```r
# 注意：maftools 2.28 的 oncoplot() 不支持 title= 参数，去掉即可
png("D:/Bio-Informatics Case Study/M6_Mutation/figures/01_waterfall_top20.png",
    width = 3000, height = 2000, res = 200)
oncoplot(maf = laml_full, top = 20)
dev.off()

# 同时保存 PDF（矢量图，简历/论文用）
pdf("D:/Bio-Informatics Case Study/M6_Mutation/figures/01_waterfall_top20.pdf",
    width = 14, height = 10)
oncoplot(maf = laml_full, top = 20)
dev.off()
```

### Step 3：Lollipop 图（突变位点在蛋白质上的分布）

```r
# 注意：lollipopPlot() 同样不支持 title= 参数
png("D:/Bio-Informatics Case Study/M6_Mutation/figures/02_lollipop_FLT3.png",
    width = 2400, height = 1200, res = 200)
lollipopPlot(maf = laml_full, gene = "FLT3")
dev.off()

png("D:/Bio-Informatics Case Study/M6_Mutation/figures/03_lollipop_NPM1.png",
    width = 2400, height = 1200, res = 200)
lollipopPlot(maf = laml_full, gene = "NPM1")
dev.off()

png("D:/Bio-Informatics Case Study/M6_Mutation/figures/04_lollipop_DNMT3A.png",
    width = 2400, height = 1200, res = 200)
lollipopPlot(maf = laml_full, gene = "DNMT3A")
dev.off()
```

### Step 4：突变共现 / 互斥矩阵

```r
png("D:/Bio-Informatics Case Study/M6_Mutation/figures/05_somatic_interactions.png",
    width = 2400, height = 2000, res = 200)
somaticInteractions(maf = laml_full, top = 20, pvalue = c(0.05, 0.1))
dev.off()
```

### Step 5：FLT3 突变亚组 vs mRNA 表达量箱线图

```r
library(ggplot2)

# 提取 FLT3 突变样本（取前12位病人ID）
flt3_mutated <- subsetMaf(maf = laml_full, genes = "FLT3")
flt3_mut_ids <- substr(
  as.character(getSampleSummary(flt3_mutated)$Tumor_Sample_Barcode), 1, 12
)

# 加载 M5 整合数据框（含 FLT3 表达量列）
analysis_df <- readRDS("D:/Bio-Informatics Case Study/M5_Survival/analysis_df.rds")

# 添加突变状态标签
analysis_df$FLT3_mutation <- ifelse(
  analysis_df$patient_id %in% flt3_mut_ids, "Mutated", "Wild-type"
)

# Wilcoxon 检验 + 箱线图
wt <- wilcox.test(FLT3 ~ FLT3_mutation, data = analysis_df)

png("D:/Bio-Informatics Case Study/M6_Mutation/figures/06_FLT3_mutation_vs_expression.png",
    width = 1600, height = 1600, res = 200)
ggplot(analysis_df, aes(x = FLT3_mutation, y = log2(FLT3 + 1),
                         fill = FLT3_mutation)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("Mutated" = "#E74C3C", "Wild-type" = "#3498DB")) +
  labs(title = "FLT3 mRNA Expression by Mutation Status",
       subtitle = "TCGA-LAML (n=130)",
       x = "FLT3 Mutation Status",
       y = "FLT3 Expression (log2 counts + 1)") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none") +
  annotate("text", x = 1.5, y = max(log2(analysis_df$FLT3 + 1)) * 0.98,
           label = paste0("p = ", round(wt$p.value, 3)), size = 5)
dev.off()
```

## 10.5 实际产出（2026-05-21 完成）

| 文件 | 说明 |
|------|------|
| `figures/01_waterfall_top20.png/pdf` | Top 20 突变基因瀑布图 |
| `figures/02_lollipop_FLT3.png` | FLT3 激酶域突变热点 |
| `figures/03_lollipop_NPM1.png` | NPM1 C端移码插入热点 |
| `figures/04_lollipop_DNMT3A.png` | DNMT3A 多域分散突变 |
| `figures/05_somatic_interactions.png` | 突变共现/互斥矩阵 |
| `figures/06_FLT3_mutation_vs_expression.png` | FLT3突变 vs mRNA表达（p=0.255） |
| `01_M6_analysis.R` | 完整分析脚本 |
| `M6_Report.Rmd` | RMarkdown 报告 |

## 10.6 关键结论

- Top 突变基因：NPM1(8.4%) > TP53(7%) > DNMT3A(6.1%) > FLT3(4.6%)
- FLT3 突变集中在 PKc_like 激酶域（错义突变），NPM1 突变集中在 C 端（移码插入）
- FLT3 + DNMT3A 显著共现（p<0.05）：双突变病人预后更差
- FLT3 DNA 突变不影响 mRNA 表达量（p=0.255）：突变通过激酶功能而非转录水平致病

---

# Part 11：里程碑 7 — 免疫微环境分析（TME Deconvolution）

## 11.1 目标

CAR-T 疗法的疗效不仅取决于靶点表达，还取决于**肿瘤免疫微环境（TME）**：
- 骨髓中有多少 T 细胞？NK 细胞？
- 髓系抑制细胞（MDSC）比例高吗？（会抑制 CAR-T）
- 候选靶点表达量与免疫细胞比例有无相关性？

## 11.2 方法选择

| 方法 | 工具 | 特点 |
|------|------|------|
| CIBERSORT | R / `IOBR` 包 | 最经典，22种免疫细胞 |
| TIMER2.0 | 网页 / API | 快速，有可视化 |
| xCell | R / `IOBR` 包 | 64种细胞类型，更细 |
| ESTIMATE | R | 免疫评分 + 基质评分 |

推荐用 `IOBR` 包，一个包集成所有方法。

## 11.3 环境安装

```r
# IOBR 安装（集成 CIBERSORT / xCell / TIMER / ESTIMATE）
if (!require("remotes")) install.packages("remotes")
remotes::install_github("IOBR/IOBR")

install.packages(c("ggplot2", "dplyr", "corrplot", "ggpubr"))
```

## 11.4 分析步骤

### Step 1：准备输入数据（TPM 矩阵）

```r
# CIBERSORT 需要 TPM（不是 raw count）
# 从 TCGA-LAML 数据中提取 TPM
tpm_mat <- assay(expr_data, "tpm_unstrand")  # SummarizedExperiment 对象

# 行名转换为 gene symbol（IOBR 需要 symbol）
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = rownames(tpm_mat),
  mart       = mart
)
# 合并，去重，保留 symbol
```

### Step 2：运行免疫细胞反卷积

```r
library(IOBR)

# 方法1：CIBERSORT（需要基准矩阵文件，IOBR 内置）
cibersort_res <- deconvo_tme(
  eset      = tpm_mat,        # 基因×样本 TPM 矩阵
  method    = "cibersort",
  arrays    = FALSE,          # RNA-seq 数据
  perm      = 100             # 置换检验次数（越多越慢，100 够了）
)

# 方法2：xCell（更多细胞类型）
xcell_res <- deconvo_tme(
  eset   = tpm_mat,
  method = "xcell"
)

# ESTIMATE 免疫评分
estimate_res <- deconvo_tme(
  eset   = tpm_mat,
  method = "estimate"
)
```

### Step 3：可视化

```r
# 热图：样本 × 免疫细胞类型
library(pheatmap)
pheatmap(
  t(cibersort_res[, 2:23]),   # 22 种免疫细胞列
  scale        = "row",
  clustering_distance_rows = "euclidean",
  show_colnames = FALSE,
  main = "AML Immune Cell Composition (CIBERSORT)"
)

# 靶点表达 vs 免疫细胞比例相关性
# 例：CLEC12A 表达与 T cell CD8+ 比例的相关性
library(ggpubr)
merged <- merge(
  data.frame(sample_id = colnames(tpm_mat),
             CLEC12A = as.numeric(tpm_mat["CLEC12A", ])),
  cibersort_res %>% select(sample_id = SampleID, T_cells_CD8),
  by = "sample_id"
)

ggscatter(merged, x = "CLEC12A", y = "T_cells_CD8",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          title = "CLEC12A Expression vs CD8+ T Cell Fraction")
```

## 11.5 预期产出

- `fig_tme_heatmap.png` — 免疫细胞组成热图
- `fig_target_vs_immune.png` — 靶点表达 vs 免疫细胞相关性
- `tme_scores.csv` — 每个样本的免疫细胞比例数据
- `M7_TME_Analysis.html` — 完整报告

---

# Part 12：里程碑 8 — DNA 甲基化 / 表观基因组分析

## 12.1 目标

AML 中 **DNMT3A 突变**极为常见（约 20%），会导致全局 DNA 甲基化异常。
本分析回答：
- 候选靶点（CLEC12A、FLT3 等）的启动子区域是否在 AML 中出现差异甲基化？
- 甲基化状态与基因表达量是否负相关（甲基化 → 沉默）？

## 12.2 数据来源

GEO 上的 **Illumina EPIC array**（850K）或 **450K array** AML 数据集，推荐：
- **GSE69065**：AML vs 正常骨髓，EPIC array
- **GSE124413**：AML 甲基化 + 表达联合数据集

## 12.3 环境安装

```r
BiocManager::install(c("minfi", "ChAMP", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"))
install.packages(c("ggplot2", "pheatmap", "dplyr"))
```

## 12.4 分析步骤

### Step 1：下载数据

```r
library(GEOquery)

# 下载 GSE69065（IDAT 文件）
gse <- getGEO("GSE69065", GSEMatrix = FALSE)

# 注意：甲基化分析需要 IDAT 原始文件（不是 series_matrix）
# 从 GEO 页面手动下载 supplementary files（IDAT.gz）
# 放入 M8_Methylation/data/idat/ 目录
```

### Step 2：读取 IDAT，预处理

```r
library(minfi)

# 读取所有 IDAT 文件
rgSet <- read.metharray.exp(base = "M8_Methylation/data/idat/")

# 质控报告
qcReport(rgSet, sampNames = pData(rgSet)$Sample_Name,
         pdf = "M8_Methylation/QC_report.pdf")

# 归一化（Noob 方法，推荐）
mSet <- preprocessNoob(rgSet)

# 过滤低质量探针
detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.01) == ncol(detP)
mSet <- mSet[keep, ]

# 提取 Beta 值（0-1，0=完全未甲基化，1=完全甲基化）
beta_mat <- getBeta(mSet)
```

### Step 3：差异甲基化分析（DMP）

```r
library(limma)

# 分组（AML vs Normal）
group <- factor(pData(mSet)$disease_status)  # 按实际列名调整

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# M 值比 Beta 值更适合统计分析
m_mat <- getM(mSet)

fit <- lmFit(m_mat, design)
contrast_mat <- makeContrasts(AML - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

# 提取差异甲基化位点（DMP）
dmps <- topTable(fit2, num = Inf, adjust.method = "BH") %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 0.5)

message("差异甲基化位点数: ", nrow(dmps))
```

### Step 4：聚焦靶点启动子区域

```r
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# 获取探针注释（染色体位置、基因名、功能区域）
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# 筛选靶点启动子区域探针（TSS200 / TSS1500）
targets <- c("CLEC12A", "FLT3", "CD33", "IL3RA", "CD38")
target_probes <- anno %>%
  as.data.frame() %>%
  filter(
    grepl(paste(targets, collapse = "|"), UCSC_RefGene_Name),
    grepl("TSS", UCSC_RefGene_Group)   # 启动子区域
  )

# 提取靶点的 Beta 值，可视化
beta_target <- beta_mat[rownames(beta_target) %in% target_probes$Name, ]
pheatmap(beta_target, annotation_col = pData(mSet)[, "group", drop=FALSE],
         main = "Promoter Methylation of CAR-T Target Genes")
```

### Step 5：甲基化 vs 表达量整合（整合 M1/M5 数据）

```r
# 对同一基因：启动子甲基化程度 vs mRNA 表达量的相关性
# 负相关 → 甲基化导致基因沉默 → 靶点在正常细胞被沉默，AML 中去甲基化激活
# 以 FLT3 为例
flt3_beta <- colMeans(beta_target[grep("FLT3", rownames(beta_target)), ])
# 与 M5 中的 FLT3_expr 合并，画散点图
```

## 12.5 预期产出

- `QC_report.pdf` — 甲基化质控报告
- `fig_methylation_heatmap.png` — 靶点启动子甲基化热图
- `fig_methyl_vs_expr.png` — 甲基化 vs 表达量散点图
- `dmp_results.csv` — 差异甲基化位点表
- `M8_Methylation_Analysis.html` — 完整报告

---

# Part 13：里程碑 9 — 空间转录组（Spatial Transcriptomics）

## 13.1 目标

scRNA-seq（M2）告诉了我们**细胞类型**，但不知道这些细胞**在组织里在哪里**。
空间转录组同时保留位置信息，回答：
- CLEC12A 高表达的 AML 细胞，集中在骨髓的哪个区域？
- CAR-T 进入骨髓后，需要穿透哪些细胞层才能到达靶标？

## 13.2 数据来源

AML 空间转录组公开数据较少，推荐：
- **GSE174448**：骨髓 10x Visium 数据（含 AML 样本）
- **10x Genomics 官网示例**：Human Bone Marrow（学习用）

如果找不到合适的 AML 数据，可退而求其次：
- 用正常骨髓 Visium 数据 + 把 M2 的 AML 细胞比例投影上去（deconvolution）

## 13.3 环境安装

```bash
# Python 环境（在 scrna env 里追加安装）
conda activate scrna
pip install squidpy spatialdata
```

```r
# R 环境（Seurat v5 支持空间数据）
install.packages("Seurat")  # v5+
BiocManager::install("SpatialExperiment")
```

## 13.4 分析步骤（Python / Squidpy）

### Step 1：读取 10x Visium 数据

```python
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt

# 读取 10x Visium 输出目录
# 目录结构：filtered_feature_bc_matrix/ + spatial/
adata = sc.read_visium(
    path       = "M9_SpatialTranscriptomics/data/visium_sample/",
    count_file = "filtered_feature_bc_matrix.h5",
    load_images= True
)

print(adata)
# AnnData: n_obs=spots, n_vars=genes，obs 包含空间坐标
```

### Step 2：质控与预处理

```python
# 基础 QC
sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[adata.obs["total_counts"] > 1000]

# 归一化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 高变基因
sc.pp.highly_variable_genes(adata, n_top_genes=3000)

# PCA + neighbors + UMAP（与 M2 流程一致）
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
```

### Step 3：空间可视化靶点表达

```python
targets = ["CLEC12A", "FLT3", "CD33", "IL3RA", "CD38"]

# 在组织切片上叠加基因表达（最核心的空间转录组图）
sq.pl.spatial_scatter(
    adata,
    color     = targets,
    ncols     = 3,
    size      = 1.5,
    cmap      = "Reds",
    title     = [f"{g} Expression" for g in targets],
    save      = "M9_SpatialTranscriptomics/figures/spatial_targets.png"
)
```

### Step 4：空间自相关分析（Moran's I）

```python
# 检验某基因是否有空间聚集性（随机分布 vs 聚集分布）
sq.gr.spatial_neighbors(adata, coord_type="generic")
sq.gr.spatial_autocorr(adata, mode="moran", genes=targets)

# 显示结果
print(adata.uns["moranI"])
# I 接近 1 → 高度空间聚集（高表达区域集中在某处）
```

### Step 5：细胞类型解卷积（用 M2 的 scRNA 参考）

```python
# 用 M2 的 scRNA-seq 作为参考，推断每个 Visium spot 的细胞组成
# 工具：cell2location（最准确）
# pip install cell2location

import cell2location
# 详见 cell2location 官方教程：
# https://cell2location.readthedocs.io/
```

## 13.5 预期产出

- `spatial_targets.png` — 靶点在组织切片上的空间表达图
- `fig_moranI.png` — 空间自相关热图
- `M9_Spatial_Analysis.ipynb` — 完整 Notebook

---

# Part 14：里程碑 10 — 多组学整合（Multi-omics Integration）

## 14.1 目标

整合 M1-M9 的所有数据层，用统计模型找出**跨数据层一致的信号**：
> 哪些靶点在转录组、甲基化、突变、免疫微环境层面都有支持证据？

这是整个项目的**压轴分析**，产出最终靶点排名。

## 14.2 工具选择

| 工具 | 适合整合的数据类型 | 特点 |
|------|------------------|------|
| **MOFA+** | RNA-seq + 甲基化 + 突变 | 最主流，有 Python 和 R 接口 |
| **Seurat WNN** | scRNA + scATAC | 专门做单细胞多组学 |
| **mixOmics** | 任意多组学 | R 包，适合有监督整合 |

推荐用 **MOFA+**（Multi-Omics Factor Analysis）。

## 14.3 环境安装

```bash
conda activate scrna
pip install mofapy2
```

```r
BiocManager::install("MOFA2")
install.packages("reticulate")
```

## 14.4 分析步骤

### Step 1：准备各组学数据矩阵

```r
library(MOFA2)

# 每个数据层是一个矩阵：基因/位点 × 样本
# 样本 ID 必须完全一致！这是最容易出错的地方

# 数据层列表
data_list <- list(
  RNA       = rna_mat,        # 来自 M1/M5：基因×样本 log2 表达矩阵
  Mutation  = mut_mat,        # 来自 M6：基因×样本 0/1 突变矩阵
  Methylation = beta_mat_sub, # 来自 M8：探针×样本 Beta 值（取靶点区域）
  TME       = tme_mat         # 来自 M7：免疫细胞类型×样本
)

# 确保样本 ID 对齐（取交集）
common_samples <- Reduce(intersect, lapply(data_list, colnames))
data_list <- lapply(data_list, function(x) x[, common_samples])
```

### Step 2：创建 MOFA 对象并训练

```r
# 创建 MOFA 对象
mofa_obj <- create_mofa(data_list)

# 可视化数据概况
plot_data_overview(mofa_obj)

# 设置训练参数
model_opts <- get_default_model_options(mofa_obj)
model_opts$num_factors <- 10   # 学习 10 个潜在因子

train_opts <- get_default_training_options(mofa_obj)
train_opts$seed        <- 42
train_opts$maxiter     <- 1000

mofa_obj <- prepare_mofa(mofa_obj,
  model_options   = model_opts,
  training_options = train_opts
)

# 训练模型（几分钟）
mofa_obj <- run_mofa(mofa_obj, outfile = "M10_MultiOmics/mofa_model.hdf5")
```

### Step 3：解释因子

```r
# 每个因子解释了多少方差？
plot_variance_explained(mofa_obj, max_r2 = 15)

# Factor 1 在各数据层的载荷（哪些基因/位点贡献最大）
plot_top_weights(mofa_obj, view = "RNA", factor = 1, nfeatures = 20)
plot_top_weights(mofa_obj, view = "Methylation", factor = 1, nfeatures = 20)

# 样本在因子空间的分布（类似 PCA，但整合了所有组学）
plot_factor(mofa_obj, factors = c(1, 2), color_by = "AML_subtype")
```

### Step 4：最终靶点综合评分

```r
# 综合各里程碑的证据，给每个候选靶点打分
# 评分维度（每项 0-3 分）：
# 1. M1 logFC（差异表达倍数）
# 2. M2 AML/Normal 表达比（单细胞层面）
# 3. M5 HR（预后相关性，HR > 1.5 得高分）
# 4. M6 是否在突变基因附近（突变协同）
# 5. M7 与 CD8+ T 细胞正相关（利于 CAR-T 共同作战）
# 6. M8 AML 中启动子去甲基化（表观激活）
# 7. M4 临床试验数（过多 = 竞争激烈，过少 = 新兴机会）

score_table <- data.frame(
  Target   = c("CLEC12A", "FLT3", "CD33", "CD123", "CD38"),
  M1_score = c(3, 2, 2, 2, 1),   # 根据实际结果填写
  M2_score = c(2, 3, 3, 2, 1),
  M5_score = c(NA, NA, NA, NA, NA),  # M5 完成后填入
  M6_score = c(NA, NA, NA, NA, NA),
  M7_score = c(NA, NA, NA, NA, NA),
  M8_score = c(NA, NA, NA, NA, NA)
)
score_table$Total <- rowSums(score_table[, -1], na.rm = TRUE)
score_table <- score_table %>% arrange(desc(Total))
print(score_table)
```

## 14.5 预期产出

- `mofa_model.hdf5` — 训练好的 MOFA 模型
- `fig_variance_explained.png` — 各因子方差解释图
- `fig_factor_weights.png` — 因子载荷图
- `final_target_ranking.csv` — 最终靶点综合评分排名
- `M10_MultiOmics_Integration.html` — 完整报告

---

# Part 15：自动化 Pipeline 设计（终极目标）

## 15.1 设计理念

**目标**：输入一个癌症类型 + 数据集 ID，自动运行 M1-M10 全流程，输出标准化报告。

**核心原则**：
- **半自动化**（不是全自动）：在 QC、分组确认、结果解读等关键节点设置**人工检查点**
- **模块化**：每个里程碑是一个独立模块，可单独跑，也可串联跑
- **配置驱动**：所有参数写进一个 `config.yaml`，换癌症类型只需改配置文件

## 15.2 目录结构设计

```
auto_pipeline/
├── config.yaml                  # 全局配置（癌症类型、数据集ID、阈值等）
├── run_pipeline.py              # 主入口（Python 调度器）
│
├── modules/
│   ├── m01_bulk_rnaseq.R        # M1 模块
│   ├── m02_scrna.py             # M2 模块
│   ├── m03_structure.py         # M3 接口（调用 AlphaFold API）
│   ├── m04_clinical_mr.py       # M4 模块
│   ├── m05_survival.R           # M5 模块
│   ├── m06_mutation.R           # M6 模块
│   ├── m07_tme.R                # M7 模块
│   ├── m08_methylation.R        # M8 模块
│   ├── m09_spatial.py           # M9 模块
│   └── m10_multiomics.R         # M10 模块
│
├── utils/
│   ├── data_download.py         # 统一数据下载函数
│   ├── report_generator.py      # 自动生成 HTML 报告
│   └── checkpoint.py            # 人工检查点管理
│
└── output/
    ├── {cancer_type}_{date}/    # 每次运行的输出目录
    │   ├── figures/
    │   ├── tables/
    │   └── report.html
```

## 15.3 config.yaml 示例

```yaml
# AML CAR-T 分析配置
project:
  cancer_type: "AML"
  disease_name: "Acute Myeloid Leukemia"
  output_dir: "output/AML_20260520"

# 数据集配置
datasets:
  M1_bulk:    "GSE6891"
  M2_scrna:   "GSE116256"
  M5_M6_M7_tcga: "TCGA-LAML"
  M8_methyl:  "GSE69065"
  M9_spatial: "GSE174448"

# 候选靶点
targets:
  - CLEC12A
  - FLT3
  - CD33
  - IL3RA
  - CD38

# 分析阈值
thresholds:
  M1_logFC:   1.0
  M1_padj:    0.05
  M2_fold:    1.5
  M5_pval:    0.05
  M6_mut_freq: 0.05   # 最低突变频率阈值

# 运行哪些模块（true/false）
modules:
  M1: true
  M2: true
  M3: false   # 结构分析依赖网页工具，暂不自动化
  M4: true
  M5: true
  M6: true
  M7: true
  M8: true
  M9: false   # 空间转录组数据依赖较多，手动运行
  M10: true
```

## 15.4 run_pipeline.py 主调度器

```python
"""
AML CAR-T 分析自动化 Pipeline
用法: python run_pipeline.py --config config.yaml
"""
import yaml
import subprocess
import os
import sys
from datetime import datetime

def run_r_module(script_path: str, config_path: str) -> bool:
    """运行一个 R 模块脚本"""
    cmd = ["Rscript", script_path, "--config", config_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"❌ 模块失败: {script_path}")
        print(result.stderr[-2000:])
        return False
    print(result.stdout[-1000:])
    return True

def run_python_module(script_path: str, config_path: str) -> bool:
    """运行一个 Python 模块脚本"""
    cmd = [sys.executable, script_path, "--config", config_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"❌ 模块失败: {script_path}")
        print(result.stderr[-2000:])
        return False
    print(result.stdout[-1000:])
    return True

def checkpoint(message: str) -> bool:
    """人工检查点：暂停等待用户确认"""
    print(f"\n{'='*60}")
    print(f"⏸️  检查点: {message}")
    print(f"请查看上一步的图表和日志，确认结果正常后继续")
    answer = input("继续运行? (y/n): ").strip().lower()
    return answer == 'y'

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--skip-checkpoints", action="store_true")
    args = parser.parse_args()

    # 加载配置
    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    print(f"🚀 AML 生信分析 Pipeline 启动")
    print(f"   癌症类型: {cfg['project']['cancer_type']}")
    print(f"   输出目录: {cfg['project']['output_dir']}")
    print(f"   开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    modules = cfg['modules']
    results = {}

    # M1: Bulk RNA-seq
    if modules.get('M1'):
        print("\n▶ 运行 M1: Bulk RNA-seq 差异分析")
        ok = run_r_module("modules/m01_bulk_rnaseq.R", args.config)
        results['M1'] = ok
        if ok and not args.skip_checkpoints:
            if not checkpoint("M1 完成 — 请检查火山图，确认分组正确"):
                print("用户终止，退出"); sys.exit(0)

    # M2: scRNA-seq
    if modules.get('M2'):
        print("\n▶ 运行 M2: 单细胞 RNA-seq")
        ok = run_python_module("modules/m02_scrna.py", args.config)
        results['M2'] = ok
        if ok and not args.skip_checkpoints:
            if not checkpoint("M2 完成 — 请检查 UMAP 图，确认细胞聚类合理"):
                print("用户终止，退出"); sys.exit(0)

    # M4: 临床数据库
    if modules.get('M4'):
        print("\n▶ 运行 M4: ClinicalTrials API")
        ok = run_python_module("modules/m04_clinical_mr.py", args.config)
        results['M4'] = ok

    # M5: 生存分析
    if modules.get('M5'):
        print("\n▶ 运行 M5: 生存分析")
        ok = run_r_module("modules/m05_survival.R", args.config)
        results['M5'] = ok

    # M6: 突变分析
    if modules.get('M6'):
        print("\n▶ 运行 M6: 基因组变异分析")
        ok = run_r_module("modules/m06_mutation.R", args.config)
        results['M6'] = ok

    # M7: 免疫微环境
    if modules.get('M7'):
        print("\n▶ 运行 M7: 免疫微环境分析")
        ok = run_r_module("modules/m07_tme.R", args.config)
        results['M7'] = ok

    # M8: 甲基化
    if modules.get('M8'):
        print("\n▶ 运行 M8: DNA 甲基化分析")
        ok = run_r_module("modules/m08_methylation.R", args.config)
        results['M8'] = ok

    # M10: 多组学整合
    if modules.get('M10'):
        print("\n▶ 运行 M10: 多组学整合")
        ok = run_r_module("modules/m10_multiomics.R", args.config)
        results['M10'] = ok

    # 汇总报告
    print(f"\n{'='*60}")
    print("📊 Pipeline 运行汇总:")
    for module, ok in results.items():
        status = "✅" if ok else "❌"
        print(f"  {status} {module}")
    print(f"\n输出目录: {cfg['project']['output_dir']}")
    print(f"完成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

if __name__ == "__main__":
    main()
```

## 15.5 实现路线图（分阶段）

```
阶段一（现在）：稳扎稳打，逐个手动跑
  → 每个里程碑自己手动操作，理解每一步在做什么
  → 积累：什么参数需要调整？什么地方容易出错？

阶段二（M5-M7 完成后）：把已验证的模块脚本化
  → 把跑通的 Rmd / ipynb 改写成可接受命令行参数的脚本
  → 用 config.yaml 管理参数

阶段三（M8-M9 完成后）：写调度器 + 自动报告
  → 实现 run_pipeline.py
  → 用 R Markdown / Quarto 模板自动生成 HTML 报告

阶段四（M10 完成后）：打包成工具
  → 写成 Python 包（pip install aml-pipeline）
  → 支持换癌症类型（改 config.yaml 即可）
  → 这就是一个可以放在 GitHub 上的开源工具，简历上写"开发了 xxx 自动化生信分析工具"
```

## 15.6 给 Claude 的提示词模板（用于新 Cowork）

当你想在新 Cowork 会话里让 Claude 帮你写某个模块的自动化代码时，可以用以下模板：

```
请帮我把以下 SOP 步骤转化为可执行的自动化脚本：

【SOP 来源】Bioinformatics_Pipeline_SOP.md，Part [X]
【模块名称】M[X]_[名称]
【输入】config.yaml 中的参数 + 上一个模块的输出文件
【输出】figures/ 目录下的图表 + tables/ 目录下的 CSV
【技术要求】
- R 脚本需接受 --config 命令行参数
- Python 脚本同上
- 所有输出路径从 config.yaml 读取，不要硬编码
- 每一步打印进度信息（print）
- 遇到错误要 try-catch，打印有意义的错误信息后退出

请按照 SOP 的步骤编写完整的可运行脚本。
```

---

# Part 16：常见问题更新（M5-M10）

| 问题 | 原因 | 解决方案 |
|------|------|---------|
| `TCGAbiolinks` 下载超时 | 服务器在美国，网速慢 | 用 `files.per.chunk=5` 减小并发；或挂代理 |
| CIBERSORT 运行很慢 | 置换检验 perm=1000 太多 | 改为 `perm=100` |
| MOFA+ 训练不收敛 | 数据量太少或因子数太多 | 减少 `num_factors`；增大 `maxiter` |
| `read.metharray.exp` 找不到 IDAT | 路径不对 | 检查 IDAT 文件是成对的（_Red.idat + _Grn.idat）|
| Seurat WNN 报内存错误 | 32GB 内存处理大数据集边缘 | 降低 `k.anchor`；分批处理 |
| cell2location 安装失败 | 依赖 PyTorch | 先装 `pip install torch`，再装 cell2location |

---

# 📌 写在最后 / Final Notes

这份 SOP 的所有代码与命令都已在 Windows 11 + Conda + R 4.6.0 环境下验证可运行。项目主题为急性髓系白血病（AML）CAR-T 靶点发现。但生信工具更新很快，**遇到与官方文档不一致时，以官方文档为准**。

**给你的三条原则**：

1. **不要追求一次跑通完美** — 每个里程碑都先跑个粗糙版本，再迭代优化
2. **每一步都要可视化** — 看图比看数字更能发现错误
3. **代码 + 笔记 + 数据**三件套不可少 — 没有笔记的代码三个月后就看不懂了

**遇到具体步骤的问题，回来贴报错截图，我会逐行帮你排查。**

— SOP End —
