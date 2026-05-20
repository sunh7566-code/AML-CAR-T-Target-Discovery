# M4: 临床数据库 API + 孟德尔随机化因果推断

## 文件清单

| 文件 | 说明 |
|------|------|
| `01_clinicaltrials_api.py` | Part1: 抓取 ClinicalTrials.gov 数据 |
| `02_analysis_viz.py` | Part1: 数据分析与可视化 |
| `M4_Part1_ClinicalTrials.ipynb` | Part1: 完整 Jupyter Notebook 版 |
| `M4_Part2_MendelianRandomization.Rmd` | Part2: MR 分析 RMarkdown |
| `aml_cart_trials.csv` | 抓取后的临床试验数据（运行后生成） |
| `raw_trials.json` | 原始 JSON 备份（运行后生成） |

---

## Part 1: ClinicalTrials.gov API 分析

### 步骤 1 — 抓取数据

在 Windows 本地的 **Anaconda Prompt** 或 **Git Bash** 中运行：

```bash
conda activate base
pip install requests pandas matplotlib seaborn

cd "D:\Bio-Informatics Case Study\M4_Clinical_MR"
python 01_clinicaltrials_api.py
```

运行成功后，目录下会出现：
- `raw_trials.json` — 原始数据备份
- `aml_cart_trials.csv` — 清洗后的 CSV，供分析用

### 步骤 2 — 生成可视化图表

```bash
python 02_analysis_viz.py
```

输出图表：
- `fig1_status_pie.png` — 试验状态饼图
- `fig2_phase_bar.png` — 试验阶段柱状图
- `fig3_target_bar.png` — 各靶点试验数量
- `fig4_timeline.png` — 试验启动时间线
- `fig5_country_top10.png` — 国家分布 Top10
- `summary_table.csv` — 靶点 × 阶段透视表

### 步骤 3 — 运行 Jupyter Notebook（可选，含解读）

```bash
conda activate scrna
jupyter notebook M4_Part1_ClinicalTrials.ipynb
```

---

## Part 2: 孟德尔随机化（MR）分析

### 步骤 1 — 安装 R 包（首次运行）

在 RStudio Console 中运行：

```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages(c("ggplot2", "dplyr", "knitr", "kableExtra"))
```

> ⚠️ `TwoSampleMR` 安装较慢（从 GitHub），耐心等待 5-10 分钟

### 步骤 2 — Knit 报告

在 RStudio 中打开 `M4_Part2_MendelianRandomization.Rmd`，点击 **Knit** 按钮。

输出：
- `M4_Part2_MendelianRandomization.html` — 完整 HTML 报告
- `fig_MR_forest_FLT3.png` — 森林图
- `fig_MR_scatter_FLT3.png` — 散点图
- `fig_MR_loo_FLT3.png` — 留一法敏感性图

### 步骤 3 — 切换为真实数据（进阶）

Rmd 中的模拟数据替换为真实 API 调用（需要访问 `api.opengwas.io`）：

1. 取消 `{r mr_online, eval=FALSE}` 中的 `eval=FALSE`
2. 查找真实 AML GWAS ID：访问 https://gwas.mrcieu.ac.uk 搜索 "leukemia"
3. 替换脚本中的 `AML_GWAS_ID <- "ieu-b-4957"`

---

## 常见问题

**Q: `01_clinicaltrials_api.py` 报 ProxyError？**  
A: 在你的 Windows 本地运行，不要在沙箱/服务器里跑。

**Q: TwoSampleMR 安装失败？**  
A: 确认 Rtools45 已安装，然后重试：
```r
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "D:/rtools45/usr/bin", sep=";"))
remotes::install_github("MRCIEU/TwoSampleMR")
```

**Q: Knit 时找不到文件路径？**  
A: 在 RStudio 中，先 Session → Set Working Directory → To Source File Location
