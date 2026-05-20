# M5 生存分析 — 运行指南

## 目录结构

```
M5_Survival/
├── 00_install_packages.R          # 第一步：安装 R 包
├── M5_Survival_Analysis.Rmd       # 第二步：主分析（Knit 出 HTML）
├── README_M5_运行指南.md           # 本文件
├── GDC_data/                      # 自动创建：TCGA 数据缓存（首次下载后保留）
└── figures/                       # 自动创建：输出图片
    ├── KM_FLT3.png
    ├── KM_CD33.png
    ├── KM_IL3RA.png
    ├── KM_CLEC12A.png
    ├── KM_CD38.png
    └── Cox_ForestPlot.png
```

---

## 运行步骤

### 第一步：安装 R 包（仅首次运行）

1. 打开 RStudio
2. 打开 `00_install_packages.R`
3. 全选（Ctrl+A）→ 运行（Ctrl+Enter）
4. 等待安装完成，确认所有包显示 ✅

> ⚠️ 如果 TCGAbiolinks 安装卡住，请检查网络（国内建议挂代理）

### 第二步：Knit 主分析报告

1. 在 RStudio 打开 `M5_Survival_Analysis.Rmd`
2. 点击 **Knit** 按钮（或 Ctrl+Shift+K）
3. 选择 **Knit to HTML**
4. 首次运行会下载 TCGA-LAML 数据（约 10-20 分钟），之后读缓存很快

---

## 常见问题

### Q: TCGAbiolinks 下载很慢？
**A:** TCGA 数据从美国 GDC 服务器下载，国内确实慢。解决方案：
- 挂代理（最推荐）
- 或耐心等待，数据只需下载一次，之后读缓存

### Q: 提示 "package not found"？
**A:** 先运行 `00_install_packages.R` 安装所有依赖

### Q: colData 里找不到 `days_to_death` 列？
**A:** 不同版本 TCGAbiolinks 列名可能略有差异，在 Console 运行：
```r
se <- readRDS("GDC_data/TCGA-LAML.rds")  # 如果已缓存
colnames(colData(se))  # 查看所有列名
```
然后在 Rmd 中对应修改列名即可

### Q: 某个基因找不到（target_idx 长度 < 5）？
**A:** STAR-Counts 数据使用 Ensembl ID，rowData 里的 gene_name 列应包含 HGNC 名称。运行：
```r
head(rowData(se_exp)[, c("gene_id", "gene_name")])
```
确认 gene_name 列存在

---

## 分析方法说明

| 分析 | 方法 | 说明 |
|------|------|------|
| 分组 | 中位数分组 | 高/低表达各 50% 样本 |
| 生存曲线 | Kaplan-Meier | 非参数法，95%CI |
| 组间比较 | Log-rank 检验 | 默认显著性阈值 p < 0.05 |
| 风险回归 | Cox 比例风险模型 | 表达量取 log2(x+1) 标准化 |
| 多变量校正 | 年龄（<60 vs ≥60）+ 性别 | AML 重要临床协变量 |

---

## 输出文件

运行成功后，RStudio 同目录会生成 `M5_Survival_Analysis.html`，包含：
- 全部 5 个靶点的 K-M 曲线（含风险表）
- 单变量 Cox 森林图
- 多变量 Cox 校正结果表
- 综合预后评分排名

---

*M5 生存分析 | AML CAR-T 靶点发现 Pipeline*
