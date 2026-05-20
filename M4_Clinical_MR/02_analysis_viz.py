"""
M4-Part1: 临床试验数据分析与可视化
=========================================================
【前置条件】先运行 01_clinicaltrials_api.py，生成 aml_cart_trials.csv

【运行方式】
    conda activate base
    pip install pandas matplotlib seaborn
    python 02_analysis_viz.py

【输出图表】（保存到同目录下）
    fig1_status_pie.png        — 试验状态饼图
    fig2_phase_bar.png         — 试验阶段柱状图
    fig3_target_bar.png        — 各靶点试验数量
    fig4_timeline.png          — 试验启动时间线
    fig5_country_top10.png     — 国家分布 Top10
    summary_table.csv          — 按靶点 × 阶段透视表
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from collections import Counter
from datetime import datetime

# ─────────────────────────────────────────────────────
# 0. 配置
# ─────────────────────────────────────────────────────

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH   = os.path.join(SCRIPT_DIR, "aml_cart_trials.csv")

# 配色（参考 Nature 风格）
PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52",
           "#8172B2", "#937860", "#DA8BC3", "#8C8C8C"]

plt.rcParams.update({
    "figure.dpi":      150,
    "font.size":       11,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


# ─────────────────────────────────────────────────────
# 1. 加载数据
# ─────────────────────────────────────────────────────

def load_data(csv_path: str) -> pd.DataFrame:
    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"找不到 {csv_path}\n"
            "请先运行 01_clinicaltrials_api.py 生成数据文件"
        )

    df = pd.read_csv(csv_path, encoding="utf-8-sig")
    print(f"✅ 加载成功：{len(df)} 条记录，{df.shape[1]} 列")

    # ── 日期列转换 ──
    for col in ["StartDate", "PrimaryCompletionDate", "CompletionDate"]:
        if col in df.columns:
            # 日期格式可能是 "2021-05-01" 或 "2021-05"，都处理
            df[col] = pd.to_datetime(df[col], format="mixed", errors="coerce")

    # ── 招募人数转数值 ──
    if "EnrollmentCount" in df.columns:
        df["EnrollmentCount"] = pd.to_numeric(df["EnrollmentCount"], errors="coerce")

    # ── 阶段标准化 ──
    # v2 API 返回的 Phase 格式：["PHASE1"] 或 ["PHASE1", "PHASE2"]
    # 已在 extract_fields 中 join 为 "PHASE1" 或 "PHASE1, PHASE2"
    # 这里做一次统一映射
    phase_map = {
        "PHASE1":             "Phase I",
        "PHASE2":             "Phase II",
        "PHASE3":             "Phase III",
        "PHASE4":             "Phase IV",
        "PHASE1, PHASE2":     "Phase I/II",
        "PHASE1,PHASE2":      "Phase I/II",
        "PHASE2, PHASE3":     "Phase II/III",
        "PHASE2,PHASE3":      "Phase II/III",
        "EARLY_PHASE1":       "Early Phase I",
        "N/A":                "N/A",
    }
    if "Phase" in df.columns:
        df["PhaseClean"] = (
            df["Phase"]
            .str.strip()
            .str.upper()
            .map(phase_map)
            .fillna(df["Phase"])   # 未匹配的保留原值
        )

    # ── 状态简化 ──
    status_map = {
        "RECRUITING":             "Recruiting",
        "ACTIVE_NOT_RECRUITING":  "Active (Not Recruiting)",
        "COMPLETED":              "Completed",
        "TERMINATED":             "Terminated",
        "WITHDRAWN":              "Withdrawn",
        "SUSPENDED":              "Suspended",
        "NOT_YET_RECRUITING":     "Not Yet Recruiting",
        "ENROLLING_BY_INVITATION":"Enrolling by Invitation",
        "UNKNOWN":                "Unknown",
    }
    if "OverallStatus" in df.columns:
        df["StatusClean"] = df["OverallStatus"].map(status_map).fillna(df["OverallStatus"])

    # ── 拆分 MatchedTargets 为多行（用于按靶点统计）──
    if "MatchedTargets" in df.columns:
        df["TargetList"] = df["MatchedTargets"].str.split("; ")
    else:
        df["TargetList"] = [["Unknown"]] * len(df)

    return df


# ─────────────────────────────────────────────────────
# 2. 各图绘制函数
# ─────────────────────────────────────────────────────

def fig1_status_pie(df: pd.DataFrame, out_dir: str):
    """图1：试验状态饼图"""
    counts = df["StatusClean"].value_counts()

    fig, ax = plt.subplots(figsize=(8, 6))
    wedges, texts, autotexts = ax.pie(
        counts.values,
        labels=counts.index,
        autopct=lambda p: f"{p:.1f}%\n(n={int(p*sum(counts.values)/100)})",
        colors=PALETTE[:len(counts)],
        startangle=140,
        pctdistance=0.75,
        wedgeprops=dict(edgecolor="white", linewidth=1.5),
    )
    for t in autotexts:
        t.set_fontsize(9)

    ax.set_title(
        "AML CAR-T Clinical Trials — Status Distribution\n"
        f"(n={len(df)} unique studies, ClinicalTrials.gov)",
        fontsize=13, fontweight="bold", pad=15
    )

    path = os.path.join(out_dir, "fig1_status_pie.png")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✅ 保存: fig1_status_pie.png")


def fig2_phase_bar(df: pd.DataFrame, out_dir: str):
    """图2：试验阶段柱状图（按状态分色堆叠）"""
    # 阶段顺序
    phase_order = [
        "Early Phase I", "Phase I", "Phase I/II",
        "Phase II", "Phase II/III", "Phase III", "Phase IV", "N/A"
    ]
    phase_order = [p for p in phase_order if p in df["PhaseClean"].values]

    pivot = (
        df.groupby(["PhaseClean", "StatusClean"])
          .size()
          .unstack(fill_value=0)
          .reindex(phase_order, fill_value=0)
    )

    fig, ax = plt.subplots(figsize=(10, 6))
    pivot.plot(kind="bar", ax=ax, stacked=True,
               color=PALETTE[:pivot.shape[1]], edgecolor="white", linewidth=0.5)

    ax.set_xlabel("Clinical Trial Phase", fontsize=12)
    ax.set_ylabel("Number of Studies", fontsize=12)
    ax.set_title(
        "AML CAR-T Clinical Trials — Phase Distribution\n"
        f"(n={len(df)}, colored by status)",
        fontsize=13, fontweight="bold"
    )
    ax.legend(title="Status", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    path = os.path.join(out_dir, "fig2_phase_bar.png")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✅ 保存: fig2_phase_bar.png")


def fig3_target_bar(df: pd.DataFrame, out_dir: str):
    """图3：各靶点临床试验数量（一条试验可匹配多个靶点）"""
    target_counter = Counter()
    for targets in df["TargetList"]:
        for t in targets:
            if t and t != "Unknown":
                target_counter[t] += 1

    if not target_counter:
        print("  ⚠️  MatchedTargets 列为空，跳过图3")
        return

    targets_sorted = sorted(target_counter.keys(), key=lambda x: -target_counter[x])
    counts = [target_counter[t] for t in targets_sorted]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(targets_sorted, counts,
                  color=PALETTE[:len(targets_sorted)], edgecolor="white", linewidth=0.5)

    # 数值标签
    for bar, cnt in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.3,
                str(cnt), ha="center", va="bottom", fontsize=10)

    ax.set_xlabel("Target", fontsize=12)
    ax.set_ylabel("Number of Studies", fontsize=12)
    ax.set_title(
        "AML CAR-T Clinical Trials — by Target\n"
        "(one study may match multiple targets)",
        fontsize=13, fontweight="bold"
    )
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    path = os.path.join(out_dir, "fig3_target_bar.png")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✅ 保存: fig3_target_bar.png")


def fig4_timeline(df: pd.DataFrame, out_dir: str):
    """图4：试验启动年份时间线"""
    df_with_date = df.dropna(subset=["StartDate"]).copy()
    if len(df_with_date) == 0:
        print("  ⚠️  StartDate 全为空，跳过图4")
        return

    df_with_date["StartYear"] = df_with_date["StartDate"].dt.year
    year_status = (
        df_with_date.groupby(["StartYear", "StatusClean"])
                    .size()
                    .unstack(fill_value=0)
    )

    # 只保留 2010 年以后（更有代表性）
    year_status = year_status[year_status.index >= 2010]

    fig, ax = plt.subplots(figsize=(12, 6))
    year_status.plot(kind="bar", ax=ax, stacked=True,
                     color=PALETTE[:year_status.shape[1]], edgecolor="white", linewidth=0.3)

    ax.set_xlabel("Start Year", fontsize=12)
    ax.set_ylabel("Number of Studies", fontsize=12)
    ax.set_title(
        "AML CAR-T Clinical Trials — Launch Timeline (2010–present)\n"
        "(colored by status)",
        fontsize=13, fontweight="bold"
    )
    ax.legend(title="Status", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    path = os.path.join(out_dir, "fig4_timeline.png")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✅ 保存: fig4_timeline.png")


def fig5_country_top10(df: pd.DataFrame, out_dir: str):
    """图5：国家分布 Top 10"""
    country_counter = Counter()
    for countries_str in df["Countries"].dropna():
        for c in str(countries_str).split("; "):
            c = c.strip()
            if c:
                country_counter[c] += 1

    if not country_counter:
        print("  ⚠️  Countries 列为空，跳过图5")
        return

    top10 = country_counter.most_common(10)
    labels = [x[0] for x in top10]
    counts = [x[1] for x in top10]

    fig, ax = plt.subplots(figsize=(9, 6))
    bars = ax.barh(labels[::-1], counts[::-1],
                   color=PALETTE[0], edgecolor="white", linewidth=0.5)

    for bar, cnt in zip(bars, counts[::-1]):
        ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height() / 2,
                str(cnt), va="center", fontsize=10)

    ax.set_xlabel("Number of Studies", fontsize=12)
    ax.set_title(
        "AML CAR-T Clinical Trials — Top 10 Countries",
        fontsize=13, fontweight="bold"
    )
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    path = os.path.join(out_dir, "fig5_country_top10.png")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✅ 保存: fig5_country_top10.png")


def save_summary_table(df: pd.DataFrame, out_dir: str):
    """透视表：靶点 × 阶段"""
    # 展开 TargetList（一行变多行）
    exploded = df.explode("TargetList")
    exploded = exploded[exploded["TargetList"].notna() & (exploded["TargetList"] != "Unknown")]

    if len(exploded) == 0:
        print("  ⚠️  没有靶点数据，跳过透视表")
        return

    pivot = (
        exploded.groupby(["TargetList", "PhaseClean"])
                .size()
                .unstack(fill_value=0)
    )
    pivot["Total"] = pivot.sum(axis=1)

    path = os.path.join(out_dir, "summary_table.csv")
    pivot.to_csv(path, encoding="utf-8-sig")
    print(f"  ✅ 保存: summary_table.csv")
    print("\n📋 靶点 × 阶段透视表：")
    print(pivot.to_string())


# ─────────────────────────────────────────────────────
# 3. 主流程
# ─────────────────────────────────────────────────────

def main():
    print("=" * 65)
    print("  M4-Part1: AML CAR-T 临床试验数据分析与可视化")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 65)

    # 加载数据
    df = load_data(CSV_PATH)

    print("\n▶ 生成图表...")
    fig1_status_pie(df, SCRIPT_DIR)
    fig2_phase_bar(df, SCRIPT_DIR)
    fig3_target_bar(df, SCRIPT_DIR)
    fig4_timeline(df, SCRIPT_DIR)
    fig5_country_top10(df, SCRIPT_DIR)

    print("\n▶ 生成透视表...")
    save_summary_table(df, SCRIPT_DIR)

    print("\n🎉 分析完成！共生成 5 张图 + 1 张透视表")
    print(f"   文件位置: {SCRIPT_DIR}")
    print("\n   下一步：打开 M4_Analysis.ipynb 查看完整 Notebook 版本")


if __name__ == "__main__":
    main()
