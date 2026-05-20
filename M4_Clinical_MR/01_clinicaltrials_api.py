"""
M4-Part1: ClinicalTrials.gov API 数据抓取脚本
=========================================================
【运行方式】在你的 Windows 本地执行，不要在沙箱里跑：

    # 激活任意 conda 环境（base 或 scrna 都行）
    conda activate base
    pip install requests
    cd "D:\Bio-Informatics Case Study\M4_Clinical_MR"
    python 01_clinicaltrials_api.py

输出文件：
    raw_trials.json        — 原始 JSON（备用）
    aml_cart_trials.csv    — 清洗后 CSV（供分析脚本使用）

靶点范围：CD33 / FLT3 / CD123 / CLEC12A(CLL-1) / CD38
API 文档：https://clinicaltrials.gov/data-api/api
"""

import requests
import json
import csv
import time
import os
from datetime import datetime
from collections import Counter

# ─────────────────────────────────────────────────────
# 1. 配置
# ─────────────────────────────────────────────────────

# 脚本所在目录（自动检测）
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

BASE_URL = "https://clinicaltrials.gov/api/v2/studies"

# 靶点 → 搜索词映射
# 每个靶点独立查询，最后按 NCTId 去重，防止重复计数
TARGETS = {
    "CD33":    "AML CAR-T CD33",
    "FLT3":    "AML CAR-T FLT3",
    "CD123":   "AML CAR-T CD123",
    "CLEC12A": "AML CAR-T CLEC12A",
    "CLL-1":   "AML CAR-T CLL-1",
    "CD38":    "AML CAR-T CD38",
}

# 每个靶点最多抓取条数（一般 200 条足够，避免超时）
MAX_PER_TARGET = 200


# ─────────────────────────────────────────────────────
# 2. API 查询
# ─────────────────────────────────────────────────────

def fetch_trials(query: str, max_results: int = 200) -> list:
    """
    调用 ClinicalTrials.gov v2 API，分页抓取所有结果

    v2 API 特点：
    - GET /api/v2/studies?query.term=xxx&pageSize=100
    - 返回 {"studies": [...], "nextPageToken": "...", "totalCount": N}
    - 用 nextPageToken 翻页，没有 nextPageToken 说明到最后一页
    """
    all_studies = []
    next_page_token = None
    page_size = 100  # v2 API 最大 1000，但 100 更稳定

    print(f"  🔍 查询: '{query}'")

    while True:
        params = {
            "query.term": query,
            "pageSize": page_size,
            "format": "json",
        }
        if next_page_token:
            params["pageToken"] = next_page_token

        try:
            resp = requests.get(BASE_URL, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
        except requests.HTTPError as e:
            print(f"  ❌ HTTP 错误: {e}")
            break
        except requests.ConnectionError:
            print("  ❌ 网络连接失败，请检查网络或代理设置")
            break
        except Exception as e:
            print(f"  ❌ 未知错误: {e}")
            break

        studies = data.get("studies", [])
        all_studies.extend(studies)

        total = data.get("totalCount", "?")
        next_page_token = data.get("nextPageToken")

        print(f"     已获取 {len(all_studies)}/{total} 条")

        # 到达上限或最后一页，停止
        if not next_page_token or len(all_studies) >= max_results:
            break

        time.sleep(0.3)  # 礼貌等待，不被限速

    return all_studies


# ─────────────────────────────────────────────────────
# 3. 字段提取（v2 API 嵌套结构解包）
# ─────────────────────────────────────────────────────

def extract_fields(study: dict) -> dict:
    """
    v2 API 的 study 对象结构很深，需要逐层解包
    返回扁平化 dict，方便写 CSV

    主要结构：
    study
    └── protocolSection
        ├── identificationModule   → NCTId, BriefTitle
        ├── statusModule           → OverallStatus, StartDate
        ├── designModule           → Phase, StudyType, Enrollment
        ├── conditionsModule       → Conditions
        ├── armsInterventionsModule → Interventions
        ├── contactsLocationsModule → Countries
        ├── sponsorCollaboratorsModule → LeadSponsor
        ├── descriptionModule      → BriefSummary
        └── outcomesModule         → PrimaryOutcomes
    """
    proto = study.get("protocolSection", {})

    id_mod       = proto.get("identificationModule", {})
    status_mod   = proto.get("statusModule", {})
    design_mod   = proto.get("designModule", {})
    cond_mod     = proto.get("conditionsModule", {})
    arms_mod     = proto.get("armsInterventionsModule", {})
    contacts_mod = proto.get("contactsLocationsModule", {})
    sponsor_mod  = proto.get("sponsorCollaboratorsModule", {})
    desc_mod     = proto.get("descriptionModule", {})
    outcomes_mod = proto.get("outcomesModule", {})

    # 基本信息
    nct_id         = id_mod.get("nctId", "")
    brief_title    = id_mod.get("briefTitle", "")
    official_title = id_mod.get("officialTitle", "")

    # 状态和日期
    status       = status_mod.get("overallStatus", "")
    start_date   = status_mod.get("startDateStruct", {}).get("date", "")
    primary_comp = status_mod.get("primaryCompletionDateStruct", {}).get("date", "")
    completion   = status_mod.get("completionDateStruct", {}).get("date", "")

    # 研究设计
    phase_list = design_mod.get("phases", [])
    phase      = ", ".join(phase_list) if phase_list else "N/A"
    study_type = design_mod.get("studyType", "")
    enrollment = design_mod.get("enrollmentInfo", {}).get("count", "")

    # 适应症
    conditions = "; ".join(cond_mod.get("conditions", []))

    # 干预措施（提取名称和类型）
    interventions = arms_mod.get("interventions", [])
    intv_names = "; ".join(i.get("name", "") for i in interventions)
    intv_types = "; ".join(i.get("type", "") for i in interventions)

    # 国家（去重排序）
    locs = contacts_mod.get("locations", [])
    countries = sorted({loc.get("country", "") for loc in locs if loc.get("country")})
    country_str = "; ".join(countries)

    # 申办方
    lead_sponsor = sponsor_mod.get("leadSponsor", {}).get("name", "")

    # 摘要（截断防止 CSV 太大）
    brief_summary = (desc_mod.get("briefSummary", "") or "")[:600]

    # 主要终点
    primary_outcomes = outcomes_mod.get("primaryOutcomes", [])
    primary_measure = "; ".join(o.get("measure", "") for o in primary_outcomes)

    return {
        "NCTId":                 nct_id,
        "BriefTitle":            brief_title,
        "OfficialTitle":         official_title,
        "OverallStatus":         status,
        "Phase":                 phase,
        "StudyType":             study_type,
        "EnrollmentCount":       enrollment,
        "Conditions":            conditions,
        "InterventionNames":     intv_names,
        "InterventionTypes":     intv_types,
        "StartDate":             start_date,
        "PrimaryCompletionDate": primary_comp,
        "CompletionDate":        completion,
        "Countries":             country_str,
        "LeadSponsor":           lead_sponsor,
        "PrimaryOutcome":        primary_measure,
        "BriefSummary":          brief_summary,
    }


# ─────────────────────────────────────────────────────
# 4. 主流程
# ─────────────────────────────────────────────────────

def main():
    print("=" * 65)
    print("  M4-Part1: ClinicalTrials.gov AML CAR-T 临床试验数据抓取")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 65)

    # Step 1: 按靶点分别查询，用 dict 去重
    all_studies_map = {}   # NCTId → raw study
    target_hit_map  = {}   # NCTId → 匹配到的靶点列表（记录来源）

    for target_name, query_str in TARGETS.items():
        print(f"\n▶ 靶点: {target_name}")
        studies = fetch_trials(query_str, max_results=MAX_PER_TARGET)

        for s in studies:
            nct_id = (s.get("protocolSection", {})
                       .get("identificationModule", {})
                       .get("nctId", "UNKNOWN"))

            if nct_id not in all_studies_map:
                all_studies_map[nct_id] = s
                target_hit_map[nct_id] = []

            # 记录这条试验匹配了哪些靶点（同一试验可能同时匹配多个）
            target_hit_map[nct_id].append(target_name)

    total_unique = len(all_studies_map)
    print(f"\n✅ 去重后共 {total_unique} 条唯一临床试验记录")

    if total_unique == 0:
        print("❌ 没有抓到任何数据，请检查网络连接后重试")
        return

    # Step 2: 保存原始 JSON（备用，可完整回溯）
    raw_json_path = os.path.join(OUTPUT_DIR, "raw_trials.json")
    with open(raw_json_path, "w", encoding="utf-8") as f:
        json.dump(list(all_studies_map.values()), f, ensure_ascii=False, indent=2)
    print(f"✅ 原始 JSON 已保存: raw_trials.json")

    # Step 3: 提取字段，生成扁平化记录
    records = []
    failed_ids = []
    for nct_id, study in all_studies_map.items():
        try:
            rec = extract_fields(study)
            rec["MatchedTargets"] = "; ".join(target_hit_map.get(nct_id, []))
            records.append(rec)
        except Exception as e:
            failed_ids.append(nct_id)
            print(f"  ⚠️  解析 {nct_id} 失败: {e}")

    if failed_ids:
        print(f"  ⚠️  {len(failed_ids)} 条解析失败: {failed_ids[:5]}")

    # Step 4: 保存 CSV（utf-8-sig = 带 BOM，Excel 打开不乱码）
    csv_path = os.path.join(OUTPUT_DIR, "aml_cart_trials.csv")
    if records:
        fieldnames = list(records[0].keys())
        with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)
        print(f"✅ 清洗 CSV 已保存: aml_cart_trials.csv")
        print(f"   {len(records)} 条记录 × {len(fieldnames)} 个字段")

    # Step 5: 控制台统计预览
    print("\n" + "─" * 50)
    print("📊 状态分布 (OverallStatus):")
    for status, cnt in Counter(r["OverallStatus"] for r in records).most_common():
        bar = "█" * min(cnt, 30)
        print(f"  {status:<30s} {cnt:>4d}  {bar}")

    print("\n📊 阶段分布 (Phase):")
    for phase, cnt in Counter(r["Phase"] for r in records).most_common():
        bar = "█" * min(cnt, 30)
        print(f"  {phase:<30s} {cnt:>4d}  {bar}")

    print("\n📊 匹配靶点分布 (MatchedTargets，一条试验可匹配多个):")
    target_counter = Counter()
    for r in records:
        for t in r["MatchedTargets"].split("; "):
            if t:
                target_counter[t] += 1
    for target, cnt in target_counter.most_common():
        bar = "█" * min(cnt, 30)
        print(f"  {target:<15s} {cnt:>4d}  {bar}")

    print("\n🎉 数据抓取完成！")
    print("   下一步：运行 02_analysis_viz.py 生成可视化图表")


if __name__ == "__main__":
    main()
