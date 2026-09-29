#!/usr/bin/env python3
"""Generate evidence-separated pilot report and path plots from actual run artifacts."""
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from mechbridge.event_graph import geometry_mol,graph_smiles,contract_atom_relays
from mechbridge.chemistry import graph_pair_key

ROOT=Path(__file__).resolve().parents[1]


def main():
    benchmark=json.loads((ROOT/"reports/graph_bridge_benchmark/metrics.json").read_text())
    inventory=json.loads((ROOT/"reports/local_download_inventory.json").read_text(encoding="utf-8"))
    events=[]; paired=[]; feedback_records=[]
    for path in sorted((ROOT/"reports/feasibility").glob("*/verification.json")):
        d=json.loads(path.read_text()); folder=path.parent
        row={k:d.get(k) for k in ["event_id","status","physical_event_verified","verified_pair","elapsed_seconds","gradient_evaluations","error","expected_endpoint_match","event_outcome"]}
        row["imaginary_frequencies_cm-1"]=d.get("ts",{}).get("frequencies_cm-1",[])
        row["arrow_hypotheses"]=d.get("symbolic",{}).get("arrows",[])
        row["orbital_low_overlap"]=d.get("symbolic",{}).get("tracking_low_overlap")
        if (folder/"orbital_ambiguity.json").exists():
            ambiguity=json.loads((folder/"orbital_ambiguity.json").read_text())
            row["near_degenerate_arrow_sets"]=len(ambiguity["last_step_alternative_hypotheses"])
            row["minimum_individual_squared_overlap"]=min(x["minimum_individual_squared_overlap"] for x in ambiguity["adjacent_frames"])
            row["minimum_occupied_subspace_squared_overlap"]=min(x["minimum_occupied_subspace_squared_overlap"] for x in ambiguity["adjacent_frames"])
        stability_path=folder/"rigid_rotation_check/result.json"
        if stability_path.exists(): row["representation_stability"]=json.loads(stability_path.read_text())
        if (folder/"ibo/electronic_report.json").exists():
            e=json.loads((folder/"ibo/electronic_report.json").read_text())
            energy=np.array([x["energy_hartree"] for x in e["frames"]])*27.211386245988
            arrays=[np.load(folder/f"ibo/frame_{i:04d}.npz") for i in range(len(energy))]
            xyz=np.stack([x["positions_A"] for x in arrays])
            distance=np.r_[0,np.cumsum(np.linalg.norm(np.diff(xyz,axis=0).reshape(len(xyz)-1,-1),axis=1))]
            populations=np.stack([x["iao_atom_populations_e"] for x in arrays])
            moving=sorted({x["orbital"] for x in row["arrow_hypotheses"]})
            fig,axes=plt.subplots(2,1,figsize=(8,6))
            axes[0].plot(distance,energy-energy[0],"o-",ms=3)
            axes[0].set_ylabel("Electronic energy relative to reactant (eV)")
            axes[0].set_title(d["event_id"]+": verified connection, automated orbital hypotheses")
            z=d["source"]["atomic_numbers"]
            # Occupied-IBO atomic population changes are evidence, not arrow labels.
            for orbital in moving[:4]:
                changes=np.abs(populations[-1,:,orbital]-populations[0,:,orbital])
                atom=int(np.argmax(changes))
                axes[1].plot(distance,populations[:,atom,orbital],"o-",ms=3,
                             label=f"IBO {orbital}, atom {atom+1} (Z={z[atom]})")
            axes[1].set_xlabel("Accumulated sampled Cartesian displacement (Angstrom)")
            axes[1].set_ylabel("IAO atomic population (electrons)")
            if moving: axes[1].legend(fontsize=8)
            fig.text(.5,.005,"Sampled IRC plus relaxed endpoints; orbital identity ambiguity requires review",ha="center",fontsize=8)
            fig.tight_layout(rect=[0,.025,1,1]); fig.savefig(folder/"path_evidence.png",dpi=180); plt.close(fig)
            if d.get("physical_event_verified") and "symbolic" in d:
                atoms=read(folder/"electronic_frames.xyz",index=":")
                r=graph_smiles(geometry_mol(atoms[0].numbers,atoms[0].positions,d["charge"]))
                p=graph_smiles(geometry_mol(atoms[-1].numbers,atoms[-1].positions,d["charge"]))
                hypotheses=d["symbolic"]
                hypothesis_set=[hypotheses]
                if (folder/"orbital_ambiguity.json").exists():
                    hypothesis_set=json.loads((folder/"orbital_ambiguity.json").read_text())["last_step_alternative_hypotheses"]
                if stability_path.exists():
                    rotated=json.loads(stability_path.read_text())["rotated_arrow_hypotheses"]
                    rotated["evidence_origin"]=str(stability_path.relative_to(ROOT))
                    hypothesis_set.append(rotated)
                for hypothesis in hypothesis_set:
                    hypothesis["atom_relay_normalization"]=contract_atom_relays(hypothesis["arrows"])
                paired.append({"schema_version":"feasibility-pilot/1","event_id":"pilot:"+d["event_id"],
                    "chemical_event_group":graph_pair_key(r,p),"path_instance":d["event_id"],
                    "system":{k:d[k] for k in ["charge","multiplicity","environment"]},
                    "symbolic":{"reactant_smiles":r,"product_smiles":p,
                                "arrow_hypotheses":hypothesis_set,"atom_index_base":0,
                                "alternative_enumeration_complete":False,
                                "equivalence_independently_certified":False,
                                "tracking_low_overlap":hypotheses.get("tracking_low_overlap",False)},
                    "physical":{"method":d["method"],"basis":d["basis"],"ts":d["ts"],
                                "endpoints":d["endpoints"],"atomic_numbers":atoms[0].numbers.tolist(),
                                "ordered_path":str((folder/"ordered_path.xyz").relative_to(ROOT)),
                                "electronic_features":str((folder/"ibo").relative_to(ROOT))},
                    "validation":{"physical_event_verified":True,"arrows":"automated_hypotheses_require_independent_review",
                                  "unique_mechanism_certified":False,"verified_pair":False,
                                  "expected_endpoint_match":d["expected_endpoint_match"]},
                    "provenance":{"source_record":d["source"],"verification_report":str(path.relative_to(ROOT))}})
                feedback_records.append({"physical_event_id":"pilot:"+d["event_id"],
                    "query_id":d["source"].get("parent_event_id",d["event_id"]),
                    "query_endpoints":[d["source"]["reactant_smiles"],d["source"]["product_smiles"]],
                    "observed_endpoints":[r,p],
                    "relation":"compatible_graph_pair" if d["expected_endpoint_match"] else "incompatible_query_event_pair",
                    "physical_event_retained":True,"not_a_negative_reaction_label":True,
                    "representation":"endpoint_graph_query; curved-arrow labels remain unreviewed",
                    "evidence":str(path.relative_to(ROOT)),"novel_chemistry_claim":False})
        events.append(row)
    (ROOT/"data/processed/verified_event_pilot.jsonl").write_text(
        "".join(json.dumps(r,ensure_ascii=False)+"\n" for r in paired),encoding="utf-8")
    (ROOT/"data/processed/physics_feedback_records.jsonl").write_text(
        "".join(json.dumps(r,ensure_ascii=False)+"\n" for r in feedback_records),encoding="utf-8")
    recovery=[]
    for p in (ROOT/"reports/seed_recovery").glob("*/comparison.json"):
        recovery.extend(json.loads(p.read_text())["results"])
    result={"core_data_all_verified":inventory["all_verified"],"graph_benchmark":benchmark,
            "generated_at_utc":datetime.now(timezone.utc).isoformat(),
            "pending_events":[e["event_id"] for e in events if e["status"] in
                              ["started","refining_source_ts","integrating_irc","computing_electronic_path"]],
            "quantum_events":events,"seed_recovery":recovery,
            "physically_verified_events":sum(bool(x["physical_event_verified"]) for x in events),
            "independently_reviewed_arrow_pairs":0,
            "physics_feedback_records":feedback_records,
            "conclusion":"Graph/geometric complementarity and a small physical-evidence pipeline are tested; full arrow-level bidirectional learning and iterative physics-feedback benefit remain unproven"}
    feedback=ROOT/"reports/physics_feedback_probe/result.json"
    if feedback.exists(): result["physics_feedback_probe"]=json.loads(feedback.read_text())
    (ROOT/"reports/feasibility_results.json").write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    m=benchmark["metrics"]
    lines=["# 可行性运行结果", "", "以下数值由实际输出自动汇总；图层代理任务不等于完整弯箭头机理学习。", "",
        f"核心原始数据全部校验：{inventory['all_verified']}；图层样本 {benchmark['events']} 条，{benchmark['formula_groups']} 个分子组成分组，5 折交叉验证。", "",
        "| 对照 | 基线 | 加入跨表示信息 |", "|---|---:|---:|",
        f"| 反应相关原子对 TS 距离 MAE / Å | {m['forward_R']['changed_pair_distance_MAE_A']:.4f} | {m['forward_R_graph_edits']['changed_pair_distance_MAE_A']:.4f} |",
        f"| 净变键识别 F1 | {m['inverse_R']['changed_pair_F1']:.4f} | {m['inverse_R_TS']['changed_pair_F1']:.4f} |",
        f"| 完整事件净变键全部预测正确比例 | {m['inverse_R']['exact_event_edit_fraction']:.4f} | {m['inverse_R_TS']['exact_event_edit_fraction']:.4f} |", "",
        "正向增加的是产物图导出的净键变化，不证明完整箭头超出 R/P 图的额外收益。反向预测净变键，不是轨道到箭头的解码。", "",
        "## 新计算的真实事件", "", "计算层级：气相闭壳层 ωB97X/6-31G(d)。重新优化 TS、检查振动、运行 Sella 双向 IRC，并验证端点图。此层级与 RGD1 原始 B3LYP-D3/TZVP 分开保存。", "",
        "| 事件 | 物理连接验证 | 预期端点匹配 | 状态 | 自动箭头数 | 低轨道重叠 |", "|---|---|---|---|---:|---|"]
    for e in events:
        lines.append(f"| {e['event_id']} | {e['physical_event_verified']} | {e['expected_endpoint_match']} | {e['status']} | {len(e['arrow_hypotheses'])} | {e['orbital_low_overlap']} |")
    lines += ["", "箭头来自 Lewis 电子对容量约束与 IBO 跟踪的自动假设；尚无独立人工真值。电子对守恒是构造约束，不能作为独立准确率。低重叠和近简并候选必须保留。", "",
        "## 留出模型初猜的物理检验", "", "参考 TS 仅用于事后评价；两种初猜均使用同一 DFT 势与梯度预算。", "",
        "| 事件 | 初猜 | 找回参考鞍点 | 梯度次数 | 状态 |", "|---|---|---|---:|---|"]
    for r in recovery:
        lines.append(f"| {r['event_id']} | {r['seed_method']} | {r.get('recovered_reference_saddle',False)} | {r['gradient_evaluations']} | {r['status']} |")
    lines += ["", "本次寻鞍点对照样本很少，不能外推为总体成功率或加速比。", "",
        "## 结论边界", "", "已运行图层双向基线、真实物理连接与电子结构分析，并检验留出模型初猜。尚未证明完整箭头表示的额外学习收益，也尚未证明新增物理证据的多轮反馈能改善独立测试事件。单个模型或轨道表示给出的自洽结果不等于实验机理真值。", "",
        "复现方法见 `docs/FEASIBILITY_RUN_zh.md`；原始 JSON、轨迹、波函数特征和模型均保留在 `reports/` 下。"]
    if feedback.exists():
        f=result["physics_feedback_probe"]
        lines += ["", "## 新物理证据的更新探针", "", f"状态：{f['status']}；独立组成组数：{f['groups']}。"]
        if f["status"]=="completed":
            lines += [f"按新量化事件留一组成验证，变键相关距离的逐事件平均 MAE：更新前 {f['event_mean_MAE_before_A']:.4f} Å，更新后 {f['event_mean_MAE_after_A']:.4f} Å。",
                "该探针只对新的 DFT 几何做残差校准，源理论层级与目标层级分别保存；样本很少，不能证明多轮主动学习或双向箭头模型的反馈收益。"]
            if not f.get("validation_improved",False):
                lines += ["本次留出误差反而增大，候选更新没有替换原模型。不能把这一结果描述为反馈学习已成功。"]
    for event in events:
        if "representation_stability" in event:
            s=event["representation_stability"]
            lines += ["", "## 表示稳定性检查", "",
                f"事件 {event['event_id']} 的整条路径做同一刚体变换后，最大能量差 {s['max_energy_difference_eV']:.8f} eV；原始箭头集合不变：{s['raw_arrow_set_unchanged']}。",
                f"变换前 {s['baseline_arrow_count']} 支箭头，变换后 {len(s['rotated_arrow_hypotheses']['arrows'])} 支。两者都满足构造性的电子对守恒。",
                "这揭示当前局域化/跟踪/解码表示的歧义，不说明物理反应发生改变。未经审核不能自动认定两组箭头在所需机理语义下等价。"]
            if 'normalized_flow_unchanged' in s:
                lines += [f"候选的受限原子位点中继收缩后，流对应一致：{s['normalized_flow_unchanged']}。它仅收缩一入一出、净电子对变化为零的原子位点，保留原始中继证据，不合并不同源—汇配对，也不消除闭环。",
                    "该约定是在观察本例后提出的，尚需在独立事件上检验，不能作为普适化学等价认证。"]
    lines += ["", "占据子空间重叠也较低的帧不能只归因于轨道置换；还需核查采样间隔、刚体对齐及电子态一致性。当前未执行电子态稳定性认证，自动箭头不能直接充当独立监督真值。"]
    (ROOT/"FEASIBILITY_RESULTS_zh.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print(json.dumps({"events":len(events),"physically_verified":result["physically_verified_events"],"core_data_all_verified":inventory["all_verified"]}))


if __name__=="__main__": main()
