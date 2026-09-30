"""Summarize performance, domain, reactive transfer, and the model choice."""
import hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];REPORT=ROOT/'reports/aimnetcentral_v8'
NAMES={'aimnet2':'AIMNet2','aimnet2-2025':'AIMNet2-2025','aimnet2-nse':'AIMNet2-NSE','aimnet2-rxn':'AIMNet2-rxn'}
URLS={'aimnet2':'https://storage.googleapis.com/aimnetcentral/aimnet2v2/AIMNet2/aimnet2_wb97m_d3_0.pt',
      'aimnet2-2025':'https://storage.googleapis.com/aimnetcentral/aimnet2v2/AIMNet2/aimnet2_2025_b973c_d3_0.pt',
      'aimnet2-nse':'https://storage.googleapis.com/aimnetcentral/aimnet2v2/AIMNet2NSE/aimnet2nse_wb97m_0.pt',
      'aimnet2-rxn':'https://storage.googleapis.com/aimnetcentral/aimnet2v2/AIMNet2rxn/aimnet2_rxn_0.pt'}


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def main():
    cpu=read(REPORT/'cpu7/results.json');gpu=read(REPORT/'gpu6/results.json')
    by_cpu={m['family']:m for m in cpu['models']};by_gpu={m['family']:m for m in gpu['models']}
    reactive={}
    for path in REPORT.glob('reactive_*/results.json'):
        rows=read(path);reactive[rows[0]['family']]=rows
    sn2={read(p)['family']:read(p) for p in REPORT.glob('sn2_*/result.json')}
    rows=[]
    for family in by_cpu:
        c,g=by_cpu[family],by_gpu[family];chno_c=c['systems']['chno_11'];chno_g=g['systems']['chno_11']
        transfer=reactive[family]
        rows.append(dict(family=family,name=NAMES[family],elements=len(c['description']['implemented_species']),
            reference={'aimnet2':'wB97M-D3','aimnet2-2025':'B97-3c','aimnet2-nse':'wB97M-D3, spin-polarized','aimnet2-rxn':'reaction-specialized'}[family],
            cpu_single_ms=1000*chno_c['single']['median_seconds'],
            cpu_batch32_per_second=chno_c['batch32_geometries_per_second'],
            gpu_single_ms=1000*chno_g['single']['median_seconds'],
            gpu_batch32_per_second=chno_g['batch32_geometries_per_second'],
            gpu_peak_MiB=g['peak_cuda_memory_MiB'],
            cpu_fd_hessian_ms=1000*c['finite_difference_hessian']['seconds'],
            cpu_native_hessian_ms=1000*c['native_hessian']['median_seconds'],
            gpu_fd_hessian_ms=1000*g['finite_difference_hessian']['seconds'],
            gpu_native_hessian_ms=1000*g['native_hessian']['median_seconds'],
            reactive_strict_passes=sum(r['status']=='validated_descents' and r['DFT_reference_pair_preserved'] for r in transfer),
            reactive_pair_preserved=sum(r['DFT_reference_pair_preserved'] for r in transfer),
            reactive_evaluations=sum(r['evaluations'] for r in transfer),
            sn2_strict_pass=sn2.get(family,{}).get('strict_self_exchange_passed'),
            sn2_evaluations=sn2.get(family,{}).get('evaluations'),
            broad_pair_supported=c['systems']['flowER_halogenated_organic_pair']['supported']))
    files=[]
    for model in cpu['model_files']:
        family=next(f for f,m in by_cpu.items() if Path(m['description']['model_path']).name==Path(model['file']).name)
        files.append(dict(family=family,url=URLS[family],**model))
    receipt=dict(files=files,registry='aimnet 0.2.0 bundled model_registry.yaml',
        registry_sha256=hashlib.sha256((ROOT/'.venv-explore/Lib/site-packages/aimnet/calculators/model_registry.yaml').read_bytes()).hexdigest())
    (REPORT/'model_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    summary=dict(rows=rows,recommendation=dict(
        broad_closed_shell='aimnet2-2025',high_accuracy_broad='aimnet2',open_shell='aimnet2-nse',
        chno_reactive_baseline='aimnet2-rxn pinned historical artifact',palladium='aimnet2-pd, separate unbenchmarked branch'),
        CPU=cpu,GPU=gpu,reactive_transfer=reactive,SN2=sn2,
        registry_vs_pinned_rxn=cpu['current_rxn_artifact_consistency'],
        pipeline_smoke=read(REPORT/'pipeline_smoke_2025/aimnetcentral_documented_sn2_reactant/geometry/network.json'),
        claims='Reference-assisted local PES checks and timing; not broad chemical accuracy, autonomous discovery, or family-generalization evidence')
    (REPORT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    labels=[r['name'] for r in rows];x=range(len(rows));width=.35
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    axes[0].bar([i-width/2 for i in x],[r['cpu_batch32_per_second'] for r in rows],width,label='CPU · 2 threads')
    axes[0].bar([i+width/2 for i in x],[r['gpu_batch32_per_second'] for r in rows],width,label='RTX 4060')
    axes[0].set(ylabel='32-geometry batch throughput / geometries s⁻¹',xticks=list(x),xticklabels=labels)
    axes[0].legend();axes[0].tick_params(axis='x',rotation=18)
    axes[1].bar(labels,[r['reactive_strict_passes'] for r in rows],color='#17808a')
    axes[1].set(ylabel='Strict passes / 3 DFT-seeded CHNO events',ylim=(0,3.35));axes[1].tick_params(axis='x',rotation=18)
    fig.suptitle('AIMNetCentral: throughput and reference-assisted reactive checks')
    fig.savefig(REPORT/'comparison.png',dpi=180);fig.savefig(REPORT/'comparison.pdf');plt.close(fig)
    lines=['# AIMNetCentral 广元素模型评测 v8','',
        '结论：广元素闭壳层探索首选 `aimnet2-2025`；需要更高参考层级时用 `aimnet2`；自由基/开壳层才启用 `aimnet2-nse`；CHNO 历史结果继续固定原 `aimnet2-rxn`。','',
        '| 模型 | 元素数 | CPU 单点/ms | CPU batch32/s | GPU 单点/ms | GPU batch32/s | GPU峰值/MiB | CHNO严格通过 | SN2 |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---|']
    for r in rows:lines.append(f"| {r['name']} | {r['elements']} | {r['cpu_single_ms']:.2f} | {r['cpu_batch32_per_second']:.0f} | {r['gpu_single_ms']:.2f} | {r['gpu_batch32_per_second']:.0f} | {r['gpu_peak_MiB']:.1f} | {r['reactive_strict_passes']}/3 | {r['sn2_strict_pass'] if r['sn2_strict_pass'] is not None else '不支持 Cl'} |")
    lines+=['','## 如何使用','',
        '- Dimer 的逐步单构型调用：对这些 11 原子体系 CPU 更快；继续用 CPU。',
        '- 批量有限差分 Hessian、NEB 图像或更大分子：GPU 吞吐明显更高。',
        '- 当前 11 原子体系上，GPU 批量有限差分 Hessian 比 GPU 原生二阶导数更快；CPU 原生 Hessian略快。当前工作流保留批量有限差分。',
        '- 四成员 ensemble 约增加四倍计算量；搜索阶段使用 member0，只有筛选/不确定性阶段再计算 ensemble。','',
        '## 反应区检查','',
        '三个模型间共有的 DFT TS/负模辅助事件是丙酮互变、甲醛水合、甲醛二聚。`aimnet2` 与 `aimnet2-2025` 都严格保留 3/3；NSE 保留图对 3/3，但二聚端点未全部通过曲率检查。该检查使用了参考 TS 和负模，不是自主发现。','',
        '官方文档的 Cl⁻ + CH₃Cl 对称 SN2 例子在三个广元素模型上都通过严格检查；AIMNet2-2025 为 210 次评估、约 1.28 s CPU。','',
        '## 重要版本差异','',
        f"当前 AIMNetCentral 注册表对 `aimnet2-rxn` 加入外部 D3；相对历史固定 HF 推理，测试几何的能量差为 {summary['registry_vs_pinned_rxn']['official_registry_default_vs_pinned']['energy_abs_difference_eV']:.4f} eV、最大力差为 {summary['registry_vs_pinned_rxn']['official_registry_default_vs_pinned']['force_max_abs_difference_eV_A']:.4f} eV/Å。关闭 D3 后两者数值完全一致。旧结果不能静默换成注册表默认。",'',
        '## 边界','',
        '- 14 元素范围：H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, I；不含碱金属/碱土金属及除 Pd 外的过渡金属。',
        '- 当前符号动作库仍以中性闭壳层 CHNO 为主。换势能扩大了物理后端，尚未自动扩大电子机理提议覆盖。',
        '- GPU 结果使用现有 CUDA 环境做工程计时；正式部署应建立干净锁定环境。未评测 `torch.compile`。',
        '- 不同模型参考理论层级不同，能量不能混在同一网络或直接拼接训练。','',
        '数据：[summary.json](summary.json)、[model_receipt.json](model_receipt.json)、[comparison.png](comparison.png)。']
    (REPORT/'RESULTS_zh.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
