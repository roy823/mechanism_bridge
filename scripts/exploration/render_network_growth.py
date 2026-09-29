"""Render actual growth counts and link molecule-coordinate viewers."""
import argparse
import html
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path)
    root=p.parse_args().root.resolve();summary=json.loads((root/'summary.json').read_text())
    names={'acetone':'丙酮','propanal':'丙醛','glycolaldehyde':'羟基乙醛','glyceraldehyde':'甘油醛',
           'acetylacetone':'乙酰丙酮','cyclobutanone':'环丁酮','acetamide':'乙酰胺','nitroethane':'硝基乙烷'}
    labels={'geometry':'几何','arrows':'符号','hybrid':'混合'}
    colors={'geometry':'#8593a6','arrows':'#3979b7','hybrid':'#098879'}
    plt.rcParams['font.sans-serif']=['Microsoft YaHei'];plt.rcParams['axes.unicode_minus']=False
    systems=list(names);x=np.arange(len(systems))
    fig,axes=plt.subplots(2,1,figsize=(12,7.8),layout='constrained')
    for i,strategy in enumerate(labels):
        rows={r['start']:r for r in summary['rows'] if r['strategy']==strategy}
        for ax,values in zip(axes,[[len(rows[s]['chemical_pairs']) for s in systems],
                                  [rows[s]['longest_demonstrated_chemical_path'] for s in systems]]):
            bars=ax.bar(x+(i-1)*.24,values,.23,label=labels[strategy],color=colors[strategy])
            ax.bar_label(bars,padding=2,fontsize=8)
    for ax in axes:
        ax.set_xticks(x,[names[s] for s in systems]);ax.set_ylim(0,5)
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
    axes[0].set_ylabel('根连通化学图对');axes[0].legend(ncols=3)
    axes[1].set_ylabel('已展示的连续化学步数')
    fig.suptitle('八体系反应网络探索：实际 MLIP 候选',fontsize=17)
    fig.supxlabel('每体系一个随机种子；基于实际极小值连接；不等同于 DFT 认证或实验产物预测',fontsize=10)
    fig.savefig(root/'growth_overview.png',dpi=180);fig.savefig(root/'growth_overview.pdf');plt.close(fig)
    audit_path=root/'DFT_chain_audit.json'
    audit=json.loads(audit_path.read_text()) if audit_path.exists() else None
    verdict=('独立 DFT 复核尚在进行' if audit is None else
        '原两步链通过独立 DFT 复核' if audit['original_MLIP_chain_preserved'] else
        '原 MLIP 两步链未获完整 DFT 认证，逐边结果见报告')
    rows=''.join('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [names[r['start']],labels[r['strategy']],
        r['attempts'],len(r['root_species']),len(r['chemical_pairs']),r['longest_demonstrated_chemical_path'],r['evaluations']])+'</tr>'
        for r in summary['rows'])
    page=f'''<!doctype html><meta charset="utf-8"><title>反应网络增长 v4</title>
<style>body{{font:16px/1.7 system-ui;max-width:1100px;margin:35px auto;padding:0 20px;background:#f6f8fc;color:#193149}}
table{{border-collapse:collapse;width:100%;background:white}}td,th{{padding:8px;border-bottom:1px solid #dbe3ef;text-align:left}}
img{{width:100%}}.status{{padding:16px;background:#fff4d8;border-radius:8px}}a{{color:#087b93}}</style>
<h1>符号与几何引导的反应网络增长</h1>
<p>{summary['attempts']} 次尝试；{summary['multistep_runs']} 组有多步 MLIP 候选；{summary['runs_expanding_new_species']} 组从新物种继续探索。未更新模型权重。</p>
<p class="status">{verdict}</p><img src="growth_overview.png">
<h2>真实分子结构与路径</h2><ul>
<li><a href="search/propanal_s17/molecules/index.html">丙醛：符号提议与新中间体的几何续探</a></li>
<li><a href="search/cyclobutanone_s17/molecules/index.html">环丁酮：多物种网络与真实 TS 构型</a></li></ul>
<h2>逐次统计汇总</h2><table><tr><th>体系</th><th>方法</th><th>尝试</th><th>物种</th><th>化学图对</th><th>连续步</th><th>几何评估</th></tr>{rows}</table>
<p><a href="RESULTS_zh.md">完整结果、效率检查和限制</a> · <a href="summary.json">统计 JSON</a> · <a href="growth_overview.pdf">图表 PDF</a></p>'''
    (root/'index.html').write_text(page,encoding='utf-8')
    print(root/'index.html')


if __name__=='__main__':main()
