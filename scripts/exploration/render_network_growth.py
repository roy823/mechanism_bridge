"""Render actual growth counts and link molecule-coordinate viewers."""
import argparse
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



if __name__=='__main__':main()
