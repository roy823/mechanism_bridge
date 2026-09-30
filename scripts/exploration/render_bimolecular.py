"""Render saved real molecular paths and a fixed-denominator campaign chart."""
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.molecular_visuals import render_visuals


def main():
    root=ROOT/'reports/bimolecular_v5'
    summary=json.loads((root/'summary.json').read_text())
    jobs=[]
    for system in dict.fromkeys(r['system'] for r in summary['rows']):
        candidates=[r for r in summary['rows'] if r['system']==system and r['heavy_events'] and r['strategy']=='arrows']
        if candidates:
            jobs.append((ROOT/candidates[0]['file']).parents[2])
    for job in jobs:print(json.dumps(render_visuals(job)),flush=True)
    font=Path('C:/Windows/Fonts/msyh.ttc')
    if font.exists():
        font_manager.fontManager.addfont(str(font));plt.rcParams['font.family']=font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams['axes.unicode_minus']=False
    systems=list(dict.fromkeys(r['system'] for r in summary['rows']))
    names=['甲醛 + 水','甲醛 + 甲醛','甲醛 + 羟基乙醛','甲醛 + 乙烯二醇','Coley 3217','Coley 3216']
    fig,ax=plt.subplots(figsize=(13,5.5),layout='constrained')
    for index,(strategy,label,color) in enumerate([('geometry','纯几何','#879ab0'),('arrows','符号引导','#087f8c'),('hybrid','混合','#e8a239')]):
        counts=[sum(r['heavy_events']>0 for r in summary['rows'] if r['system']==s and r['strategy']==strategy) for s in systems]
        bars=ax.bar(np.arange(len(systems))+(index-1)*.24,counts,width=.23,label=label,color=color)
        ax.bar_label(bars,padding=3)
    ax.set_xticks(np.arange(len(systems)),names)
    ax.set_yticks([0,1,2]);ax.set_ylim(0,2.45)
    ax.set_ylabel('找到跨分子重原子成键的运行数 / 2 个计划取向')
    ax.set_title('双分子探索：实际 MLIP 鞍点与双侧下降事件',fontsize=17,pad=20)
    ax.legend(loc='upper right');ax.spines[['top','right']].set_visible(False)
    fig.supxlabel('初始极小值失败计入计划分母；不等于参考产物命中，且这些事件尚未与初始相遇盆地严格连通。',fontsize=10)
    fig.savefig(root/'bimolecular_overview.png',dpi=180);fig.savefig(root/'bimolecular_overview.pdf');plt.close(fig)


if __name__=='__main__':main()
