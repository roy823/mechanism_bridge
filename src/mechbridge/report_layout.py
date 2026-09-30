"""Shared offline page chrome for both research reports and molecular viewers."""
import html,json,os,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
NAVIGATION=[('overview','研究总览','index.html'),('experiments','实验结果','index.html#experiments'),
            ('events','物理事件','events.html'),('benchmark','文献与基准','benchmark.html'),('data','数据与复现','data.html')]


def relative_link(target,folder):
    return os.path.relpath(target,folder).replace('\\','/')


def prepare_shared_assets():
    out=ROOT/'reports/_site';out.mkdir(exist_ok=True)
    shutil.copy2(ROOT/'assets/report_site/site.css',out/'site.css')
    shutil.copy2(ROOT/'assets/molecular_viewer/3Dmol-min.js',out/'3Dmol-min.js')
    shutil.copy2(ROOT/'assets/molecular_viewer/LICENSE',out/'3Dmol-LICENSE.txt')


def attach_checks(events,folder,checks=None):
    from .report_catalog import quantum_checks
    checks=quantum_checks() if checks is None else checks
    for event in events:
        event['source_href']=relative_link(ROOT/event['source_result'],folder)
        event['related_checks']=[dict(label=q['id']+' · '+q['label'],href=relative_link(ROOT/q['file'],folder))
            for q in checks if q['parent']==event['source_result'] or q['file']==event['source_result']]


def molecular_document(payload,folder):
    base=relative_link(ROOT/'reports',folder)+'/'
    links=''.join(f'<a href="{html.escape(base+url)}">{label}</a>' for _,label,url in NAVIGATION)
    nav=f'<header class="topbar"><a class="brand" href="{base}index.html">Mechanism Bridge<small>反应事件 · 证据 · 基准</small></a><nav class="nav" aria-label="全站导航">{links}</nav></header>'
    template=(ROOT/'assets/molecular_viewer/viewer.html').read_text(encoding='utf-8')
    return (template.replace('__MOLECULAR_DATA__',json.dumps(payload,ensure_ascii=False).replace('</','<\\/'))
        .replace('__REPORT_NAV__',nav).replace('__SITE_CSS__',base+'_site/site.css')
        .replace('__VIEWER_SCRIPT__',base+'_site/3Dmol-min.js')
        .replace('__LICENSE_LINK__',base+'_site/3Dmol-LICENSE.txt'))
