"""Check local report links and actual browser interactions with the installed Edge."""
import json,hashlib
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
from playwright.sync_api import sync_playwright,expect
ROOT=Path(__file__).resolve().parents[2];REPORTS=ROOT/'reports'


class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        for key,value in attrs:
            if key in ('href','src') and value:self.links.append(value)


def main():
    out=REPORTS/'site_checks';out.mkdir(exist_ok=True)
    missing=[];count=0
    pages=list(REPORTS.rglob('*.html'))
    for path in sorted(set(pages)):
        parser=Links();parser.feed(path.read_text(encoding='utf-8'))
        for href in parser.links:
            parts=urlsplit(href)
            if parts.scheme or not parts.path:continue
            target=(path.parent/unquote(parts.path)).resolve();count+=1
            if not target.exists():missing.append(dict(page=path.relative_to(ROOT).as_posix(),href=href))
    manifest=json.loads((REPORTS/'site_manifest.json').read_text())
    sources_ok=all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==sha for p,sha in manifest['sources'].items())
    ui=[];errors=[]
    with sync_playwright() as p:
        browser=p.chromium.launch(executable_path='C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe',
            headless=True,args=['--use-angle=swiftshader','--enable-unsafe-swiftshader'])
        page=browser.new_page(viewport={'width':1440,'height':1050},device_scale_factor=1)
        page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto((REPORTS/'index.html').as_uri());expect(page.get_by_role('heading',level=1)).to_be_visible()
        page.screenshot(path=str(out/'overview.png'),full_page=False)
        (out/'overview_dom.txt').write_text(page.locator('body').aria_snapshot(),encoding='utf-8')
        page.get_by_role('link',name='文献与基准',exact=True).click()
        expect(page.get_by_role('heading',name='现在应该比较什么？',exact=True)).to_be_visible()
        expect(page.get_by_role('heading',name='AIMNetCentral 广元素势：实际速度与选择',exact=True)).to_be_visible()
        page.screenshot(path=str(out/'benchmark.png'),full_page=False);ui.append('overview_to_benchmark_navigation')
        page.get_by_role('link',name='物理事件',exact=True).click()
        page.get_by_label('筛选证据等级',exact=True).select_option('DFT_IRC')
        expect(page.locator('#eventCount')).to_have_text(f"{manifest['DFT_IRC_events']} 条事件记录")
        frame=page.frame_locator('#eventFrame')
        expect(frame.locator('#evidenceBadge')).to_have_text('DFT / IRC 实际坐标')
        expect(frame.locator('#status')).not_to_contain_text('符号')
        expect(frame.locator('#event option')).not_to_contain_text('非起点连通')
        expect(frame.locator('#pathNote')).to_contain_text('双向 IRC')
        expect(frame.locator('#forward')).to_contain_text('eV')
        if frame.locator('canvas').count()<4:raise AssertionError('Missing 3D canvases')
        page.screenshot(path=str(out/'events_dft.png'),full_page=False);ui.append('DFT_filter_actual_coordinates_and_WebGL')
        frame.get_by_role('button',name='跳到 TS',exact=True).click()
        expect(frame.locator('#frameLabel')).to_contain_text('过渡态 TS')
        previous=frame.locator('#frame').input_value()
        frame.get_by_role('button',name='播放构型',exact=True).click()
        expect(frame.locator('#frame')).not_to_have_value(previous)
        frame.get_by_role('button',name='暂停',exact=True).click();ui.append('TS_jump_and_actual_frame_playback')
        page.get_by_label('筛选证据等级',exact=True).select_option('MLIP_DESCENT')
        page.get_by_label('筛选实验阶段',exact=True).select_option('bimolecular_v5')
        page.get_by_label('搜索物理事件',exact=True).fill('formaldehyde_water_o0')
        button=page.locator('.event-button').filter(has_text='完整箭头')
        expect(button).to_have_count(1);button.click()
        expect(frame.locator('#evidenceBadge')).to_have_text('MLIP 双侧下降坐标')
        expect(frame.locator('#relatedChecks')).to_contain_text('原图对通过')
        ui.append('MLIP_coordinates_keep_MLIP_label_with_linked_DFT')
        page.get_by_label('搜索物理事件',exact=True).fill('no_such_reaction_12345')
        expect(page.locator('#emptyState')).to_be_visible();expect(page.locator('#eventFrame')).to_be_hidden()
        ui.append('empty_filter_state')
        page.goto((REPORTS/'bimolecular_v5/search/formaldehyde_water_o0_s17/molecules/index.html').as_uri())
        expect(page.get_by_role('navigation',name='全站导航')).to_be_visible()
        expect(page.locator('#relatedChecks')).to_contain_text('原图对通过');ui.append('legacy_molecular_page_uses_shared_navigation_and_evidence')
        page.set_viewport_size({'width':390,'height':844});page.goto((REPORTS/'index.html').as_uri())
        width=page.evaluate('({page:document.documentElement.scrollWidth,viewport:innerWidth})')
        assert width['page']<=width['viewport']+2,width
        page.screenshot(path=str(out/'mobile.png'),full_page=False);ui.append('mobile_layout_no_horizontal_overflow')
        browser.close()
    report=dict(passed=not missing and sources_ok and not errors,html_pages=len(set(pages)),local_links_checked=count,
        missing_links=missing,source_hashes_unchanged=sources_ok,browser='Installed Edge, headless; isolated Playwright environment',
        built_in_browser='unavailable: runtime sandbox metadata error before browser setup',ui_checks=ui,page_errors=errors)
    (out/'validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=True,indent=2))
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':main()
