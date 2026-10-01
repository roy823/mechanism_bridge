"""Downsample embedded trajectory previews while preserving external full paths."""
import argparse
import json
from pathlib import Path


def spaced(start, stop, count):
    if count >= stop-start+1:
        return list(range(start, stop+1))
    if count <= 1:
        return [start]
    return sorted({round(start+i*(stop-start)/(count-1)) for i in range(count)})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('html', type=Path)
    parser.add_argument('--max-frames', type=int, default=61)
    args = parser.parse_args()
    if args.max_frames < 3:
        raise ValueError('At least three preview frames are required')
    path = args.html.resolve()
    text = path.read_text(encoding='utf-8')
    marker = 'const DATA='
    start = text.index(marker)+len(marker)
    payload, consumed = json.JSONDecoder().raw_decode(text[start:])
    before = path.stat().st_size
    removed = 0
    for event in payload['events']:
        frames = event['frames']
        if len(frames) <= args.max_frames:
            continue
        old_ts = event['ts_index']
        left_count = args.max_frames//2+1
        right_count = args.max_frames-left_count+1
        indices = sorted(set(spaced(0, old_ts, left_count)+
                             spaced(old_ts, len(frames)-1, right_count)))
        event['frames'] = [frames[i] for i in indices]
        event['ts_index'] = indices.index(old_ts)
        for change in event['changes']:
            change['distances'] = [change['distances'][i] for i in indices]
        removed += len(frames)-len(indices)
    encoded = json.dumps(payload, ensure_ascii=False, separators=(',', ':'))
    path.write_text(text[:start]+encoded+text[start+consumed:], encoding='utf-8')
    after = path.stat().st_size
    print(json.dumps(dict(events=len(payload['events']), removed_preview_frames=removed,
        before_MB=round(before/1048576, 1), after_MB=round(after/1048576, 1))))


if __name__ == '__main__':
    main()
