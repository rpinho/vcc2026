#!/usr/bin/env python3
"""Public leaderboard snapshots and incremental UTC-day reports. No credentials."""
import argparse
import datetime as dt
import fcntl
import json
import math
import os
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parent
URL = 'https://virtualcellchallenge.org/api/leaderboard?get_final=false'
COMPONENTS = ('pds', 'mse', 'nmae', 'fid', 'reach', 'jac')


def entries(payload):
    rows = payload['entries'] if isinstance(payload, dict) else payload
    if not isinstance(rows, list) or not rows:
        raise ValueError('empty or malformed leaderboard')
    if any(not isinstance(r, dict) or not r.get('id') for r in rows):
        raise ValueError('entry without submission id')
    return rows


def score(row):
    value = row.get('score_avg')
    return float(value) if isinstance(value, (float, int)) and math.isfinite(value) else -math.inf


def compatible(row):
    return (row.get('partition') == 'val' and row.get('panel_id') == 'vcc2026-val-1'
            and row.get('anchor_version') == 'vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4')


def top_teams(rows):
    result = {}
    for row in sorted(filter(compatible, rows), key=score, reverse=True)[:100]:
        result.setdefault(row.get('team_id') or row.get('team_name'), row)
    return result


def clean(value):
    return str(value or '').replace('\n', ' ').replace('|', '\\|').replace('`', "'")[:2000]


def describe(row):
    return f"{clean(row.get('team_name'))} / {clean(row.get('model_name'))} (`{row['id']}`, {score(row):.5f})"


def diff_report(now, old, new, previous):
    before = {r['id']: r for r in old}
    after = {r['id']: r for r in new}
    added = [r for key, r in after.items() if key not in before]
    changed = [(before[key], r) for key, r in after.items() if key in before and
               any(before[key].get(k) != r.get(k) for k in ('description', 'model_name', 'score_avg', *(f'score_{c}' for c in COMPONENTS)))]
    old_best = {}
    for row in old:
        if compatible(row):
            key = row.get('team_id') or row.get('team_name')
            old_best[key] = max(old_best.get(key, -math.inf), score(row))
    new_top = top_teams(new)
    old_top = top_teams(old)
    jumps = []
    for row in new:
        key = row.get('team_id') or row.get('team_name')
        baseline = old_best.get(key)
        if compatible(row) and baseline is not None and abs(score(row) - baseline) >= .02:
            jumps.append((score(row) - baseline, row))
    ranked = sorted(filter(compatible, new), key=score, reverse=True)
    lines = [f"## {now.isoformat()} UTC", '', f"Source: {URL}",
             f"Previous snapshot: {previous.name if previous else 'none (initial baseline)'}.",
             f"Entries: {len(new)}; added: {len(added)}; changed: {len(changed)}; removed: {len(set(before)-set(after))}.",
             f"Comparable panel/anchor entries: {len(ranked)}. Top 100 means entries, then unique teams."]
    if ranked:
        lines += [f"Leader: {describe(ranked[0])}."]
    lines += ['', '### New or changed descriptions', '']
    descriptions = [(None, r) for r in added if r.get('description')] + [(a,b) for a,b in changed if a.get('description') != b.get('description')]
    lines += [f"- {describe(b)}: {clean(b.get('description')) or '(removed)'}" +
              (f"; previous: {clean(a.get('description'))}" if a else '') for a,b in descriptions] or ['None.']
    lines += ['', '### Score changes (absolute delta >= 0.02 vs previous team best)', '']
    lines += [f"- {delta:+.5f}: {describe(r)}" for delta,r in sorted(jumps, key=lambda x: -abs(x[0]))] or ['None.']
    lines += ['', '### Teams entering the top 100 entries', '']
    lines += [f"- {describe(row)}" for key,row in new_top.items() if key not in old_top] or ['None.']
    lines += ['', '### New submissions', '']
    lines += [f"- {describe(r)}" for r in sorted(added, key=score, reverse=True)] or ['None.']
    lines += ['', '### Changed existing submissions', '']
    lines += [f"- {describe(b)}; previous score {score(a):.5f}" for a,b in changed] or ['None.']
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, help='offline replay instead of public fetch')
    parser.add_argument('--output-dir', type=Path, default=ROOT)
    args = parser.parse_args()
    root = args.output_dir
    raw = root / 'raw'
    raw.mkdir(parents=True, exist_ok=True)
    with (root / '.monitor.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print('monitor already running', flush=True)
            return
        now = dt.datetime.now(dt.timezone.utc).replace(microsecond=0)
        if args.input:
            payload = json.loads(args.input.read_text())
        else:
            request = urllib.request.Request(URL, headers={'User-Agent':'VCC-public-leaderboard-monitor/1.0'})
            with urllib.request.urlopen(request, timeout=45) as response:
                payload = json.load(response)
        new = entries(payload)  # validate before publishing anything
        candidates = sorted(raw.glob('lb-*.json'), key=lambda p: p.name)
        previous = candidates[-1] if candidates else None
        old = entries(json.loads(previous.read_text())) if previous else []
        stamp = now.strftime('%Y-%m-%dT%H%M%SZ')
        target = raw / f'lb-{stamp}.json'
        temporary = target.with_suffix('.tmp')
        temporary.write_text(json.dumps(payload, ensure_ascii=False, separators=(',',':')))
        temporary.replace(target)
        report = root / f'DAILY-{now.date()}.md'
        content = report.read_text() if report.exists() else f'# Public leaderboard monitor — {now.date()}\n\n'
        temporary = report.with_suffix('.tmp')
        temporary.write_text(content + '\n' + diff_report(now, old, new, previous))
        temporary.replace(report)
        print(f'{now.isoformat()} saved {target.name}: {len(new)} entries; report {report.name}', flush=True)


if __name__ == '__main__':
    try:
        os.nice(10)
    except OSError:
        pass
    main()
