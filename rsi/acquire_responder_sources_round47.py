"""Bounded public-code inspection; no imports or execution of downloaded code."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request

from transfer_budget import TransferBudget

ROOT = Path('artifacts/responder_sources_round47')
LIMIT = 5_000_000
ALLOCATION = 'Responder_primary_code_round47'
REPO = 'urwahnawaz/cd4_buffering'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def blob(data):
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def fetch(budget, url, cap):
    budget.charge(cap + 1, 'GET ' + url)
    req = urllib.request.Request(url, headers={'Accept-Encoding': 'identity',
        'User-Agent': 'VCC-RSI-primary-source-audit', 'Accept': 'application/vnd.github+json'})
    with urllib.request.build_opener(NoRedirect).open(req, timeout=30) as response:
        assert response.status == 200 and response.geturl() == url
        data = response.read(cap + 1)
    assert len(data) <= cap, 'Response exceeds charged cap'
    return data


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as f:
        f.write(data)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', action='store_true')
    parser.add_argument('--files', nargs='*', default=[])
    args = parser.parse_args()
    allocation = json.loads(Path('artifacts/download_allocations.json').read_text())
    assert allocation['new_allocations'][ALLOCATION] == LIMIT
    assert sum(allocation['new_allocations'].values()) == allocation['total_allocated_upper_bound_bytes'] <= allocation['limit_bytes']
    budget = TransferBudget('artifacts/Responder_primary_code_round47_charges.jsonl', LIMIT)
    budget.bind('repo', REPO)
    if args.index:
        commit_path = ROOT / 'commit.json'
        if not commit_path.exists():
            save(commit_path, fetch(budget, f'https://api.github.com/repos/{REPO}/commits/HEAD', 200_000))
        commit = json.loads(commit_path.read_text())['sha']
        assert len(commit) == 40
        budget.bind('commit', commit)
        tree_path = ROOT / 'tree.json'
        if not tree_path.exists():
            save(tree_path, fetch(budget, f'https://api.github.com/repos/{REPO}/git/trees/{commit}?recursive=1', 600_000))
        tree = json.loads(tree_path.read_text())
        assert tree['sha'] == commit and tree['truncated'] is False
        print('PIN', commit, 'blobs', sum(x['type'] == 'blob' for x in tree['tree']))
        for item in tree['tree']:
            if item['type'] == 'blob':
                print(item['size'], item['path'])
    if args.files:
        commit = json.loads((ROOT / 'commit.json').read_text())['sha']
        tree = json.loads((ROOT / 'tree.json').read_text())
        assert tree['sha'] == commit and tree['truncated'] is False
        budget.bind('commit', commit)
        budget.bind('tree_sha256', sha((ROOT / 'tree.json').read_bytes()))
        entries = {x['path']: x for x in tree['tree'] if x['type'] == 'blob'}
        for name in args.files:
            assert not Path(name).is_absolute() and '..' not in Path(name).parts
            item = entries[name]
            assert item['size'] <= 1_000_000, 'Inspect large notebook separately; no unbounded downloads'
            out = ROOT / 'release' / name
            if out.exists():
                data = out.read_bytes()
            else:
                data = fetch(budget, f'https://raw.githubusercontent.com/{REPO}/{commit}/{name}', item['size'])
                assert len(data) == item['size'] and blob(data) == item['sha'], name
                save(out, data)
            assert len(data) == item['size'] and blob(data) == item['sha'], name
            print('VERIFIED', name, len(data))
        files = {}
        for out in sorted((ROOT / 'release').rglob('*')):
            if out.is_file():
                name = out.relative_to(ROOT / 'release').as_posix()
                data = out.read_bytes()
                assert blob(data) == entries[name]['sha']
                files[name] = dict(bytes=len(data), git_blob=blob(data), sha256=sha(data))
        receipt = dict(utc=datetime.now(timezone.utc).isoformat(), repo=REPO, commit=commit,
            tree_sha256=sha((ROOT / 'tree.json').read_bytes()), files=files,
            charged_bytes=budget.total, allocation_bytes=LIMIT,
            status='PASS immutable source-file identities only; scientific usability requires inspection')
        # Version each acquisition receipt; preserve earlier selections.
        save(ROOT / f'receipt_{len(files)}.json', (json.dumps(receipt, indent=2) + '\n').encode())
    budget.close()


if __name__ == '__main__':
    main()
