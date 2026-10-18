"""Recover selected code/docs from the already inventoried immutable Sidechain tree.

No model data, account access, submission commands or mutable branch reads. Each
attempt is charged before I/O, capped at the tree's blob size, and git-blob checked.
"""
import hashlib
import json
from pathlib import Path
import urllib.request

from transfer_budget import TransferBudget

COMMIT = '2d8f3d17042030ac866c4d568584dc8e5b2bbcc5'
TREE = Path('../research/sources/sidechain_tree.json')
DEST = Path('vendor/sidechain_release')
ALLOCATION = 'Sidechain_pinned_code'
LIMIT = 10_000_000
DOCS = {'LICENSE', 'README.md', 'CLAUDE.md', 'CLAUDE-SCIENCE.md',
        'pyproject.toml', 'configs/datasets.yaml', 'configs/data_sources.yaml',
        'challenges/vcc2026/config.yaml', 'challenges/vcc2026/CLAUDE.md',
        'site/data/submissions.json', 'site/content/models.md'}
TESTS = {'test_count_emitters.py', 'test_dual_moment.py', 'test_shrink_rules.py',
         'test_neighbour_arm.py', 'test_submit_build_alpha_bulk.py'}


def git_blob(data):
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()


def main():
    allocations = json.loads(Path('artifacts/download_allocations.json').read_text())
    assert allocations['new_allocations'][ALLOCATION] == LIMIT
    assert sum(allocations['new_allocations'].values()) == allocations['total_allocated_upper_bound_bytes']
    assert allocations['total_allocated_upper_bound_bytes'] <= allocations['limit_bytes']
    tree = json.loads(TREE.read_text())
    assert tree['sha'] == COMMIT and tree['truncated'] is False
    selected = [x for x in tree['tree'] if x['type'] == 'blob' and
                (x['path'] in DOCS or x['path'].startswith('src/sidechain/') or
                 (x['path'].startswith('tests/') and Path(x['path']).name in TESTS))]
    budget = TransferBudget('artifacts/Sidechain_pinned_code_charges.jsonl', LIMIT)
    budget.bind('commit', COMMIT)
    budget.bind('tree_sha256', hashlib.sha256(TREE.read_bytes()).hexdigest())
    budget.bind('selection', selected)
    receipt = {'commit': COMMIT, 'scope': 'Pinned public code and docs only; not source data or an exact submission reproduction', 'files': {}}
    for item in selected:
        path = item['path']
        assert not Path(path).is_absolute() and '..' not in Path(path).parts
        out = DEST / path
        if out.exists():
            data = out.read_bytes()
            assert len(data) == item['size'] and git_blob(data) == item['sha'], path
        else:
            if item['size']:
                budget.charge(item['size'] + 1, 'GET pinned raw blob ' + path)
                url = f'https://raw.githubusercontent.com/saberhq/sidechain/{COMMIT}/{path}'
                req = urllib.request.Request(url, headers={'Accept-Encoding': 'identity', 'User-Agent': 'VCC-RSI-public-reproduction'})
                # urllib follows HTTP redirects: forbid them so no uncharged second
                # response can be read. raw.githubusercontent.com needs none.
                class NoRedirect(urllib.request.HTTPRedirectHandler):
                    def redirect_request(self, req, fp, code, msg, headers, newurl):
                        return None
                with urllib.request.build_opener(NoRedirect).open(req, timeout=45) as response:
                    assert response.status == 200 and response.geturl() == url
                    data = response.read(item['size'] + 1)
                assert len(data) == item['size'] and git_blob(data) == item['sha'], path
            else:
                data = b''
                assert git_blob(data) == item['sha']
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_bytes(data)
        receipt['files'][path] = {'bytes': len(data), 'git_blob': item['sha'], 'sha256': hashlib.sha256(data).hexdigest()}
        print('verified', path, len(data), flush=True)
    receipt['charged_bytes'] = budget.total
    receipt['allocation_bytes'] = LIMIT
    receipt['status'] = 'PASS pinned code identities only'
    Path('artifacts/sidechain_release_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    budget.close()


if __name__ == '__main__':
    main()
