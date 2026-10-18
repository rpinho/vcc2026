"""Isolated Phase1 source-adequacy revision; historical recipes stay immutable.

Only HCT116's source population changes (437 -> 509). Templates, pooling,
variance propagation, amplitudes, promoter rules and emitters are archived v4.
"""
from contextlib import contextmanager
from pathlib import Path
import runtime as rt
import numpy as np
from scipy import sparse
import recipes

ROOT = Path('artifacts/HCT116_revision_round54')
PLAN = ROOT/'plan.json'
TASKS = ['H1_rich_search', 'CD4_depth192_search']
RECIPES = ['v4a', 'v4b']
VARIANTS = {'v4a': 'anchor_single_mean', 'v4b': 'iv_ash_dual_guide'}


def load(path):
    return rt.json.loads(Path(path).read_text())


def check_plan():
    plan = load(PLAN)
    for path, digest in plan['sha256'].items():
        assert rt.sha(path) == digest, path
    for task, entries in plan['sources'].items():
        for name, record in entries.items():
            path = ROOT/'banks'/task/(name+'.npz')
            assert path.is_symlink() and path.resolve() == Path(record['path']).resolve(), path
            assert rt.sha(path) == record['sha256'], path
    index = load('artifacts/frozen_metadata_index.json')
    assert rt.sha(index['archive']) == index['archive_sha256']
    for p, record in index['entries'].items():
        assert rt.sha(p) == record['sha256'], p
    assert load('artifacts/calibration_gate.json')['status'] == 'CLOSED'
    return plan


@contextmanager
def bank_scope(destination, folder):
    """Process-local bindings; clear pooled effects at both namespace boundaries."""
    v4 = recipes.v4
    previous = (v4.BANKS, v4.load_banks, v4.cd4_family, v4.atlas.aligned_cd4)
    old_load = v4.load_banks
    v4._POOL_CACHE.clear()
    try:
        v4.BANKS = Path(folder)
        v4.load_banks = lambda *a, **k: [n for n in old_load(*a, **k)
            if recipes.allowed(n, destination) and n != 'H1']
        if destination == 'CD4':
            v4.cd4_family = lambda *a, **k: None
            # Preserve the archived absent-family float32 roundtrip exactly.
            v4.atlas.aligned_cd4 = lambda path, targets, genes, *a, **k: (
                np.zeros((len(targets), len(genes)), np.float32),
                np.zeros((len(targets), len(genes)), bool))
        yield
    finally:
        v4.BANKS, v4.load_banks, v4.cd4_family, v4.atlas.aligned_cd4 = previous
        v4._POOL_CACHE.clear()


def heads(ctrl, targets, destination, recipe, folder):
    assert recipe in RECIPES and destination in ['H1', 'CD4']
    genes = ctrl.var_names.to_numpy(str)
    with bank_scope(destination, folder):
        q, depths, mean, bulk = recipes.v4.template(ctrl, genes, seed=20260910)
        desired, coverage = recipes.v4.heads(np.array(targets), genes, mean, bulk,
            VARIANTS[recipe], exclude=destination, stats={})
    return q, depths, desired, coverage


def predictions(ctrl, task, recipe, plan):
    entry = plan['tasks'][task]
    q, depths, desired, _ = heads(ctrl, entry['targets'], entry['context'], recipe,
                                 ROOT/'banks'/task)
    for i, target in enumerate(entry['targets']):
        seed = recipes.seed_for(entry['context']+':'+target, 20260910)
        x = recipes.v4.emit(q, desired[0][i], desired[1][i], depths, seed, VARIANTS[recipe])
        assert x.shape == (400, ctrl.n_vars) and np.array_equal(x.sum(1), depths)
        yield target, sparse.csr_matrix(x)
