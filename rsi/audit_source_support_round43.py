"""Response-blind source-support comparability certificate; no expression access.

Only NPZ target/count/quality metadata members are decompressed. Whole-file hashes
bind the containers but do not decode response arrays. Counts describe potential
same-target evidence, not effective gene coverage, causal explanations or scores.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import io
import json
from pathlib import Path
import re
import zipfile

import runtime as rt
import numpy as np
import pandas as pd

PLAN = Path('artifacts/source_support_plan_round43_retry.json')
OUT = Path('artifacts/source_support_audit_round43.json')
TABLE = Path('artifacts/source_support_targets_round43.tsv')
NAMES = ['H1', 'HepG2', 'Jurkat', 'K562', 'RPE1', 'K562_GWPS_bulk']
SHARED = '../public_data/model_v4_signature_banks/'
RECIPES = ['v1', 'v2b', 'v3b', 'v4a', 'v4b']
READS = set()


def load(p):
    return json.loads(Path(p).read_text())


def source(p, minimum, *, alias=False, exclude_sealed=False, de=False):
    return dict(path=p, minimum_cells=minimum, alias=alias,
                exclude_sealed=exclude_sealed, de=de,
                family='CD4' if Path(p).stem.startswith('CD4_') else
                       'K562' if Path(p).stem.startswith('K562') else Path(p).stem)


def pools():
    out = {}
    for context in ['H1', 'CD4', 'public_ABC']:
        public = context == 'public_ABC'
        folder = 'artifacts/cd4_sources/' if context == 'CD4' else 'artifacts/sources/'
        legacy = [n for n in NAMES if public or n not in [context, 'RPE1']]
        for recipe in RECIPES:
            if recipe in ['v1', 'v2b']:
                prefix = ('../model/artifacts/signatures/' if recipe == 'v1' else
                          '../model_v2/artifacts/fine_signatures/') if public else folder
                entries = [source(prefix+n+'.npz', 30, alias=True,
                                  exclude_sealed=not public and n == 'H1') for n in legacy]
            elif recipe == 'v3b':
                # The shared prior is not target-specific evidence. This is only
                # residual *candidacy*; z/sign/abundance tests need expression.
                entries = [source('../public_data/K562_GWPS_bulk.h5ad', 30)]
            else:
                base = 'artifacts/cd4_banks/' if context == 'CD4' else SHARED
                names = ['K562_GWPS_raw', 'HCT116']
                if public:
                    names += ['H1', 'HEK293T']
                entries = [source(base+n+'.npz', 20) for n in names]
                if context != 'CD4':
                    if recipe == 'v4b':
                        entries += [source(SHARED+'CD4_'+c+'.npz', 20)
                                    for c in ['Rest', 'Stim8hr', 'Stim48hr']]
                    entries += [source(SHARED+'CD4_DE_statistics.npz', 20, de=True)]
            out[context+'/'+recipe] = entries
    return out


def freeze():
    assert not PLAN.exists() and not OUT.exists() and not TABLE.exists()
    ps = pools()
    files = {s['path'] for entries in ps.values() for s in entries}
    files.update([__file__, 'runtime.py', 'recipes.py', 'vendor/legacy_pipeline.py',
                  'vendor/legacy_v2.py', 'vendor/legacy_profiles.py',
                  'vendor/v4_model.py', 'vendor/atlas_model.py',
                  '../model/pipeline.py', '../model_v2/v2_model.py',
                  '../model_v2/build_v2.py', '../model_v3/build.py',
                  '../model_v3/artifacts/input_identity.json',
                  '../model_v4/model.py', '../model_v4/build.py',
                  '../model_v4/out/prediction_v4a.manifest.json',
                  '../model_v4/out/prediction_v4b.manifest.json',
                  '../data/pert_counts.csv', '../model_v4/artifacts/targets.json',
                  'artifacts/gene_aliases.json',
                  '../model/artifacts/gene_aliases.json',
                  '../model_v2/artifacts/gene_aliases.json',
                  'artifacts/tasks/manifest.json', 'artifacts/benchmark_roles.json',
                  'artifacts/anchor_policy.json', 'artifacts/calibration_gate.json',
                  'artifacts/download_allocations.json',
                  'artifacts/frozen_metadata_index.json',
                  'artifacts/zsp_results_audit_round42.json'])
    files.update(str(p) for task in ['H1_rich_search', 'CD4_depth192_search']
                 for p in Path('artifacts/scores', task).glob('*_s42.json'))
    rt.save(PLAN, dict(utc=datetime.now(timezone.utc).isoformat(), status='FROZEN',
        pools=ps, preserve_sha256={p: rt.sha(p) for p in sorted(files)},
        panels={'H1': 'H1_rich_search', 'CD4': 'CD4_depth192_search',
                'public_ABC': '../data/pert_counts.csv'},
        allowed_npz_members=['targets', 'n_cells', 'available', 'quality_pass'],
        allowed_hdf5_elements=['obs/_index (construct labels)', 'obs/num_cells_filtered'],
        target_matching='v1/v2b exact symbol first, else first same annotation-ID row; others exact symbols. Family count collapses K562 datasets and CD4 conditions. H1 sealed-target rows removed only from local legacy H1, matching recipes.py.',
        counts='Existing recipe thresholds only: legacy30, v4_20. All target-level DE quality flags retained; any qualifying condition is metadata eligibility, not proof of a finite gene effect. No expression-derived robust v3b, IV or shrinkage gates are evaluated.',
        estimands='For all five historical recipes: per-target source-file metadata counts and eligible independent families; full family-count histogram, each family coverage and any/multifamily totals. Compare local panel/local pools, public300/local pools, public300/declared public pools. Never choose a new panel, source, weight or threshold.',
        public_scope='Cached source pools declared by archived builders, not byte-certified submitted historical source banks. Current public-source metadata need not equal submission-time identities. ABC share one public300 list; no hidden context identity is assumed. Sealed sources are count metadata only in the public-pool comparison; never activated locally.',
        verification='NPZ member allowlist; independent pandas-index/vector and scalar row lookup/sets/count histograms agree for every target; fixture exercises aliases, threshold boundary, family collapse and sealed exclusion. Bind all prior ten historical records and frozen metadata; recheck files before output.',
        limits='Descriptive support certificate only. Independent family means distinct declared cell-line family, not statistical independence. Does not explain rank reversal, repair gate, show biological improvement or certify effective gene coverage. v3b is same-target residual candidacy only; shared priors and promoter priors do not count as measured same-target evidence.',
        expression_payloads_decoded=False, network_reads=False, holdout_scored=False,
        phase2_allowed=False, new_task_registered=False))
    print('FROZEN', PLAN, flush=True)


def bindings(plan):
    for p, digest in plan['preserve_sha256'].items():
        assert rt.sha(p) == digest, p


def member(path, name):
    assert name in ['targets', 'n_cells', 'available', 'quality_pass']
    with zipfile.ZipFile(path) as z:
        data = np.load(io.BytesIO(z.read(name+'.npy')), allow_pickle=False)
    READS.add((path, name))
    return data


def metadata(entry, sealed):
    p = entry['path']
    if p.endswith('.h5ad'):
        import h5py
        import anndata as ad
        with h5py.File(p, 'r') as f:
            index = f['obs'].attrs['_index']
            names = ad.io.read_elem(f['obs'][index]).astype(str)
            ns = ad.io.read_elem(f['obs/num_cells_filtered'])
        READS.update([(p, 'obs/'+index), (p, 'obs/num_cells_filtered')])
        labels = np.array(['non-targeting' if 'non-targeting' in s else
                           re.match(r'^[^_]+_(.+)_[^_]+_([^_]+)$', s).group(1)
                           for s in names])
        ns = np.nan_to_num(ns, nan=0)
        assert np.isfinite(ns).all() and (ns >= 0).all()
        keep = (ns >= 30) & (labels != 'non-targeting')
        ts = np.array(sorted(set(labels[keep])))
        ix = pd.Index(ts).get_indexer(labels[keep])
        n = np.bincount(ix, weights=ns[keep], minlength=len(ts))
        direct = Counter()
        for label, count in zip(labels, ns):
            if count >= 30 and label != 'non-targeting':
                direct[str(label)] += float(count)
        assert n.tolist() == [direct[t] for t in ts]
        return ts, n, n >= entry['minimum_cells']
    ts = member(p, 'targets').astype(str)
    n = member(p, 'n_cells')
    assert len(set(ts)) == len(ts) and np.isfinite(n).all() and (n >= 0).all()
    if entry['de']:
        a, q = member(p, 'available'), member(p, 'quality_pass')
        assert a.shape == q.shape == n.shape and n.shape[1] == len(ts)
        valid = (a & q & (n >= entry['minimum_cells'])).any(axis=0)
    else:
        assert n.shape == (len(ts),)
        valid = n >= entry['minimum_cells']
    if entry['exclude_sealed']:
        keep = ~np.isin(ts, list(sealed))
        ts, valid = ts[keep], valid[keep]
        n = n[:, keep] if entry['de'] else n[keep]
    return ts, n, valid


def vector_indices(ts, targets, aliases, use_alias):
    ix = pd.Index(ts).get_indexer(targets)
    if use_alias:
        ids = np.array([aliases.get(t, '') for t in ts])
        for j in np.flatnonzero(ix < 0):
            tid = aliases.get(targets[j], '')
            matches = np.flatnonzero(ids == tid) if tid else []
            if len(matches):
                ix[j] = matches[0]
    return ix


def scalar_index(ts, t, aliases, use_alias):
    for i, s in enumerate(ts):
        if s == t:
            return i
    if use_alias and aliases.get(t):
        for i, s in enumerate(ts):
            if aliases.get(s) == aliases[t]:
                return i
    return -1


def fixture():
    ts = np.array(['EXACT', 'OTHER', 'LOW', 'SEALED'])
    targets = ['EXACT', 'ALIAS', 'LOW', 'ABSENT', 'SEALED']
    aliases = {'EXACT': '1', 'OTHER': '2', 'ALIAS': '2'}
    ix = vector_indices(ts, targets, aliases, True)
    assert ix.tolist() == [scalar_index(ts, t, aliases, True) for t in targets] == [0, 1, 2, -1, 3]
    ns = np.array([30, 31, 29, 100])
    assert [int(ns[i] >= 30) if i >= 0 else 0 for i in ix] == [1, 1, 0, 0, 1]
    kept = ts[ts != 'SEALED']
    assert scalar_index(kept, 'SEALED', aliases, True) == -1
    assert len({'K562', 'K562', 'CD4', 'CD4'}) == 2


def main():
    assert not OUT.exists() and not TABLE.exists(), 'Preserve prior evidence'
    plan = load(PLAN)
    bindings(plan)
    fixture()
    manifest = load('artifacts/tasks/manifest.json')
    sealed = set(manifest['H1_holdout']['targets'])
    public = pd.read_csv('../data/pert_counts.csv').target_gene.astype(str).tolist()
    assert len(public) == len(set(public)) == 300
    preparation_targets = load('../model_v4/artifacts/targets.json')
    assert preparation_targets[:len(public)] == public
    assert len(preparation_targets) == len(set(preparation_targets))
    assert rt.sha('../data/pert_counts.csv') == load('../model_v3/artifacts/input_identity.json')['input_sha256']['data/pert_counts.csv']
    for arm, variant in [('a', 'anchor_single_mean'), ('b', 'iv_ash_dual_guide')]:
        assert load(f'../model_v4/out/prediction_v4{arm}.manifest.json')['variant'] == variant
    frozen = load('artifacts/frozen_metadata_index.json')
    for p, rec in frozen['entries'].items():
        assert rt.sha(p) == rec['sha256'], p
    assert rt.sha(frozen['archive']) == frozen['archive_sha256']
    gate = load('artifacts/calibration_gate.json')
    assert gate['status'] == 'CLOSED' and not gate['phase2_allowed']
    # Verify counts independently without importing recipes/source loaders,
    # because those eagerly decompress expression arrays.
    summaries, rows = {}, []
    for pool_id, entries in plan['pools'].items():
        context, recipe = pool_id.split('/')
        panels = [('public300', public)]
        if context != 'public_ABC':
            task = plan['panels'][context]
            assert manifest[task]['split'] == 'search'
            panels.insert(0, (task, manifest[task]['targets']))
        alias_path = ('../model/artifacts/gene_aliases.json' if recipe == 'v1' else
                      '../model_v2/artifacts/gene_aliases.json') if context == 'public_ABC' else 'artifacts/gene_aliases.json'
        aliases = load(alias_path)
        data = [(entry, metadata(entry, sealed)) for entry in entries]
        de_flags = {e['path']: (member(e['path'], 'available'),
                               member(e['path'], 'quality_pass'))
                    for e in entries if e['de']}
        for panel, targets in panels:
            assert len(targets) == len(set(targets))
            family_arrays = {}
            scalar_families = [set() for _ in targets]
            details = [[] for _ in targets]
            for entry, (ts, n, valid) in data:
                ix = vector_indices(ts, targets, aliases, entry['alias'])
                eligible = np.zeros(len(targets), bool)
                present = ix >= 0
                eligible[present] = valid[ix[present]]
                family = entry['family']
                family_arrays.setdefault(family, np.zeros(len(targets), bool))
                family_arrays[family] |= eligible
                for j, t in enumerate(targets):
                    i = scalar_index(ts.tolist(), t, aliases, entry['alias'])
                    assert i == int(ix[j])
                    if i < 0:
                        ok = False
                        count, flags = None, None
                    elif entry['de']:
                        a, q = de_flags[entry['path']]
                        # No sealed DE branch exists in these fixed pools.
                        flags = [bool(a[c, i]) and bool(q[c, i]) and int(n[c, i]) >= entry['minimum_cells'] for c in range(n.shape[0])]
                        ok, count = any(flags), n[:, i].tolist()
                    else:
                        count, flags = int(n[i]), None
                        ok = count >= entry['minimum_cells']
                    assert ok == bool(eligible[j])
                    if ok:
                        scalar_families[j].add(family)
                    details[j].append(dict(path=entry['path'], family=family,
                        matched_target=str(ts[i]) if i >= 0 else None,
                        n_cells=count, condition_eligibility=flags, eligible=ok))
            counts = np.sum(list(family_arrays.values()), axis=0).astype(int)
            assert counts.tolist() == [len(x) for x in scalar_families]
            histogram = dict(Counter(map(int, counts)))
            assert histogram == {i: sum(len(x) == i for x in scalar_families) for i in set(map(len, scalar_families))}
            by_family = {f: int(v.sum()) for f, v in family_arrays.items()}
            assert by_family == {f: sum(f in x for x in scalar_families) for f in family_arrays}
            key = pool_id+'/'+panel
            summaries[key] = dict(n_targets=len(targets), any_source=int((counts > 0).sum()),
                multiple_families=int((counts >= 2).sum()),
                mean_families=float(counts.mean()), family_histogram=histogram,
                by_family=by_family,
                eligibility_only=True, v3b_residual_candidacy_only=recipe == 'v3b')
            for j, t in enumerate(targets):
                rows.append(dict(pool=pool_id, panel=panel, recipe=recipe, target=t,
                    n_families=int(counts[j]), families=','.join(sorted(scalar_families[j])),
                    source_metadata=json.dumps(details[j], separators=(',', ':'))))
            print(key, summaries[key], flush=True)
    bindings(plan)
    pd.DataFrame(rows).to_csv(TABLE, sep='\t', index=False)
    rt.save(OUT, dict(status='PASS metadata comparability only',
        utc=datetime.now(timezone.utc).isoformat(), plan_sha256=rt.sha(PLAN),
        verifier_sha256=rt.sha(__file__), table_sha256=rt.sha(TABLE),
        rows=len(rows), summaries=summaries,
        decoded_metadata_members=[dict(path=p, member=m) for p, m in sorted(READS)],
        frozen_metadata_members_verified=len(frozen['entries']),
        historical_score_records_unchanged=10,
        expression_payloads_decoded=False, phase2_allowed=False,
        holdout_scored=False, new_task_registered=False, network_reads=False))
    print('PASS', OUT, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    args = parser.parse_args()
    freeze() if args.freeze else main()
