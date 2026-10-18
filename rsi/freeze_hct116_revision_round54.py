"""Freeze count-justified Phase1 source revision before any revised prediction."""
from datetime import datetime, timezone
from pathlib import Path
import runtime as rt
import numpy as np
import hct116_revision_round54 as revision


def main():
    assert not revision.PLAN.exists()
    prior = Path('artifacts/HCT116_extension_round53')
    audit = revision.load(prior/'audit.json')
    assert audit['status'].startswith('PASS') and not audit['source_activated']
    assert audit['plan_sha256'] == rt.sha(prior/'plan.json')
    assert audit['receipt_sha256'] == rt.sha(prior/'receipt.json')
    assert audit['bank_sha256'] == rt.sha(prior/'bank.npz')
    assert audit['blocks_sha256'] == rt.sha(prior/'extension_blocks.npz')
    assert audit['auditor_sha256'] == rt.sha('audit_hct116_extension_round53.py')
    for p, h in revision.load(prior/'plan.json')['sha256'].items():
        assert rt.sha(p) == h, p
    manifest = revision.load('artifacts/tasks/manifest.json')
    shared = rt.ROOT.parent/'public_data/model_v4_signature_banks'
    sources, tasks, support = {}, {}, {}
    for task in revision.TASKS:
        entry = manifest[task]
        assert entry['split'] == 'search' and len(entry['targets']) == 79
        tasks[task] = entry
        original = Path('artifacts/cd4_banks') if entry['context'] == 'CD4' else shared
        files = {'K562_GWPS_raw': original/'K562_GWPS_raw.npz', 'HCT116': prior/'bank.npz'}
        if entry['context'] == 'H1':
            files.update({p.stem: p for p in sorted(shared.glob('CD4_*.npz'))})
        sources[task] = {}
        for name, path in files.items():
            dest = revision.ROOT/'banks'/task/(name+'.npz')
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.is_symlink():
                assert dest.resolve() == path.resolve()
            else:
                assert not dest.exists()
                dest.symlink_to(path.resolve())
            sources[task][name] = dict(path=str(path), sha256=rt.sha(path))
        support[task] = {}
        for name in ['K562_GWPS_raw', 'HCT116']:
            with np.load(files[name]) as d:
                eligible = set(d['targets'][d['n_cells'] >= 20])
                support[task][name] = dict(eligible_source_population=len(eligible),
                    panel_supported=sum(t in eligible for t in entry['targets']),
                    panel_abstentions=[t for t in entry['targets'] if t not in eligible])
    assert support['CD4_depth192_search']['HCT116']['panel_supported'] == 78
    assert support['CD4_depth192_search']['HCT116']['panel_abstentions'] == ['HMGXB4']
    files = set(revision.load(prior/'plan.json')['sha256'])
    files.update(str(prior/p) for p in ['plan.json','receipt.json','audit.json','bank.npz','extension_blocks.npz'])
    files.update(['hct116_revision_round54.py','freeze_hct116_revision_round54.py',
        'verify_hct116_revision_round54.py','run_hct116_revision_round54.py',
        'recipes.py','provenance.py','benchmark.py','raw_calibration.py','official_scoring.py',
        'bounded_score.py','bounded_de.py','stat_cache.py','task_io.py','baseline_protocol.py',
        'audit_phase1_calibration.py','artifacts/anchor_policy.json'])
    files.update(str(p) for p in Path('vendor').glob('*.py'))
    for task in revision.TASKS:
        files.update(str(Path('artifacts/scores')/task/p) for p in
            ['anchors.json','anchor_splits.csv','anchor_summary.csv',
             'baseline_mean_response.json','baseline_mean_response_protocol.json'])
        files.add(str(Path('artifacts/tasks')/(task+'.json')))
    for entries in sources.values():
        files.update(r['path'] for r in entries.values())
    # Legacy inputs are bound by their original completed-score identity receipts.
    signatures = {t: __import__('provenance').identity(t) for t in revision.TASKS}
    for signature in signatures.values():
        files.add('artifacts/identities/'+signature+'.json')
    plan = dict(utc=datetime.now(timezone.utc).isoformat(), phase='1',
        revision_id='HCT116_source509_r54', status='FROZEN before revised predictions',
        basis='Round43 count-only6/79 source deficit; round48 full metadata and rounds49-53 frozen original-rule extension. Not chosen by revised scores.',
        scope='Isolated source-adequacy/calibration revision only; no global manifest, bank, original score or gate mutation.',
        source_change='Replace HCT116437 with verified509 union in BOTH primary SEARCH contexts. All original437 rows retained;72 appended. No row selected by response magnitude.',
        centering='Unchanged AtlasShift/v4 source-wide equal eligible-target centering, excluding self-target readout in the common mean. Use all509 union rows (subject to original>=20 eligibility), never destination79 only. Both log2fc and bulk_delta branches.',
        invariant='Same79-target panels, controls, axes, references,>=20 source threshold, family/sealed exclusions,weights,K5622:HCT1161,IV variances/ash,CD4 vote,seed20260910,400 emissions,amplitudes0.6/0.3,promoter rules,five anchors,scorer0.18.0.',
        reuse='Legacy v1/v2b/v3b original results referenced explicitly, never renamed or overwritten. Recompute v4a/v4b for both contexts under distinct score/cache namespaces. New source-population recipe ports carry no new official label.',
        prerequisites='Independent archived-path bitwise head/template and fixed emission audit; source-wide centering and threshold/abstention checks before revised scores.',
        decision='Report revised ranks once. Preserve original failed gate regardless. All-pair ranks alone cannot resolve missing exact external AtlasShift/Sidechain preparation; no automatic OPEN, promotion or forecast.',
        tasks=tasks, sources=sources, support_metadata=support,
        original_implementation_sha256=signatures,
        official_environment=__import__('official_scoring').provenance(),
        sha256={p:rt.sha(p) for p in sorted(files)},
        holdout_scoring_allowed=False, phase2_allowed=False, network_reads=False)
    rt.save(revision.PLAN, plan)
    print('FROZEN',revision.PLAN,rt.sha(revision.PLAN),support,flush=True)


if __name__ == '__main__':
    main()
