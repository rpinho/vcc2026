"""Build isolated raw HCT116 moments from certified local records; no model use."""
import argparse
from datetime import datetime, timezone
from pathlib import Path

import runtime as rt
import numpy as np
import pandas as pd

ROOT = Path('artifacts/HCT116_extension_round53')
STAGE = Path('artifacts/HCT116_payload_round52')
PRE = Path('artifacts/HCT116_extension_preflight_round49')
BASE = Path('artifacts/cd4_banks/HCT116.npz')
GENES = Path('artifacts/HCT116_expression_identity_round50/genes.parquet')
PLAN = ROOT/'plan.json'
ROW_KEYS = ['n_cells', 'target_count_sums', 'target_mean_cpm',
            'matched_control_probability', 'matched_control_mean_cpm',
            'delta_cpm_variance', 'target_mean_cpm_variance',
            'matched_control_mean_cpm_variance']


def load(path):
    return rt.json.loads(Path(path).read_text())


def freeze():
    ROOT.mkdir(exist_ok=True)
    assert not PLAN.exists() and not (ROOT/'bank.npz').exists()
    paths = {Path(__file__), Path('audit_hct116_extension_round53.py'), Path('runtime.py'),
             STAGE/'plan.json', STAGE/'receipt.json', STAGE/'audit.json',
             Path('audit_hct116_payload_round52.py'), GENES,
             PRE/'preflight.json', PRE/'targets.tsv', PRE/'target_sample_counts.tsv',
             PRE/'controls.tsv', PRE/'selection.parquet', PRE/'logical_record_ranges.parquet',
             Path('../model_v4/prepare_orion.py'), Path('../model_v4/correct_orion_controls.py'),
             Path('artifacts/frozen_metadata.tar.gz')}
    paths |= {Path(p) for p in load(STAGE/'plan.json')['sha256']}
    paths |= set(Path('artifacts').glob('*banks/*.npz'))
    paths |= set(Path('artifacts/scores').glob('*/*.json'))
    assert load(STAGE/'audit.json')['status'].startswith('PASS')
    rt.save(PLAN, dict(utc=datetime.now(timezone.utc).isoformat(),
        scope='Local-only isolated raw source extension. No source activation, recentering, prediction, metric or HOLDOUT response reads.',
        targets='Append exactly round49 72 targets in fixed panel order to original437; HMGXB4 stays below20 and excluded from additions. Destination79 unchanged.',
        moments='Original prepare_orion/correct_orion_controls: raw count sums; per-record value/full metadata cell total; sum squared per-record proportions; across-cell conservative target variance; first250 controls/sample before missing-block exclusion; cell-count weights for CPM, mapped-library weights for bulk; sum squared cell-count weights times sample mean variance. Publish float32.',
        identity='Preserve all original437 target rows and every non-row bank array byte-for-byte. Verify fresh global control moments against originals, then retain original bytes. New sample weights are separate artifacts. New namespace only; source centering/provenance revision required before use.',
        audit_tolerance='Independent sparse reduction: integer counts exact; float32 moments rtol2e-6 atol1e-7, probabilities rtol2e-6 atol1e-12. Old rows/non-row arrays exact dtype/shape/bytes. Full-depth and n_genes exact.',
        minimum_source_cells=20, controls_per_sample=250,
        sha256={str(p):rt.sha(p) for p in sorted(paths)}))
    print('FROZEN', PLAN, flush=True)


def check():
    plan = load(PLAN)
    for p,h in plan['sha256'].items():
        assert rt.sha(p) == h, p
    audit = load(STAGE/'audit.json')
    assert audit['receipt_sha256'] == rt.sha(STAGE/'receipt.json')
    assert audit['plan_sha256'] == rt.sha(STAGE/'plan.json')
    assert audit['verifier_sha256'] == rt.sha('audit_hct116_payload_round52.py')
    assert audit['selected_cells'] == 41769 and audit['logical_records'] == 215800779
    assert audit['independent_official_redecode_all_selected_records']
    for p,h in load(STAGE/'plan.json')['sha256'].items():
        assert rt.sha(p) == h, p
    index = load('artifacts/frozen_metadata_index.json')
    assert rt.sha(index['archive']) == index['archive_sha256']
    for p,d in index['entries'].items():
        assert rt.sha(p) == d['sha256'], p
    assert load('artifacts/calibration_gate.json')['status'] == 'CLOSED'
    return plan


def publish(path, **arrays):
    tmp = path.with_suffix('.tmp.npz')
    np.savez_compressed(tmp, **arrays)
    tmp.replace(path)


def main():
    check()
    assert not (ROOT/'receipt.json').exists(), 'Preserve completed evidence'
    assert not (ROOT/'bank.npz').exists(), 'Preserve partial publication for review'
    source = load(STAGE/'receipt.json')
    pre = load(PRE/'preflight.json')
    additions = np.array(pre['append_targets'])
    assert len(additions) == len(set(additions)) == 72
    with np.load(BASE) as d:
        old = {k:d[k] for k in d.files}
    assert len(old['targets']) == 437
    assert not set(additions) & set(old['targets'])
    assert 'HMGXB4' not in additions
    manifest = load('artifacts/tasks/manifest.json')
    assert not set(additions) & set(manifest['H1_holdout']['targets'])
    selected = pd.read_parquet(PRE/'selection.parquet')
    assert selected.cell_integer_id.is_unique and selected.pass_guide_filter.all()
    samples = np.array(sorted(selected['sample'].astype(str).unique()))
    nt, ns, ng = len(additions), len(samples), len(old['genes'])
    assert (nt,ns,ng) == (72,109,18536)
    ti = pd.Index(additions).get_indexer(selected.gene_target.astype(str))
    bi = pd.Index(samples).get_indexer(selected['sample'].astype(str))
    assert np.all((ti >= 0) | (selected.gene_target == 'Non-Targeting'))
    group = np.where(ti >= 0, ti, nt+bi)
    n = np.bincount(group, minlength=nt+ns)
    mapping = pd.Index(old['genes']).get_indexer(pd.read_parquet(GENES).gene_id.astype(str))
    measured = np.isin(np.arange(ng), mapping[mapping >= 0])
    np.testing.assert_array_equal(measured, old['measured_genes'])
    shape = (nt+ns,ng)
    sums, qp, q2 = (np.zeros(shape) for _ in range(3))
    dep = np.zeros(nt+ns)
    bn = np.zeros((nt,ns))
    bl = np.zeros_like(bn)
    valid = ti >= 0
    np.add.at(bn,(ti[valid],bi[valid]),1)
    ids = selected.cell_integer_id.to_numpy(np.int64)
    depths = selected.total_counts.to_numpy(float)
    counts_seen = np.zeros(len(selected),dtype=np.int64)
    depth_seen = np.zeros(len(selected))
    records = 0
    for i,(fid,item) in enumerate(sorted(source['decoded'].items(),key=lambda p:int(p[0]))):
        rt.guard(512*1024**2)
        path = STAGE/'decoded'/f'{fid}.npz'
        assert rt.sha(path) == item['sha256']
        with np.load(path) as d:
            cells, genes, values = d['cell_integer_id'],d['gene_integer_id'],d['value'].astype(float)
            local = np.searchsorted(ids,cells)
            assert np.all(local < len(ids)) and np.array_equal(ids[local],cells)
            assert len(cells) == len(genes) == len(values) == item['records']
            assert np.all(values > 0) and np.all(genes < len(mapping))
            np.add.at(counts_seen,local,1)
            np.add.at(depth_seen,local,values)
            col = mapping[genes]
            keep = col >= 0
            local, col, values = local[keep], col[keep], values[keep]
            g = group[local]
            q = values/depths[local]
            np.add.at(sums,(g,col),values)
            np.add.at(qp,(g,col),q)
            np.add.at(q2,(g,col),q*q)
            np.add.at(dep,g,values)
            ok = ti[local] >= 0
            np.add.at(bl,(ti[local][ok],bi[local][ok]),values[ok])
            records += len(cells)
        if i%100 == 0:
            print('moments',i+1,'/',len(source['decoded']),'records',records,flush=True)
    np.testing.assert_array_equal(counts_seen,selected.n_genes_by_counts)
    np.testing.assert_array_equal(depth_seen,depths)
    assert records == 215800779
    cn = n[nt:]
    control, control_q, control_q2 = sums[nt:], qp[nt:], q2[nt:]
    gp = control.sum(0)/control.sum()
    gm = control_q.sum(0)/cn.sum()
    gv = np.maximum(control_q2.sum(0)/cn.sum()-gm*gm,0)/max(cn.sum()-1,1)
    cp = np.divide(control,dep[nt:,None],out=np.tile(gp,(ns,1)),where=dep[nt:,None]>0)
    cm = np.divide(control_q,cn[:,None],out=np.tile(gm,(ns,1)),where=cn[:,None]>0)
    cv = np.maximum(np.divide(control_q2,cn[:,None],out=np.zeros_like(control_q2),where=cn[:,None]>0)-cm*cm,0)/np.maximum(cn[:,None]-1,1)
    cv[cn == 0] = gv
    tn = n[:nt].astype(np.int64)
    assert np.all(tn >= 20)
    tm = qp[:nt]/tn[:,None]*1e6
    pv = np.maximum(q2[:nt]/tn[:,None]-(tm/1e6)**2,0)/np.maximum(tn[:,None]-1,1)*1e12
    cw = bn/tn[:,None]
    lw = bl/dep[:nt,None]
    mp, mm, vc = lw@cp, cw@cm*1e6, (cw*cw)@cv*1e12
    fresh = dict(n_cells=tn,target_count_sums=sums[:nt].astype('f4'),
        target_mean_cpm=tm.astype('f4'),matched_control_probability=mp.astype('f4'),
        matched_control_mean_cpm=mm.astype('f4'),delta_cpm_variance=(pv+vc).astype('f4'),
        target_mean_cpm_variance=pv.astype('f4'),matched_control_mean_cpm_variance=vc.astype('f4'))
    for k,v in fresh.items():
        assert np.isfinite(v).all() and np.all(v >= 0),k
    np.testing.assert_allclose(gp.astype('f4'),old['global_control_probability'],rtol=2e-6,atol=1e-12)
    np.testing.assert_allclose((gm*1e6).astype('f4'),old['global_control_mean_cpm'],rtol=2e-6,atol=1e-7)
    union = dict(old)
    union['targets'] = np.concatenate([old['targets'],additions])
    assert len(union['targets']) == 509
    for k in ROW_KEYS:
        union[k] = np.concatenate([old[k],fresh[k]],axis=0)
    publish(ROOT/'bank.npz',**union)
    publish(ROOT/'extension_blocks.npz',targets=additions,samples=samples,
        target_batch_n=bn,target_batch_library=bl,control_n=cn,
        control_count_sums=control,control_sum_proportion=control_q,
        control_sum_squared_proportion=control_q2,
        fresh_global_control_probability=gp.astype('f4'),fresh_global_control_mean_cpm=(gm*1e6).astype('f4'))
    check()
    rt.save(ROOT/'receipt.json',dict(status='BUILT isolated raw bank; independent moment audit pending',
        utc=datetime.now(timezone.utc).isoformat(),plan_sha256=rt.sha(PLAN),
        bank_sha256=rt.sha(ROOT/'bank.npz'),blocks_sha256=rt.sha(ROOT/'extension_blocks.npz'),
        original_targets=437,appended_targets=72,union_targets=509,selected_cells=len(selected),
        control_cells=int(cn.sum()),target_cells=int(tn.sum()),logical_records=records,
        source_activated=False,centering_revised=False,holdout_scored=False,phase2_allowed=False))
    print('BUILT isolated509-target bank; independent audit pending',flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--freeze',action='store_true')
    args=parser.parse_args()
    freeze() if args.freeze else main()
