"""Independent grouped-bincount replay of isolated HCT116 source moments."""
from datetime import datetime, timezone
from pathlib import Path

import runtime as rt
import numpy as np
import pandas as pd

ROOT = Path('artifacts/HCT116_extension_round53')
PRE = Path('artifacts/HCT116_extension_preflight_round49')
STAGE = Path('artifacts/HCT116_payload_round52')


def load(p):
    return rt.json.loads(Path(p).read_text())


def bindings():
    plan = load(ROOT/'plan.json')
    for p,h in plan['sha256'].items():
        assert rt.sha(p) == h,p
    for p,h in load(STAGE/'plan.json')['sha256'].items():
        assert rt.sha(p) == h,p
    audit = load(STAGE/'audit.json')
    assert audit['status'].startswith('PASS')
    assert audit['receipt_sha256'] == rt.sha(STAGE/'receipt.json')
    assert audit['plan_sha256'] == rt.sha(STAGE/'plan.json')
    assert audit['verifier_sha256'] == rt.sha('audit_hct116_payload_round52.py')
    assert audit['selected_cell_ids_counts_and_full_depths_verified']
    assert audit['independent_official_redecode_all_selected_records']
    index = load('artifacts/frozen_metadata_index.json')
    assert rt.sha(index['archive']) == index['archive_sha256']
    for p,d in index['entries'].items():
        assert rt.sha(p) == d['sha256'],p
    return plan


def main():
    assert not (ROOT/'audit.json').exists()
    bindings()
    receipt = load(ROOT/'receipt.json')
    assert receipt['plan_sha256'] == rt.sha(ROOT/'plan.json')
    assert receipt['bank_sha256'] == rt.sha(ROOT/'bank.npz')
    assert receipt['blocks_sha256'] == rt.sha(ROOT/'extension_blocks.npz')
    with np.load('artifacts/cd4_banks/HCT116.npz') as d:
        old = {k:d[k] for k in d.files}
    with np.load(ROOT/'bank.npz') as d:
        bank = {k:d[k] for k in d.files}
    with np.load(ROOT/'extension_blocks.npz') as d:
        blocks = {k:d[k] for k in d.files}
    assert set(bank) == set(old)
    row_keys = {'n_cells','target_count_sums','target_mean_cpm','matched_control_probability',
                'matched_control_mean_cpm','delta_cpm_variance','target_mean_cpm_variance',
                'matched_control_mean_cpm_variance'}
    for key,original in old.items():
        observed = bank[key][:437] if key in row_keys or key == 'targets' else bank[key]
        assert observed.dtype == original.dtype and observed.shape == original.shape,key
        assert observed.tobytes() == original.tobytes(),key
    target_table = pd.read_csv(PRE/'targets.tsv',sep='\t')
    targets = target_table.loc[target_table.action == 'append','target'].to_numpy(str)
    np.testing.assert_array_equal(bank['targets'][437:],targets)
    np.testing.assert_array_equal(blocks['targets'],targets)
    assert len(targets) == 72 and len(bank['targets']) == 509
    assert not set(targets) & set(old['targets'])
    abstain = target_table.loc[target_table.action == 'abstain_below20']
    assert abstain.target.tolist() == ['HMGXB4'] and abstain.block_present_cells.tolist() == [2]
    manifest = load('artifacts/tasks/manifest.json')
    assert len(manifest['CD4_depth192_search']['targets']) == 79
    assert 'HMGXB4' in manifest['CD4_depth192_search']['targets']
    assert not set(targets) & set(manifest['H1_holdout']['targets'])
    selected = pd.read_parquet(PRE/'selection.parquet')
    assert selected.cell_integer_id.is_unique and selected.pass_guide_filter.all()
    selected = selected.sort_values('cell_integer_id').reset_index(drop=True)
    samples = sorted(selected['sample'].astype(str).unique())
    np.testing.assert_array_equal(blocks['samples'],samples)
    tmap = {t:i for i,t in enumerate(targets)}
    smap = {s:i for i,s in enumerate(samples)}
    nt,ns,ng = 72,len(samples),len(old['genes'])
    assert (ns,ng) == (109,18536)
    group = np.array([tmap[r.gene_target] if r.gene_target != 'Non-Targeting' else nt+smap[r.sample]
                      for r in selected.itertuples(index=False)])
    n = np.bincount(group,minlength=nt+ns)
    depth = selected.total_counts.to_numpy(float)
    ids = selected.cell_integer_id.to_numpy(np.int64)
    mapping_dict = {g:i for i,g in enumerate(old['genes'])}
    raw_genes = pd.read_parquet('artifacts/HCT116_expression_identity_round50/genes.parquet').gene_id.astype(str)
    gene_map = np.array([mapping_dict.get(g,-1) for g in raw_genes])
    shape = (nt+ns,ng)
    raw = np.zeros(shape)
    prop = np.zeros(shape)
    prop2 = np.zeros(shape)
    libraries = np.zeros((nt,ns))
    batch_n = np.zeros((nt,ns))
    for r in selected.itertuples(index=False):
        if r.gene_target != 'Non-Targeting':
            batch_n[tmap[r.gene_target],smap[r.sample]] += 1
    sample_codes = np.array([smap[s] for s in selected['sample']])
    counts_seen = np.zeros(len(selected),dtype=np.int64)
    depths_seen = np.zeros(len(selected))
    records = 0
    source = load(STAGE/'receipt.json')
    for i,(fid,item) in enumerate(sorted(source['decoded'].items(),key=lambda p:int(p[0]))):
        rt.guard(512*1024**2)
        path = STAGE/'decoded'/f'{fid}.npz'
        assert rt.sha(path) == item['sha256']
        with np.load(path) as d:
            cells, genes, values = d['cell_integer_id'],d['gene_integer_id'],d['value'].astype(float)
            assert len(cells) == len(genes) == len(values) == item['records']
            at = np.searchsorted(ids,cells)
            assert np.all(at < len(ids)) and np.array_equal(ids[at],cells)
            counts_seen += np.bincount(at,minlength=len(ids))
            depths_seen += np.bincount(at,weights=values,minlength=len(ids))
            cols = gene_map[genes]
            keep = cols >= 0
            at,cols,values = at[keep],cols[keep],values[keep]
            groups = group[at]
            flat = groups*ng+cols
            q = values/depth[at]
            for out,weights in [(raw,values),(prop,q),(prop2,q*q)]:
                out += np.bincount(flat,weights=weights,minlength=(nt+ns)*ng).reshape(shape)
            ok = groups < nt
            libraries += np.bincount(groups[ok]*ns+sample_codes[at[ok]],weights=values[ok],
                                    minlength=nt*ns).reshape(nt,ns)
            records += len(cells)
        if i%200 == 0:
            print('independent moments',i+1,'records',records,flush=True)
    assert records == 215800779 and len(source['decoded']) == 2018
    np.testing.assert_array_equal(counts_seen,selected.n_genes_by_counts)
    np.testing.assert_array_equal(depths_seen,depth)
    np.testing.assert_array_equal(n[:nt],target_table.loc[target_table.action=='append','block_present_cells'])
    np.testing.assert_array_equal(bank['n_cells'][437:],n[:nt])
    assert np.all(n[:nt] >= 20) and int(n[:nt].sum()) == 14601 and int(n[nt:].sum()) == 27168
    np.testing.assert_array_equal(blocks['target_batch_n'],batch_n)
    np.testing.assert_array_equal(blocks['target_batch_library'],libraries)
    np.testing.assert_array_equal(blocks['control_n'],n[nt:])
    np.testing.assert_array_equal(blocks['control_count_sums'],raw[nt:])
    np.testing.assert_allclose(blocks['control_sum_proportion'],prop[nt:],rtol=1e-12,atol=1e-13)
    np.testing.assert_allclose(blocks['control_sum_squared_proportion'],prop2[nt:],rtol=1e-12,atol=1e-15)
    np.testing.assert_array_equal(bank['target_count_sums'][437:],raw[:nt].astype('f4'))
    control_n = n[nt:]
    assert np.all(control_n > 1)
    control_mean = prop[nt:]/control_n[:,None]
    control_probability = raw[nt:]/raw[nt:].sum(axis=1)[:,None]
    control_mean_var = np.maximum(prop2[nt:]/control_n[:,None]-control_mean**2,0)/(control_n[:,None]-1)
    means = prop[:nt]/n[:nt,None]*1e6
    target_var = np.maximum(prop2[:nt]/n[:nt,None]-(means/1e6)**2,0)/(n[:nt,None]-1)*1e12
    mp,mm,vc = (np.zeros((nt,ng)) for _ in range(3))
    # Scalar sample-weight traversal is independent of producer BLAS products.
    for t in range(nt):
        for s in range(ns):
            weight = batch_n[t,s]/n[t]
            lib_weight = libraries[t,s]/libraries[t].sum()
            mp[t] += lib_weight*control_probability[s]
            mm[t] += weight*control_mean[s]*1e6
            vc[t] += weight**2*control_mean_var[s]*1e12
    expected = dict(target_mean_cpm=means,matched_control_probability=mp,
        matched_control_mean_cpm=mm,target_mean_cpm_variance=target_var,
        matched_control_mean_cpm_variance=vc,delta_cpm_variance=target_var+vc)
    errors = {}
    for key,value in expected.items():
        observed = bank[key][437:]
        assert observed.dtype == np.dtype('float32') and np.isfinite(observed).all()
        assert np.all(observed >= 0)
        atol = 1e-12 if key == 'matched_control_probability' else 1e-7
        np.testing.assert_allclose(observed,value.astype('f4'),rtol=2e-6,atol=atol,err_msg=key)
        errors[key] = float(np.max(np.abs(observed.astype(float)-value)))
    globals_ = {'global_control_probability':raw[nt:].sum(axis=0)/raw[nt:].sum(),
                'global_control_mean_cpm':prop[nt:].sum(axis=0)/control_n.sum()*1e6}
    for key,value in globals_.items():
        atol = 1e-12 if key.endswith('probability') else 1e-7
        np.testing.assert_allclose(bank[key],value.astype('f4'),rtol=2e-6,atol=atol)
        np.testing.assert_allclose(blocks['fresh_'+key],value.astype('f4'),rtol=2e-6,atol=atol)
    # Check the frozen sample-count/control summaries by a separate table path.
    desired = pd.read_csv(PRE/'target_sample_counts.tsv',sep='\t',index_col=0).reindex(index=targets,columns=samples)
    np.testing.assert_array_equal(batch_n,desired.to_numpy())
    controls = pd.read_csv(PRE/'controls.tsv',sep='\t',index_col=0).reindex(samples)
    np.testing.assert_array_equal(control_n,controls.cells)
    bindings()
    assert load('artifacts/calibration_gate.json')['status'] == 'CLOSED'
    rt.save(ROOT/'audit.json',dict(status='PASS independent raw moment/source isolation audit; NOT ACTIVATED',
        utc=datetime.now(timezone.utc).isoformat(),plan_sha256=rt.sha(ROOT/'plan.json'),
        receipt_sha256=rt.sha(ROOT/'receipt.json'),auditor_sha256=rt.sha(__file__),
        bank_sha256=rt.sha(ROOT/'bank.npz'),blocks_sha256=rt.sha(ROOT/'extension_blocks.npz'),
        original437_rows_and_nonrow_arrays_byte_identical=True,
        appended_targets=72,union_targets=509,logical_records=records,selected_cells=len(selected),
        full_depths_and_counts_exact=True,sample_count_and_mapped_library_weights_exact=True,
        all_moments_independently_reduced=True,maximum_absolute_float64_reference_errors=errors,
        frozen_metadata_members=38,source_activated=False,centering_revised=False,
        holdout_scored=False,phase2_allowed=False))
    print('PASS independent509-target raw moment audit; NOT ACTIVATED',flush=True)


if __name__ == '__main__':
    main()
