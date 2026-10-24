import common as rt
import io_guard; io_guard.install()
import re, hashlib, numpy as np, pandas as pd, anndata as ad

def main():
    dest=rt.ROOT/'artifacts/expansion_targets.json'
    if dest.exists():
        r=rt.json.loads(dest.read_text())
        for file,digest in r['files'].items(): assert rt.sha(rt.ROOT/file)==digest
        print('Frozen E unchanged',r['remote_targets']); return
    path=rt.REPO/'public_data/K562_GWPS_bulk.h5ad'
    a=ad.read_h5ad(path,backed='r')
    lab=np.array(['non-targeting' if 'non-targeting' in s else re.match(r'^[^_]+_(.+)_[^_]+_([^_]+)$',s).group(1) for s in a.obs_names])
    eligible=sorted(set(lab[a.obs.num_cells_filtered.fillna(0).to_numpy()>=30])-{'non-targeting'})
    a.file.close()
    essential=set()
    for name in ['K562','Jurkat','HepG2','RPE1']:
        with np.load(rt.REPO/f'model_billm/data/{name}.npz') as z: essential.update(z['targets'].astype(str))
    mock=pd.read_csv(rt.REPO/'model_v6/mock_final/panel.csv').target_gene.astype(str).tolist()
    order=sorted(set(eligible)-essential,key=lambda t:hashlib.sha256(('vcc_billmB_1021:'+t).encode()).hexdigest())
    selected=set(mock)
    for t in order:
        if len(selected)>=1500:break
        selected.add(t)
    remote=sorted(selected,key=lambda t:hashlib.sha256(('vcc_billmB_1021:'+t).encode()).hexdigest())
    for file,targets in [('targets_local.csv',eligible),('targets_remote.csv',remote)]:
        pd.DataFrame({'target_gene':targets}).to_csv(rt.ROOT/file,index=False)
    rt.save(dest,dict(rule='SHA256(vcc_billmB_1021:SYMBOL) order, non-essential eligible GWPS; mandatory union mock panel; total 1500',
        response_blind=True,local_targets=len(eligible),remote_targets=len(remote),mandatory_mock=len(mock),
        mock_sha256=rt.sha(rt.REPO/'model_v6/mock_final/panel.csv'),bulk_sha256=rt.sha(path),
        essential_definition='union of >=30-cell targets in four original billmA essential summaries',
        essential_summary_hashes={n:rt.sha(rt.REPO/f'model_billm/data/{n}.npz') for n in ['K562','Jurkat','HepG2','RPE1']},
        files={f:rt.sha(rt.ROOT/f) for f in ['targets_local.csv','targets_remote.csv']},
        budget='1500 initial; unchanged unless throughput-only SHA-prefix extension is frozen before training; <=40GB total planned, ~8h remote wall',
        acceptance='strictly exceed frozen v5b six-metric mean separately on H1_rich_search and CD4_depth192_search; sealed HOLDOUT forbidden'))
    print('Frozen local',len(eligible),'remote',len(remote),'mock',len(mock),flush=True)
if __name__=='__main__': main()
