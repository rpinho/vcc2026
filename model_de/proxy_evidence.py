"""Read the feedback job's existing official-code measurements, never modify them."""
import pathlib,json
P=pathlib.Path(__file__).resolve().parent
names={'pds':'pds_cosine','mse':'expr_mse_unbiased_capped_norm','nmae':'de_wilcoxon_lfc_nmae',
       'fid':'de_wilcoxon_direction_fidelity_yield_raw','reach':'de_wilcoxon_direction_reach_raw','jac':'de_wilcoxon_sig_jaccard'}
out={}
for dest in ['H1','HEK293T','K562','CD4','HCT116']:
    rows={}
    for arm in ['v4a','v4b','v5a','v5b','v5b_pie']:
        path=P.parent/f'proxy/artifacts/scores/{dest}/{arm}.json'
        if path.exists():
            data=json.loads(path.read_text());rows[arm]={k:data['scaled'].get(v) for k,v in names.items()}
    if rows:out[dest]=rows
(P/'artifacts/proxy_existing.json').write_text(json.dumps(out,indent=2))
for dest,rows in out.items():
    for arm,values in rows.items():print(dest,arm,' '.join(f'{k}={v:.4f}' for k,v in values.items() if v is not None))
