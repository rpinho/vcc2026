"""AMENDMENT_R9 instrument: stratum-reweighted official-scale projection vs B0, seed floor, validity. Writes results/r9_gate.{json,md}."""
import json,csv,sys,itertools,numpy as np
from pathlib import Path
R=Path(__file__).resolve().parent;S=R/'results'
L=lambda n:json.loads((S/n).read_text())
st={r['target']:r['stratum'] for r in csv.DictReader(open(R.parent/'layout_cd4_1025/artifacts/per_target.csv')) if r['recipient']=='H1'}
SN=['<10','10-999','>=1000'];W0=(0.315,0.65,0.035)
C=dict(fid='de_wilcoxon_direction_fidelity_yield_raw',jac='de_wilcoxon_sig_jaccard',reach='de_wilcoxon_direction_reach_raw',nmae='de_wilcoxon_lfc_nmae')
b=L('official_B0.json')
# local scaled per raw unit for reach and nmae, from B0's own scaling (raw->scaled is affine per component)
def slope(k):
  a=b['anchors'] if isinstance(b.get('anchors'),dict) else None
  return None
def tm(a):
  z=np.load(S/f'target_metrics_{a}.npz',allow_pickle=True);cols=list(z['columns']);return {t:dict(zip(cols,v)) for t,v in zip(z['targets'],z['values'])}
TB=tm('B0');T=sorted(TB)
# local affine slopes for reach / nmae: scaled delta per raw delta, taken from the D1770 official pair (exact affine scaling)
d=L('official_D1770.json')
ls={k:(d['scaled'][C[k]]-b['scaled'][C[k]])/(d['raw'][C[k]]-b['raw'][C[k]]) for k in ['reach','nmae'] if d['raw'][C[k]]!=b['raw'][C[k]]}
OFF=dict(fid=3.30,jac=2.88,reach=ls['reach'],nmae=ls.get('nmae',0.0))
def project(a):
  TA=tm(a);out={};per={}
  for k,col in C.items():
    dd={t:float(TA[t][col]-TB[t][col]) for t in T}
    # nmae is omitted by the scorer for targets with <10 gated genes; an all-omitted stratum contributes 0.
    per[k]={s:(lambda v:float(np.nanmean(v)) if np.isfinite(v).any() else 0.0)(np.array([dd[t] for t in T if st[t]==s])) for s in SN}
  def comp(w):return sum(OFF[k]*sum(w[i]*per[k][s] for i,s in enumerate(SN)) for k in C)/6
  sweep=[]
  for ps,pw in itertools.product([0.02,0.035,0.05,0.08,0.15],[0.20,0.315,0.45]):sweep.append(comp((pw,1-ps-pw,ps)))
  wm=[t for t in T if st[t]!='>=1000'];fd={t:float(TA[t][C['fid']]-TB[t][C['fid']]) for t in wm}
  return dict(per_stratum_raw=per,central=comp(W0),conservative=min(sweep),optimistic=max(sweep),
    weakmid_fid_improved_share=float(np.mean([v>0 for v in fd.values()])),weakmid_fid_worse_share=float(np.mean([v<0 for v in fd.values()])))
seeds={a:project(a)['central'] for a in ['B0s1','B0s2']};floor=max(abs(v) for v in seeds.values())
res={'slopes':OFF,'strata_counts':{s:sum(st[t]==s for t in T) for s in SN},'seed_floor':floor,'seed_central':seeds,'arms':{}}
for a in [x for x in ['RH20','K250','DET885'] if (S/f'official_{x}.json').exists()]:
  o=L(f'official_{a}.json');p=project(a)
  with np.load(S/f'rendered_bulk_{a}.npz') as x,np.load(S/'rendered_bulk_B0.npz') as y:bulk=np.array_equal(x['count_sums'],y['count_sums'])
  K='expr_mse_unbiased_capped_norm';mr=o['raw'][K]/b['raw'][K];nm=o['raw'][C['nmae']]-b['raw'][C['nmae']]
  valid=bulk and o['raw']['pds_cosine']==b['raw']['pds_cosine'] and mr<=1.001 and abs(nm)<=0.001
  pb=L(f'paired_boot_{a}_B0.json')['comparison']['score_avg'] if (S/f'paired_boot_{a}_B0.json').exists() else {}
  rd=L(f'render_{a}.json');calls=[len(r['selected']) for r in rd['records']]
  screen=valid and p['central']>=0.005 and p['central']>2*floor
  res['arms'][a]=dict(local_composite=o['score_avg']-b['score_avg'],local_ci=pb.get('ci95'),projection=p,valid=valid,bulk_identical=bulk,
    pds_identical=o['raw']['pds_cosine']==b['raw']['pds_cosine'],raw_mse_ratio=mr,nmae_raw_delta=nm,layout_calls_median=float(np.median(calls)),
    screen='PASS' if screen else ('VOID' if not valid else 'FAIL'))
(S/'r9_gate.json').write_text(json.dumps(res,indent=1))
md=['# R9 gate (AMENDMENT_R9.md)',f"strata {res['strata_counts']}; slopes {OFF}; seed floor {floor:+.5f} (B0s1 {seeds['B0s1']:+.5f}, B0s2 {seeds['B0s2']:+.5f})",'',
 '| arm | local d | central proj | conservative | optimistic | wm fid better/worse | valid | MSE ratio | nmae d | screen |','|---|---:|---:|---:|---:|---|---|---:|---:|---|']
for a,v in res['arms'].items():
  p=v['projection'];md.append(f"| {a} | {v['local_composite']:+.5f} | {p['central']:+.5f} | {p['conservative']:+.5f} | {p['optimistic']:+.5f} | {p['weakmid_fid_improved_share']:.2f}/{p['weakmid_fid_worse_share']:.2f} | {v['valid']} | {v['raw_mse_ratio']:.5f} | {v['nmae_raw_delta']:+.5f} | **{v['screen']}** |")
md.append('');md+= [f"- {a} per-stratum raw deltas: "+'; '.join(f"{k} "+'/'.join(f"{v['projection']['per_stratum_raw'][k][s]:+.4f}" for s in SN) for k in C) for a,v in res['arms'].items()]
(S/'r9_gate.md').write_text('\n'.join(md)+'\n');print('\n'.join(md))
