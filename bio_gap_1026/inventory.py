"""Inventory only. Existing metadata; no network/data downloads."""
import runtime as rt
rt.enforce()
import pandas as pd,numpy as np
targets=rt.json.loads((rt.ROOT/'results/training_metadata.json').read_text())['targets'];T=set(targets)
md=rt.json.loads((rt.ROOT/'results/training_metadata.json').read_text());covered={t for t,r in zip(targets,md['coverage']) if r};missing=T-covered
source=rt.REPO/'research/sources/forrestsheldon__virtual-cell-challenge-2026/metadata/datasets/source_targets'
rows=[];pins={}
errors=pd.read_csv(rt.ROOT/'results/per_target.csv').set_index('target').sse
def add(name,path,family,url,scope='complete',query=None,quality=None):
    pins[str(path)]=rt.sha(path);df=pd.read_csv(path,sep=r'\s+' if path.suffix=='.tsv' else ',')
    if query is not None:df=query(df)
    gene=next(c for c in ['source_target_gene','gene','target','pert_iname'] if c in df)
    obs=set(df[gene].astype(str));eligible=obs if quality is None else set(quality(df)[gene].astype(str))
    C=T&eligible;new=C-covered;existing_family={t for t,r in zip(targets,md['coverage']) if any(v['family']==family for v in r)}
    st=dict(source=name,family=family,training168_covered=len(C),training168_fraction=len(C)/168,new_to_union=len(new),new_to_union_fraction=len(new)/168,new_to_this_family=len(C-existing_family),new_context_fraction=len(C-existing_family)/168,missing_row_perfect_replacement_SSE_ceiling_pct=float(100*errors.reindex(sorted(new)).sum()/errors.sum()),covered_targets=sorted(C),new_targets=sorted(new),scope=scope,metadata=str(path),primary_url=url)
    rows.append(st)
add('KOLF pan-genome',source/'nourreddine2024_kolf_counts.csv','iPSC','https://doi.org/10.25452/figshare.plus.27261219',quality=lambda d:d[d.n_cells>=30])
add('HipSci GW',rt.REPO/'data/expansion/references/Feng_GW_target_line_counts.csv.gz','iPSC','https://doi.org/10.1016/j.xgen.2025.101076',scope='LOWER BOUND ONLY: retained table was restricted to the 2026 validation panel; absent training genes are unknown',query=lambda d:d.groupby('target',as_index=False)['0'].sum(),quality=lambda d:d[d['0']>=30])
add('VIPerturb genome-wide',rt.REPO/'data/expansion/references/vip_manifest.tsv','K562','https://zenodo.org/records/18460279',scope='Full manifest; >=30 cells, fingerprint-filter status ignored',quality=lambda d:d[d.cell_count>=30])
add('VIPerturb fingerprint-filtered',rt.REPO/'data/expansion/references/vip_manifest.tsv','K562','https://zenodo.org/records/18460279',scope='Full manifest; author fingerprint filter plus >=30 cells',quality=lambda d:d[(d.cell_count>=30)&d.is_in_filtered_obj.astype(str).str.upper().eq('TRUE')])
add('LINCS genetic KD/KO',rt.REPO/'data/expansion/lincs_target_cellline_counts.csv.gz','LINCS','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742',scope='LOWER BOUND ONLY: signature table restricted to 2026 validation targets; bulk L1000, not full scRNA CRISPRi',query=lambda d:d[d.pert_type.isin(['trt_sh','trt_xpr'])])
add('Marson/Pritchard genome-wide CD4',source/'zhu2026_cd4t_counts.csv','CD4','https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq',scope='Full metadata, cells pooled across 4 donors and 3 states; >=30',quality=lambda d:d[d.n_cells>=30])
for name,file,fam in [('Replogle K562 essential','replogle_k562_essential_counts.csv','K562'),('Replogle K562 genome-wide','replogle_k562_gwps_counts.csv','K562'),('Replogle RPE1 essential','replogle_rpe1_counts.csv','RPE1'),('Nadig Jurkat essential','nadig_jurkat_counts.csv','Jurkat'),('Nadig HepG2 essential','nadig_hepg2_counts.csv','HepG2')]:
    add(name,source/file,fam,'https://zenodo.org/records/13350497',scope='Complete source metadata; >=30 cells',quality=lambda d:d[d.n_cells>=30])
pd.DataFrame([{k:v for k,v in a.items() if k not in ['covered_targets','new_targets']} for a in rows]).to_csv(rt.ROOT/'results/source_inventory.csv',index=False)
rt.save('results/source_inventory.json',dict(training_count=168,covered_current=len(covered),uncovered_targets=sorted(missing),rows=rows,pins=pins,qualification='Metadata inventory; coverage is not efficacy, observed readout compatibility, or measured transfer gain. CD4/Replogle/Nadig are already present in B0, but retained per-family target sets can be incomplete.'))
for r in rows:print(r['source'],r['training168_covered'],r['new_to_union'],r['new_to_this_family'],r['new_targets'])
rt.status('Missing-family inventory completed using existing local metadata only; HipSci/LINCS fractions are explicit lower bounds.')
