"""Audit actual source coverage and both model routes' destination exclusions."""
import runtime as rt
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
import recipes
folder=Path('artifacts/tasks');manifest=rt.json.loads((folder/'manifest.json').read_text());m=manifest['CD4_matched_search'];m['selection_policy']='300 fixed-hash perturbations among 334 eligible (non-targeting excluded); minimum 200 raw GWPS source cells; count-only reference quota fallback';m['budget_fallback']='Before extraction, reduce target quota to 4, input to 80 and reference to 32 per lane, preserving target identities';m['replicate_splits']=3;manifest['CD4_matched_search']=m;rt.save(folder/'manifest.json',manifest)
d=np.load('artifacts/cd4_banks/K562_GWPS_raw.npz');lookup=pd.Index(d['targets']).get_indexer(m['targets']);assert (lookup>=0).all() and (d['n_cells'][lookup]>=200).all();coverage={}
for file in Path('artifacts/cd4_sources').glob('*.npz'):
 z=np.load(file);coverage[file.stem]=int(np.isin(m['targets'],z['targets']).sum())
v4=recipes.v4;before=(v4.BANKS,v4.load_banks,v4.cd4_family,v4.atlas.aligned_cd4);heads=v4.heads;emit=v4.emit;cases=[]
def spy(targets,genes,mean,bulk,variant,exclude=None,stats=None):
 assert exclude=='CD4' and v4.BANKS==rt.ROOT/'artifacts/cd4_banks';assert v4.cd4_family(targets,genes,'log2fc') is None;effect,available=v4.atlas.aligned_cd4('MUST_NOT_OPEN',targets,genes);assert not effect.any() and not available.any();assert set(v4.load_banks(exclude='CD4'))=={'K562_GWPS_raw','HCT116'};cases.append(variant);return [np.tile(mean,(len(targets),1)),np.tile(bulk,(len(targets),1))],[np.zeros((len(targets),len(genes)),bool)]*2
v4.heads=spy;v4.emit=lambda q,p,b,depth,seed,variant:np.ones((400,len(p)),np.int32)
ctrl=ad.AnnData(sparse.csr_matrix(np.ones((1600,4),np.int32)),obs=pd.DataFrame({'target_gene':['non-targeting']*1600},index=[str(i) for i in range(1600)]),var=pd.DataFrame(index=['a','b','c','d']))
try:
 for recipe in ['v4a','v4b']:
  list(recipes.predict(ctrl,['x'],'CD4',recipe,42));assert before==(v4.BANKS,v4.load_banks,v4.cd4_family,v4.atlas.aligned_cd4)
finally:v4.heads=heads;v4.emit=emit
rt.save('artifacts/CD4_source_audit.json',dict(status='PASS matched-context source coverage and both destination-family exclusion routes',v4_raw_covered=int((d['n_cells'][lookup]>=200).sum()),legacy_coverage=coverage,route_tests=cases,static_promoter_annotations='Public DNA annotations; not response-bank evidence',HOLDOUT='No scores accessed'));print('PASS CD4 source audit',coverage,cases)
