"""GENCODE promoter pairs and independent public CD4 guide-neighbour support.
Guide coordinates are not assumed to be the challenge's undisclosed guides.
"""
import runtime,gzip,re,requests
import pandas as pd,numpy as np
from pathlib import Path
rt=runtime.ROOT;p=rt/'artifacts/CD4_guide_metadata.csv'
if not p.exists():
 s=requests.Session();s.trust_env=False
 with s.get('https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/suppl_tables/sgrna_library_metadata.suppl_table.csv',stream=True,timeout=60) as r:
  r.raise_for_status();value=r.raw.read(15_000_001);assert len(value)<15_000_001;p.write_bytes(value)
guides=pd.read_csv(p);print('guide columns',guides.columns.tolist(),flush=True)
rows=[]
with gzip.open(rt.parent/'public_data/gencode.v47.basic.annotation.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'):continue
  z=line.rstrip().split('\t')
  if z[2]!='gene':continue
  attrs=dict(re.findall(r'(\w+) "([^"]+)"',z[8]));name=attrs.get('gene_name')
  if name:rows.append(dict(gene=name,chromosome=z[0],tss=int(z[3] if z[6]=='+' else z[4]),strand=z[6]))
a=pd.DataFrame(rows);a=a[~a.gene.duplicated(keep=False)].set_index('gene');targets=runtime.json.loads((rt/'artifacts/targets.json').read_text());records=[]
for target in targets:
 if target not in a.index:continue
 t=a.loc[target];near=a[(a.chromosome==t.chromosome)&((a.tss-t.tss).abs()<=5000)]
 for neighbor,b in near.iterrows():
  if neighbor==target:continue
  divergent=t.strand!=b.strand and ((t.strand=='+' and t.tss>b.tss) or (t.strand=='-' and t.tss<b.tss))
  # Use publisher gene-labelled neighbour assignments rather than mixing genome builds.
  supported=guides[(guides.target_gene_name==target)&(guides.nearest_within2kb_nontarget_gene_name==neighbor)]
  distances=pd.to_numeric(supported.nearest_within2kb_nontarget_gene_dist,errors='coerce').abs().dropna();n_guides=supported.sgRNA.nunique() if 'sgRNA' in supported else len(supported)
  records.append(dict(target=target,neighbor=neighbor,distance=abs(int(t.tss)-int(b.tss)),divergent=bool(divergent),guide_supported=bool(n_guides>=2 and len(distances)>=2),guide_distance=float(distances.median()) if len(distances) else float('nan'),guide_count=int(n_guides),guide_provenance='public CD4 source library; challenge guides unavailable'))
pairs=pd.DataFrame(records);pairs.to_csv(rt/'artifacts/promoter_pairs.csv',index=False)
runtime.save(rt/'artifacts/promoter_summary.json',dict(pairs=len(pairs),targets=int(pairs.target.nunique()),guide_supported_pairs=int(pairs.guide_supported.sum()),source_guide_bytes=p.stat().st_size,limitation='Challenge guide identities/positions are unavailable. This is a source-guide-supported prior; do not claim exact destination cis effects.'))
print(pairs[pairs.guide_supported].head(10).to_string(index=False),flush=True)
