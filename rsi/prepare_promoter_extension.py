"""Extend the unchanged promoter recipe to new genes; preserve every old row."""
import runtime as rt
from pathlib import Path
import pandas as pd,numpy as np,gzip,re,requests,shutil
root=Path('artifacts');original=root/'promoter_pairs_original.csv'
if not original.exists():shutil.copyfile(root/'promoter_pairs.csv',original)
old=pd.read_csv(original);canonical=set(rt.json.loads((rt.ROOT.parent/'model_v4/artifacts/targets.json').read_text()));wanted=set(rt.json.loads((root/'CD4_wide_targets.json').read_text()))-canonical;manifest=rt.json.loads((root/'tasks/manifest.json').read_text());h1=set(t for m in manifest.values() if m['context']=='H1' for t in m['targets']);assert not(wanted&h1),'New annotations would alter frozen H1 targets'
file=root/'CD4_guide_metadata.csv';downloaded=0
if not file.exists():
 s=requests.Session();s.trust_env=False
 with s.get('https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/suppl_tables/sgrna_library_metadata.suppl_table.csv',stream=True,timeout=60) as response:
  response.raise_for_status();data=response.raw.read(15_000_001);downloaded=len(data);rt.save(root/'promoter_extension_transfer.json',dict(bytes=downloaded,cap=15_000_001));assert downloaded<=15_000_000;file.write_bytes(data)
guides=pd.read_csv(file);rows=[]
with gzip.open(rt.ROOT.parent/'public_data/gencode.v47.basic.annotation.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'):continue
  z=line.rstrip().split('\t')
  if z[2]!='gene':continue
  attrs=dict(re.findall(r'(\w+) "([^"]+)"',z[8]));name=attrs.get('gene_name')
  if name:rows.append(dict(gene=name,chromosome=z[0],tss=int(z[3] if z[6]=='+' else z[4]),strand=z[6]))
a=pd.DataFrame(rows);a=a[~a.gene.duplicated(keep=False)].set_index('gene');records=[]
for target in sorted(wanted):
 if target not in a.index:continue
 t=a.loc[target];near=a[(a.chromosome==t.chromosome)&((a.tss-t.tss).abs()<=5000)]
 for neighbor,b in near.iterrows():
  if neighbor==target:continue
  divergent=t.strand!=b.strand and ((t.strand=='+' and t.tss>b.tss) or (t.strand=='-' and t.tss<b.tss));support=guides[(guides.target_gene_name==target)&(guides.nearest_within2kb_nontarget_gene_name==neighbor)];dist=pd.to_numeric(support.nearest_within2kb_nontarget_gene_dist,errors='coerce').abs().dropna();n=support.sgRNA.nunique() if 'sgRNA' in support else len(support);records.append(dict(target=target,neighbor=neighbor,distance=abs(int(t.tss)-int(b.tss)),divergent=bool(divergent),guide_supported=bool(n>=2 and len(dist)>=2),guide_distance=float(dist.median()) if len(dist) else float('nan'),guide_count=int(n),guide_provenance='public CD4 source library; challenge guides unavailable'))
new=pd.concat([old,pd.DataFrame(records)],ignore_index=True);temp=root/'promoter_pairs.extended.csv';new.to_csv(temp,index=False);roundtrip=pd.read_csv(temp);pd.testing.assert_frame_equal(roundtrip.iloc[:len(old)].reset_index(drop=True),old);temp.replace(root/'promoter_pairs.csv');rt.save(root/'promoter_extension_certificate.json',dict(status='PASS unchanged original annotation rows; only new target genes appended',old_rows=len(old),new_rows=len(records),original_sha256=rt.sha(original),extended_sha256=rt.sha(root/'promoter_pairs.csv'),new_targets_disjoint_from_H1=True,guide_download_bytes=downloaded,qualification='Same released lookup rule; static source-guide annotations are not proof of hidden-guide cis effects'));print('Promoter extension',len(records),'rows; original rows exact',len(old),flush=True)
