import runtime,sys,os
os.environ['RAYON_NUM_THREADS']='1';sys.path.insert(0,'vendor/deps')
import lance
from public_store import PublicStore
p=PublicStore();e=lance.dataset('s3://public/data/HCT116/expression.lance',version=682,storage_options=p.options())
f=e.get_fragments();print('fragments',len(f),'rows',e.count_rows(),flush=True)
rows=[]
for x in f:rows.append(dict(id=x.fragment_id,physical_rows=x.metadata.physical_rows,path=x.metadata.files[0].path))
runtime.save('artifacts/HCT116_fragments.json',rows)
print(rows[:3]);print('tail',rows[-3:]);print('near',[r for r in rows if 368<=r['id']<=386]);p.close()
