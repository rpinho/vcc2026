"""Restore byte-exact frozen metadata; no network or bulk count acquisition."""
from pathlib import Path
import hashlib,json,os,shutil,tarfile
root=Path(__file__).resolve().parent
index=json.loads((root/'artifacts/frozen_metadata_index.json').read_text())
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
archive=root/index['archive'];assert sha(archive)==index['archive_sha256']
assert shutil.disk_usage(root).free-sum(x['bytes'] for x in index['entries'].values())>40*1024**3
with tarfile.open(archive,'r:gz') as a:
 for member in a:
  assert member.isfile() and member.name in index['entries'];path=(root/member.name).resolve();assert path.is_relative_to(root)
  expected=index['entries'][member.name]
  if path.exists():assert sha(path)==expected['sha256'];continue
  path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_name(path.name+'.'+str(os.getpid())+'.tmp')
  with a.extractfile(member) as source,temp.open('wb') as dest:shutil.copyfileobj(source,dest,1048576)
  assert temp.stat().st_size==expected['bytes'] and sha(temp)==expected['sha256'];temp.replace(path)
print('PASS frozen metadata byte identities restored; bulk inputs remain absent')
