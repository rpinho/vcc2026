"""R9 official scoring (AMENDMENT_R9.md). Re-uses score_run.py's confidence_all/run with three textual
substitutions only (confidence cache path; layout budget per arm). Generated source saved for audit."""
import os,sys,ast,hashlib,json
import score_run as sr
j=sr.j
assert j.sha(j.ROOT/'AMENDMENT_R9.md')==(j.ROOT/'results/amendment_r9.sha256').read_text().split()[0]
src=(j.ROOT/'score_run.py').read_text();tree=ast.parse(src)
seg={n.name:ast.get_source_segment(src,n) for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in('confidence_all','run')}
gen='\n'.join([seg['confidence_all'],seg['run']])
subs=[("j.ROOT/'scratch/confidences.npz'","j.ROOT/'scratch/confidences_R9.npz'"),("int(base['renderer']['budget']),own)","BUDGET.get(arm,int(base['renderer']['budget'])),own)"),("budget=885,","budget=BUDGET.get(arm,int(base['renderer']['budget'])),")]
sr.BUDGET={'K250':250}
for old,new in subs:
  assert gen.count(old)>=1,old;gen=gen.replace(old,new)
(j.ROOT/'score_r9_generated.py').write_text(gen+'\n')
exec(compile(gen,'score_r9_generated.py','exec'),sr.__dict__)
if __name__=='__main__':
  assert os.getpriority(os.PRIO_PROCESS,0)>=10
  free=j.shutil.disk_usage(j.ROOT).free;assert free>=31_000_000_000,('Launch requires 31 GB decimal',free)
  assert list((j.ROOT/'artifacts/de_cache').glob('*.parquet')),'Warm DE cache missing'
  ra=json.loads((j.ROOT/'results/r9_arms.json').read_text())
  for a,h in ra['sha256'].items():assert j.sha(j.ROOT/f'results/predictions_{a}.npz')==h,a
  pins=json.loads((j.ROOT/'results/operational_resume.json').read_text())['predictions'];assert j.sha(j.ROOT/'results/predictions_B0.npz')==pins['B0']
  j.save(f'results/authorized_start_R9_{os.getpid()}.json',dict(utc=j.utc(),pid=os.getpid(),free_bytes=free,amendment_sha256=j.sha(j.ROOT/'AMENDMENT_R9.md')))
  finish=sr.memory_policy.install()
  try:sr.run(['RH20','K250','DET885'])
  finally:finish()
