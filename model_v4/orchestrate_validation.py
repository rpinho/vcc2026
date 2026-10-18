"""Wait for bounded read-only acquisitions, then execute frozen attribution arms."""
import runtime,time,subprocess,os
rt=runtime.ROOT;start=time.time();last_checkpoint=0
required=['HCT116_transfer.json','HEK293T_transfer.json','CD4_transfer.json','K562_raw_transfer.json','HEK293T_control_correction.json']
while not all((rt/'artifacts'/p).exists() for p in required):
 runtime.guard()
 if time.time()-last_checkpoint>900:
  done=[p for p in required if (rt/'artifacts'/p).exists()];(rt/'STATUS.md').write_text('# v4 checkpoint\n\n'+time.strftime('%Y-%m-%d %H:%M UTC',time.gmtime())+'\n\nBounded public source acquisition active. Completed: '+str(done)+'. Code, independent exact certificates, and frozen attribution protocol are present. Data-independent AtlasShift fixture and estimator checks have run. Biological quality remains UNVERIFIED. Next: published-configuration anchor, isolated ablations in native-depth and 20k H1 views, then select two distinct predictions. No submission or credentials. Free GiB: '+str(round(runtime.shutil.disk_usage('/').free/1024**3,2))+'\n');last_checkpoint=time.time()
 time.sleep(5)
print('All acquisition receipts present',flush=True)
python=os.path.abspath('../../.venv/bin/python');children=[]
def run(args,log):
 with open(rt/'logs'/log,'w') as f:
  p=subprocess.Popen([python,*args],stdout=f,stderr=subprocess.STDOUT);print('child PID',p.pid,args,flush=True);children.append(p)
  try:code=p.wait()
  finally:
   if p.poll() is None:p.terminate();p.wait()
 print('child result',code,args,flush=True);return code
assert run(['correct_orion_controls.py'],'correct_orion_controls.log')==0
assert run(['source_provenance.py'],'source_provenance.log')==0
assert run(['prepare_public_validation.py'],'prepare_public_validation.log')==0
# Import only after all acquisitions end; no numeric work while waiting.
from evaluate import VARIANTS,eval_folder
for fold,view in [('H1','native'),('HCT116','native'),('HEK293T','native'),('H1','depth20k')]:
 if fold!='H1' and not (rt/f'artifacts/validation_{fold}.h5ad').exists():continue
 if (eval_folder(view,fold)/'reference_invalid.json').exists():
  print('Excluded undefined reference',fold,view,flush=True);continue
 for variant in VARIANTS:
  code=run(['evaluate.py','--fold',fold,'--view',view,'--variants',variant],f'eval_{fold}_{view}_{variant}.log')
  (rt/'STATUS.md').write_text('# v4 checkpoint\n\n'+time.strftime('%Y-%m-%d %H:%M UTC',time.gmtime())+'\n\nAttribution active: '+fold+' / '+view+' / '+variant+' exit '+str(code)+'. Per-arm logs and JSON preserve completed results. Entire destination family excluded; quality UNVERIFIED. Next: inspect failures and attribution, choose best and distinct runner-up, build and verify archives. No submission or credentials.\n')
  if code and (eval_folder(view,fold)/'reference_invalid.json').exists():
   print('Excluded undefined reference',fold,view,flush=True);break
  if code and variant=='anchor' and fold=='H1' and view=='native':raise RuntimeError('Native AtlasShift sanity anchor failed; inspect before further attribution')
print('VALIDATION STAGE COMPLETE',flush=True)
