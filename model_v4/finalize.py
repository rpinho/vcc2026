"""Sequential local build, installed prep, exact decode verification. No submissions."""
import runtime,os,subprocess,time
rt=runtime.ROOT;py=os.path.abspath('../../.venv/bin/python')
# Explicit known acquisition/validation processes; do not inspect or stop another job.
for pid in []:
 while True:
  try:os.kill(pid,0)
  except ProcessLookupError:break
  runtime.guard();time.sleep(5)
def run(args,log):
 (rt/'STATUS.md').write_text('# v4 checkpoint\n\n'+time.strftime('%Y-%m-%d %H:%M UTC',time.gmtime())+'\n\nFinalization active: '+str(args)+'. Native public biological quality UNVERIFIED. No submission or credentials.\n')
 with (rt/'logs'/log).open('w') as f:
  p=subprocess.Popen([py,*args],stdout=f,stderr=subprocess.STDOUT);print('worker',p.pid,args,flush=True)
  try:code=p.wait()
  finally:
   if p.poll() is None:p.terminate();p.wait()
 assert code==0,(args,code)
run(['select_models.py'],'selection.log');run(['verify_selection.py'],'selection_verify.log')
s=runtime.json.loads((rt/'artifacts/selection.json').read_text())
for arm in ['a','b']:
 if (rt/f'artifacts/prediction_v4{arm}_package.json').exists():continue
 if not ((rt/f'out/prediction_v4{arm}.h5ad').exists() and (rt/f'artifacts/prediction_v4{arm}_structure.json').exists()):
  run(['build.py',arm,'--variant',s[arm]['variant']],f'build_{arm}.log')
 run(['prep.py',arm],f'prep_{arm}.log')
 run(['verify_package.py',arm],f'package_verify_{arm}.log')
run(['report_results.py'],'report.log')
print('FINAL ARCHIVES AND REPORT VERIFIED',flush=True)
