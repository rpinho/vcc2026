"""Final read-only process/file checks and completion receipts."""
import runtime,subprocess,os,time
rt=runtime.ROOT;root_pid=38256;me=os.getpid();rows=[]
for line in subprocess.check_output(['ps','-axo','pid,ppid,pcpu,command'],text=True).splitlines()[1:]:
 q=line.strip().split(None,3)
 if len(q)==4:rows.append((int(q[0]),int(q[1]),q[2],q[3]))
def children(parent):
 result={parent}
 while True:
  more=result|{p for p,pp,_,_ in rows if pp in result}
  if more==result:return result
  result=more
own=children(root_pid)-{root_pid}-children(me)
if os.getppid()!=root_pid:own.discard(os.getppid())
infra={p for p,pp,cpu,cmd in rows if cmd.endswith('/bin/codex-code-mode-host')}
workers=[p for p,pp,cpu,cmd in rows if p in own-infra];assert not workers,workers
orphans=[p for p,pp,cpu,cmd in rows if p not in children(me) and p!=root_pid and str(rt)+'/' in cmd and ('/.local/bin/vcc prep ' in cmd or 'spawn_main' in cmd)];assert not orphans,orphans
stamp=time.strftime('%Y-%m-%d %H:%M UTC',time.gmtime());runtime.save(rt/'artifacts/shutdown_verification.json',dict(checked_utc=stamp,command='ps -axo pid,ppid,pcpu,command',own_scientific_workers=workers,own_path_orphan_workers=orphans,excluded='Checker/tool process and application-owned Codex code-mode host',status='VERIFIED no own scientific workers remain'))
files=sorted(str(p.relative_to(rt)) for p in rt.rglob('*.vcc'));assert files==['out/prediction_v4a.vcc','out/prediction_v4b.vcc'];size=sum(p.stat().st_size for p in (rt/'artifacts').rglob('*') if p.is_file());assert size<50_000_000;free=runtime.shutil.disk_usage('/').free;assert free>40*1024**3
(rt/'STATUS.md').write_text('# v4 complete\n\n'+stamp+'\n\nBoth archives passed installed vcc prep (exit 0) and independent full decoded CSR/gene/context/target equality audits. A=anchor_single_mean (best public proxy); B=iv_ash_dual_guide (full-integration ablation). Biological/hidden score UNVERIFIED; published Atlas score not locally reproduced. Full attribution, controls/source coverage, licensing credits and resource/process deviations are in REPORT.md. Only two final .vcc files remain; large intermediates removed, reusable public summaries retained under the authorized shared directory. Retained artifact data '+str(size)+' bytes. Free disk '+str(round(free/1024**3,2))+' GiB. ps check: no own scientific workers remain. No submission or credentials.\n')
with (rt/'REPORT.md').open('a') as f:f.write('\n**Final cleanup and shutdown ('+stamp+').** Only the two requested `.vcc` files remain; regenerable matrices, head caches, scorer caches and job-local dependencies were removed. Retained certificate/data files total '+f'{size:,}'+' bytes (<50 MB), excluding the explicitly retained shared public summaries. Free space at completion was '+f'{free/1024**3:.2f}'+' GiB. The required process-list check found no own scientific workers remaining (`artifacts/shutdown_verification.json`); application tool infrastructure is excluded.\n')
print(stamp,'files',files,'artifact bytes',size,'free GiB',round(free/1024**3,2),'workers',workers,flush=True)
