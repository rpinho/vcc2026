import runtime as rt
rt.enforce();rt.guard(True)
import os,subprocess,shlex
provenance=rt.json.loads((rt.ROOT/'results/provenance.json').read_text())
for path,want in provenance['pins'].items():assert rt.sha(path)==want,('input changed',path)
snapshot=subprocess.run(['ps','-axo','pid,ppid,pcpu,command'],capture_output=True,text=True,check=True).stdout
(rt.ROOT/'completion_ps.txt').write_text(snapshot)
active=[]
scripts={'reconstruct.py','render_bulk.py','analyze.py','inventory.py','cpu_baselines.py','benchmark_cpu.py','verify.py','reproduce.py'}
for line in snapshot.splitlines()[1:]:
    fields=line.split(None,3)
    if len(fields)!=4:continue
    pid=int(fields[0]);command=fields[3].split()
    if pid==os.getpid() or not command or 'python' not in rt.Path(command[0]).name:continue
    if not any(rt.Path(c).name in scripts for c in command):continue
    cwd=subprocess.run(['lsof','-a','-p',str(pid),'-d','cwd','-Fn'],capture_output=True,text=True).stdout
    if '\nn'+str(rt.ROOT)+'\n' in cwd:active.append(dict(pid=pid,command=command))
assert not active,active
report=(rt.ROOT/'REPORT.md').read_text();assert len(report.split())<=600
files={str(p.relative_to(rt.ROOT)):dict(bytes=p.stat().st_size,sha256=rt.sha(p)) for p in sorted(rt.ROOT.rglob('*')) if p.is_file() and p.name not in ['completion_ps.txt','final_audit.log','completion_audit.json']}
retained=sum(v['bytes'] for v in files.values());assert retained<50_000_000,retained
rt.save('results/completion_audit.json',dict(status='PASS',active_owned_workers=active,all_frozen_inference_inputs_rehashed=True,input_pin_count=len(provenance['pins']),retained_bytes=retained,free_gib=rt.guard(),report_words=len(report.split()),files=files,no_submission_package_upload_or_promotion=True,unresolved=['Exact full-training Wilcoxon n_conf not recoverable from local moments.','Essentiality uses screen membership, not H1 fitness.','HipSci/LINCS complete coverage unknown from restricted local metadata.','Official MSE effect awaits instrumentation; no official gain claimed.','Requested nice10 inherits effective nice20; priority increase denied.']))
rt.status('Complete within wall cap: exact integer B0 decomposition on 168; source inventory; fixed LOO fixes; CPU benchmark; independent verification PASS. REPORT.md 599 words. Exact n_conf unavailable; no official MSE claim. Frozen inputs rehashed; no owned workers remain; retained <50 MB; forbidden package/cache untouched.')
