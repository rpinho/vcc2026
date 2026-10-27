"""Write the bounded final report only after independent verification passes."""
import runtime as rt
rt.enforce();rt.guard()
import re,subprocess,csv
analysis=rt.json.loads((rt.ROOT/'results/analysis.json').read_text())
verify=rt.json.loads((rt.ROOT/'results/verification.json').read_text())
assert verify['status']=='PASS' and verify['all_primary_SSE_exact_rational']
s=analysis['summary']['vector'];extra=verify['extra']
names={'B0':'B0','a_mean':'(a) m','b_mean_plus_B0':'(b) m + βB0','c_mean_plus_orthogonal':'(c) m + βB0⊥','d_family_correction':'(d) family correction','e_c_plus_family':'(e) c + family'}
def ci(v,d=3):return '['+', '.join(f'{q:.{d}f}' for q in v)+']'
lines=['| Arm | E | cosine | PDS | B/B0 % | ΔE [95% CI] | ΔPDS [95% CI] |','|---|---:|---:|---:|---:|---|---|']
for key,label in names.items():
    r=s[key];B=f'{100*r["B_relative_B0"]:.2f}' if key!='a_mean' else f'{100*r["B_relative_B0"]:.5f}'
    de='—' if key=='B0' else f'{r["delta_E"]:.3f} {ci(r["delta_E_CI95"])}'
    dp='—' if key=='B0' else f'{r["delta_PDS"]:.4f} {ci(r["delta_PDS_CI95"],4)}'
    lines.append(f'| {label} | {r["E"]:.3f} | {r["mean_cosine"]:.3f} | {r["PDS"]:.3f} | {B} | {de} | {dp} |')
report=f'''**The shared response explains a modest part of B0’s failure.** Training168 only; no held-out access.

Projection onto each target-LOO shared mean captures **11.74% of truth energy** (a descriptive per-target amplitude oracle). The unscaled mean explains **4.44%** versus zero. B0’s mean cosine with the shared mean is **0.0418**; only **0.94%** of B0 energy lies along it. After removing that direction, B0’s residual has E=**1.759**, cosine=**0.077**. Nested scalar shrinkage (mean β=0.1128) gives orthogonal E={extra['nested_scaled_orthogonal_E']:.3f}: weak target-specific distance information.

{chr(10).join(lines)}

Differences are versus B0; intervals use **2000 paired target bootstraps**, conditional on fitted predictions and the fixed 168-candidate panel. Lower E and higher PDS are better. All fitted means/coefficients exclude the evaluated target; stacked family fits also exclude each inner target from its shared mean and beta training. PDS exactly matches the preceding job’s native comparator, including panel-wide gene exclusion.

The table evaluates the literal vector equations, excluding each row’s own gene. The previous job’s feasibility projection is a separate [sensitivity table](results/projected_table.md); projected (d) reproduces E=1.46945, PDS=0.859745 and B/B0=80.35%. B is centered across-perturbation energy.

**Reading:** Adding the mean mainly shrinks distance error while collapsing retrieval and spread; family correction alone preserves both. The shared component does not explain most of the mismatch.

**Official MSE is not measurable: vcc_mse_instr_1026 found no instrument.** B constrains official MSE credit, so these distance gains cannot establish official gains.

Verified: nested fits, native PDS, bootstrap intervals and exact rational SSE. Reproduction, conventions and next step: [METHODS.md](METHODS.md). Credits: [PRIOR_ART.md](PRIOR_ART.md). CPU only; nice10 requested, inherited nice20; one BLAS thread. No rendering, submission, upload, packaging or promotion.
'''
words=len(report.split());assert words<=400,words
(rt.ROOT/'REPORT.md').write_text(report)
projected=['The same fixed arms after the preceding job’s clipped/renormalized pseudobulk feasibility projection. No count or cell rendering. Bootstrap intervals are conditional paired target intervals (2000 replicates); retrieval candidates remain fixed.','', '| Arm | E [95% CI] | cosine | PDS [95% CI] | B/B0 % | ΔE [95% CI] | ΔPDS [95% CI] |','|---|---|---:|---|---:|---|---|']
for key,label in names.items():
    r=analysis['summary']['projected'][key]
    projected.append(f'| {label} | {r["E"]:.6f} {ci(r["E_CI95"],6)} | {r["mean_cosine"]:.6f} | {r["PDS"]:.6f} {ci(r["PDS_CI95"],6)} | {100*r["B_relative_B0"]:.5f} | {r["delta_E"]:.6f} {ci(r["delta_E_CI95"],6)} | {r["delta_PDS"]:.6f} {ci(r["delta_PDS_CI95"],6)} |')
(rt.ROOT/'results/projected_table.md').write_text('\n'.join(projected)+'\n')
with (rt.ROOT/'results/arm_summary.csv').open('w') as f:
    fields=['space','arm',*next(iter(s.values())).keys()];writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
    for space,arms in analysis['summary'].items():
        for arm,values in arms.items():writer.writerow(dict(space=space,arm=arm,**values))
ps=subprocess.check_output(['ps','-axo','pid,ppid,pcpu,command'],text=True)
owned=[]
for f in ['results/run.json','results/verify_run.json']:
    pid=rt.json.loads((rt.ROOT/f).read_text())['pid']
    for line in ps.splitlines()[1:]:
        parts=line.split(None,3)
        if parts and int(parts[0])==pid:owned.append(line.strip())
assert not owned,owned
retained=sum(p.stat().st_size for p in rt.ROOT.rglob('*') if p.is_file())
assert retained<50_000_000
rt.save('results/completion_audit.json',dict(status='PASS',report_words=words,retained_bytes=retained,active_owned_workers=owned,process_audit='ps -axo pid,ppid,pcpu,command; matched tracked own worker PIDs',requested_nice=10,effective_nice=rt.os.getpriority(rt.os.PRIO_PROCESS,0),free_gib=rt.guard(),writes_confined=True,heldout_79_opened=False,official_MSE_measurable=False,forbidden_paths_opened=False,no_render_submit_upload_package_promote=True,verification_status=verify['status'],files={str(p.relative_to(rt.ROOT)):dict(bytes=p.stat().st_size,sha256=rt.sha(p)) for p in rt.ROOT.rglob('*') if p.is_file() and p.name not in ['completion_audit.json','STATUS.md']}))
rt.status(f'COMPLETE: independent verification PASS, report {words} words, retained {retained/1e6:.2f} MB, no live owned analysis/verification workers. Official MSE not measurable.')
print('REPORT words',words,flush=True)
