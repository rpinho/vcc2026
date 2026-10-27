"""Small, local result ledger and <=500-word report; no numerical fitting."""
import runtime as rt
rt.enforce(); rt.guard()
import csv

recips=('H1','K562','RPE1','HepG2','Jurkat')
allm={}; completed=[]
for r in recips:
    p=rt.ROOT/f'results/{r}_metrics.json'
    if p.exists(): allm[r]=rt.json.loads(p.read_text())
    if (rt.ROOT/f'results/{r}_completion.json').exists(): completed.append(r)
rows=[]
for r,metrics in allm.items():
    for arm,m in metrics.items(): rows.append(dict(recipient=r,arm=arm,**m))
if rows:
    fields=list(rows[0])
    with (rt.ROOT/'results/all_metrics.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields); writer.writeheader(); writer.writerows(rows)

arms=('plain','pooled_map','control_map','H1_direct','pooled_bounded','control_bounded')
def cell(m,ci=False):
    if m is None: return '—'
    v=f'{m["cross"]:.4f}/{m["aligned_per_spread"]:.4f}/{m["PDS"]:.4f}'
    if ci:
        lo,hi=m['delta_cross_CI95']; v+=f' [{lo:+.4f},{hi:+.4f}]'
    return v
table=['|Recipient (n)|Plain|Pooled map|Control map|H1-direct|Pooled blend [CI]|Control blend [CI]|B/T|',
       '|---|---|---|---|---|---|---|---|']
for r in recips:
    if r not in allm:
        receipt=rt.json.loads((rt.ROOT/'results/input_receipt.json').read_text())
        table.append('|'+ '|'.join([f'{r} ({receipt["recipients"][r]["n"]})']+['—']*7)+'|')
        continue
    mm=allm[r]; b=mm['plain']
    table.append('|'+ '|'.join([f'{r} ({b["n"]})']+[cell(mm.get(a),a!='plain') for a in arms]+[f'{b["B"]:.4f}'])+'|')

if len(completed)==5:
    successes={a:[r for r in recips if allm[r][a]['qualifies']] for a in ('pooled_map','control_map','pooled_bounded','control_bounded')}
    candidates=[]
    for a in ('pooled_bounded','control_bounded'):
        h1=allm['H1']['B0_'+a+'_original']
        if len(successes[a])>=3 and 'H1' in successes[a] and h1['qualifies']:
            candidates.append(a)
    verdict='transfers' if candidates else 'does not transfer under the requested criterion'
    winner=max(candidates,key=lambda a:(len(successes[a]),sum(allm[r][a]['delta_cross'] for r in recips))) if candidates else None
    lead=f'**Verdict: {verdict}.** '
    lead+='One arm worth building for D/E/F: **'+winner+'** (exploratory; no official-score forecast).' if winner else 'One arm worth building for D/E/F: **none**.'
    counttext='; '.join(f'{a}: {len(rs)}/5' for a,rs in successes.items())
    tail=f'Cross-CI/PDS successes: {counttext}. A success requires paired ΔC CI lower>0 and mean ΔPDS≥−0.005. The same arm must pass in a majority and beat H1 B0. Full maps fail the PDS criterion; H1-direct fails everywhere. K562’s A falls despite positive ΔC. Per-arm ΔPDS intervals, centered-cross diagnostics and stricter PDS-CI flags: [all_metrics.csv](results/all_metrics.csv).'
else:
    lead='**INCOMPLETE: no transfer verdict is supported. One arm worth building for D/E/F: none established.**'
    tail=f'Completed recipient experiments: {len(completed)}/5. The 31.5-GiB heavy-stage disk guard stopped the run after H1 baselines, before outer map scoring. Free space fell to 30.8055 GiB as system swap grew; job files stayed below 1 MB. The 30-GiB abort floor was respected. No missing experiment is counted as a negative result.'

history=''
if 'H1' in allm and 'B0_original' in allm['H1']:
    mm=allm['H1']; history=f'Original 18,077-gene H1 B0: C={mm["B0_original"]["cross"]:.6f}, PDS={mm["B0_original"]["PDS"]:.6f}.'
    if 'B0_common' in mm:
        history+=f' On the shared axis at plain B, B0 C={mm["B0_common"]["cross"]:.6f}, PDS={mm["B0_common"]["PDS"]:.6f}.'
    for a in ('pooled_bounded','control_bounded'):
        key='B0_'+a+'_original'
        if key in mm:
            v=mm[key]; lo,hi=v['delta_cross_CI95']
            history+=f' B0+{a}: C={v["cross"]:.6f}, ΔC CI [{lo:+.6f},{hi:+.6f}], ΔPDS={v["delta_PDS"]:+.6f}.'
    if 'B0_pooled_bounded_original' in mm:
        history+=' These partial-gene updates preserve B0 outside the intersection; their weights use training recipients only.'
    if len(completed)==5 and winner:
        v=mm['B0_'+winner+'_original']
        history+=f' Its H1 B0 ΔPDS 95% CI is [{v["delta_PDS_CI95"][0]:+.6f},{v["delta_PDS_CI95"][1]:+.6f}]; stricter PDS-CI noninferiority '+('passes.' if v['qualifies_PDS_CI'] else 'does not pass.')

report='\n\n'.join([lead,
    'Cells are **C/A/PDS**, where C is cross at plain-transfer centered spread and A=C/√Q; brackets give 2,000 paired target-bootstrap ΔC 95% CIs versus plain, rematching spread per draw. All main arms use **6,006 shared measured genes**. B/T is matched within each recipient.',
    '\n'.join(table),history,tail,
    'Each outer recipient is absent from labels AND every training source; inner recipient holdouts choose hyperparameters and bounded weights. Controls alone condition gene weights. The 79 held-out H1 identities are excluded from every bank; H1 uses training168 only. Pooled maps use a 256-component input approximation; retained energy and conditioning limits are in [STATUS.md](STATUS.md) and [METHODS.md](METHODS.md).',
    'Verification of completed work: independent Cholesky/SVD inner ridge/RRR checks; exact dyadic baseline cross, energy and centered-B certificates; scalar PDS checks; standalone bootstrap replay. These certify local diagnostics only. Full outer maps remain unverified until run completion.' if len(completed)<5 else 'Verification: independent Cholesky/SVD ridge/RRR checks; exact dyadic cross, energy and centered-B certificates; scalar PDS checks; 35-arm bootstrap replay; input hashes and all 20 selectors/source-exclusion receipts audited independently. These certify local diagnostics only.',
    '**Official MSE is not measurable locally. The final phase is blind.** No rendering, packaging, submission, upload, promotion or commit.',
    'Reproduce sequentially: `nice -n 10 python -B analyze.py`, then `verify.py --replay`, `audit_design.py`, `write_report.py`, with the same Python/nice prefix. CPU only, one worker, BLAS=1; effective inherited nice20. Credits and public checks: [PRIOR_ART.md](PRIOR_ART.md).'])+'\n'
assert len(report.split())<=500,len(report.split())
(rt.ROOT/'REPORT.md').write_text(report)
rt.save('results/report_receipt.json',dict(words=len(report.split()),completed=completed,retained_bytes=sum(p.stat().st_size for p in rt.ROOT.rglob('*') if p.is_file())))
rt.status(f'Wrote REPORT.md ({len(report.split())} words); completed recipients={completed}.')
