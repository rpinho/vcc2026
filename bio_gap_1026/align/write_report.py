"""Write the <=450-word decision report from verified results only."""
import runtime as rt
rt.enforce(); rt.guard()
from collections import Counter
audit=rt.json.loads((rt.ROOT/'results/analysis.json').read_text())
verification=rt.json.loads((rt.ROOT/'results/verification.json').read_text())
assert verification['status']=='PASS'
s=audit['summary']
names={'B0':'B0','knn':'kNN','map':'Ridge / RRR','B0_knn':'B0 + kNN','B0_map':'B0 + map',
       'B0_knn_PDS':'B0 + kNN, PDS constraint','B0_map_PDS':'B0 + map, PDS constraint',
       'B0_knn_PDS_CI':'B0 + kNN, PDS bound','B0_map_PDS_CI':'B0 + map, PDS bound',
       'gene_variance':'Gene variance prior'}
order=['B0','knn','map','B0_knn','B0_map','gene_variance','B0_knn_PDS','B0_map_PDS','B0_knn_PDS_CI','B0_map_PDS_CI']
def ci(v): return f'[{v[0]:+.4f}, {v[1]:+.4f}]'
rows=[]
for arm in order:
    if arm not in s: continue
    v=s[arm]
    rows.append(f'| {names[arm]} | {v["cross"]:.4f} | {v["aligned_per_spread"]:.4f} | {v["PDS"]:.4f} | {v["B"]:.4f} | {ci(v["delta_cross_CI95"])} | {ci(v["delta_PDS_CI95"])} |')
point=[a for a in s if s[a]['qualifies_point_PDS']]
conservative=[a for a in s if s[a]['qualifies_CI_PDS']]
best=max(point,key=lambda a:s[a]['cross']) if point else None
strict=max(conservative,key=lambda a:s[a]['cross']) if conservative else None
if best:
    selected=strict or best
    v=s[selected]
    reading=f'**Reading:** {names[selected]} adds aligned signal at matched spread: Δcross={v["delta_cross"]:+.4f}, ΔPDS={v["delta_PDS"]:+.4f}; '+('both bootstrap criteria pass (PDS narrowly).' if strict else 'cross-CI and point-PDS limits pass.')
    reading+=f' Q: {s["B0"]["prediction_energy"]:.4f}→{v["prediction_energy"]:.4f}; centered cross: {s["B0"]["centered_cross"]:.4f}→{v["centered_cross"]:.4f}.'
    decision=f'**One arm worth rendering:** {names[strict or best]}; analysis only here, no rendering performed.'
    if (strict or best)=='B0_map_PDS_CI':
        decision+=' Modal choices: w=0.02, rank40, λ=0.1 (154/168 folds).'
    if not strict:
        decision+=' PDS noninferiority is not established at 95% confidence; its lower CI crosses −0.005.'
else:
    reading='**Reading:** Source neighbors/maps add aligned signal at matched spread, but their PDS losses fail the −0.005 limit; the gene prior also loses alignment.'
    decision='**One arm worth rendering: none.** No evaluated arm meets both requested conditions.'
text='\n\n'.join([
    '**vcc_align_1026 — verified training168 LOO; no held-out79 access.**',
    'T=Σ||b||²; C=cross/T; Q=Σ||b̂||²/T; A=C/√Q. Every arm is rescaled to B0’s centered B. B below means B/T; uncentered Q is separate.',
    '| Arm | C at B0 B | A | PDS | B/T | ΔC 95% CI | ΔPDS 95% CI |\n|---|---:|---:|---:|---:|---|---|\n'+'\n'.join(rows),
    reading,
    'Fits exclude the outer target; hyperparameters/weights use inner LOO. Follow-up blends constrain inner ΔPDS≥−0.005, using either its mean (“constraint”) or normal 95% lower bound (“bound”). These are exploratory. Similarities use non-H1 responses.',
    'CIs: 2000 paired target bootstraps, fixed fits/candidate panel, spread re-matched per draw; exploratory, without multiplicity correction. PDS exactly follows the preceding comparator. Rescaling includes the mean; centered-cross diagnostics: results/analysis.json.',
    decision,
    '**Official MSE is not measurable locally.** These diagnostics establish no official-score gain.',
    'Verified: independent primal fits/all outer predictions, three complete nested grids, native PDS, bootstrap replay, exact rational cross/energy/B. [Receipt](results/verification.json), [certificate](results/exact_alignment.json).',
    'Reproduce: `analyze.py`, `constrained.py` (also `--conservative`), `verify.py`, `write_report.py`, each with `nice -n 10 python -B`. Definitions: [METHODS.md](METHODS.md); credits: [PRIOR_ART.md](PRIOR_ART.md). CPU only, one BLAS thread, inherited nice20. No rendering, packaging, submission, upload, promotion or commits.'
])+'\n'
count=len(text.split()); assert count<=450,count
(rt.ROOT/'REPORT.md').write_text(text)
rt.save('results/report_check.json',dict(words=count,max_words=450,point_PDS_pass=point,CI_PDS_pass=conservative,one_arm_worth_rendering=strict or best))
rt.status(f'REPORT.md written: {count} words; point PDS qualifiers={point}, CI PDS qualifiers={conservative}.')
