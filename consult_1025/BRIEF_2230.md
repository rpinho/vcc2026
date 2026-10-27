# Consultation 2026-10-24 22:30 ET — what to run after the expression-unmask factorial

Read-only; answer <=900 words; do not run jobs, download, or submit.

## Hard measured facts (do not re-derive)
- Official best billmA 0.16438: pds .706, mse .061, nmae .116, fid -.066, reach .186, jac -.016.
  Official overall = unweighted mean of the six scaled components (verified, residual 1.7e-4).
- Final phase is BLIND (D/E/F, 300 targets); only the last entry counts; deadline 2026-11-05.
  Safety entry 5NeeoWn96lfm23hgExcm stays; nothing is submitted tonight.
- Local-to-official transfer for layout: fidelity attenuates ~10x, Jaccard ~6x, reach amplifies ~1.5x.
  Local MSE is layout-blind (ratio_of_sums invariance). Local composite has inverted official ordering twice.
- The local harness anchors (b=0.9587, rep=0.0049) differ from official (b=0.9896, rep=0.0366).

## Closed (measured; do not propose retesting)
Global amplitude/rho sweep, pool_k, global bulk multipliers, B/claim, PC projection alone,
learned-blend coefficient, source breadth, per-target call-budget shrinkage, dense layout
(replicates locally but loses officially), shared/specific deflation at TS=50000.

## Currently running
Expression-leg 2x2 factorial (incumbent / no_mask / no_pc / unmasked) on H1 then CD4:
does removing the DE mask and/or the 50-control-PC deflation restore usable response direction?
Incumbent replay already reproduced bitwise (profile SHA b03e05c3..., max_abs_error 0.0).

## The two leads I want ranked
**Lead A — family-quality source fusion (already measured on H1, unpromoted by a gate misfit).**
In mse_direction_1024 (H1 SEARCH, 79 targets, same run, official scorer), `source_scorer_best`
(= family_quality fusion: sqrt(min(n_cells,1000)) weight per family) gave:
raw PDS .75365 -> .78562, raw MSE .98354 -> .96807, NMAE +.00337 raw, fidelity +.0191 scaled,
reach -0.0095 scaled, jac -0.0022 scaled; anchor-scaled deltas summed to **+0.1002 scaled = +0.0167 local overall**.
It was rejected because the pre-registered gate demanded scaled MSE >= +0.030 (it got +0.0162).
Directional rho at the scorer coordinate rose 0.121141 -> 0.143296 (+0.0222, paired 95% CI +0.0137..+0.0308).
A CD4 check showed -0.0019 rho, but **CD4's source pool contains exactly ONE family** (CD4 Rest/Stim8hr/Stim48hr
are one correlated family), so a *family-level* fusion change is a structural no-op there — the "fails to
replicate" verdict may be an artifact of the recipient choice, not evidence of non-transfer.
Available H1/HCT116 task source pools contain 5-6 families (CD4, HCT116, HEK293T, HepG2, Jurkat, K562,
K562_GWPS, RPE1), so a genuine second recipient exists.

**Lead B — per-gene-expression-decile amplitude at fixed call set.** Mechanism: v8a->v8b global x1.5 bought
+0.079 scaled MSE and +0.028 PDS but lost -0.161 NMAE, net -0.006 overall. NMAE is scored on the reference's
own significant set with a fixed denominator, so the correct per-gene amplitude, not a global one, minimises it.
Untested: an amplitude that is larger where MSE is decided (high-expression/high-|lfc| genes) and unchanged or
smaller on the reference-DE set. Hard part: it must be label-free (no truth), e.g. deciles of control
expression or deciles of predicted |lfc|, with the same B=claim solve per arm.

## Constraints
Local Mac Mini CPU only, one heavy job at a time, ~4 workers max. ~8.5 days to deadline.
Each arm costs ~1 h per recipient for 79 targets (H1) and considerably more for 300-target tasks.

## Questions
1. Rank Leads A and B by expected official gain per hour, and say which single one to run next once the
   current job finishes. If neither, say what beats them with the specific mechanism and metric component.
2. For the winner, state the exact frozen arm set, the honest local gate that would not repeat the MSE-only
   gate misfit, and which second recipient is legitimate given the family-structure argument above.
3. Name the specific trap in each lead where our own history says a local measurement will mislead.
