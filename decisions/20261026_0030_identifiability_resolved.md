# 2026-10-25 22:30Z — Blocker resolved unanimously: recipient-CONTROL-conditioned map (option B), with A nested

## The blocker (found by `dirmap_1025`, not guessed by me)
Fable's original spec said to fit the map onto "the target line's own top perturbation-response PCs". But the
held-out family is excluded from all treated responses, so its response basis cannot be estimated — eligible
own-line treated rows = 0 in all 8 folds. The job correctly refused to silently impute it and stopped at
n=0/8 scored folds (`dirmap_1025/REPORT.md`). That is the right behaviour and I kept it.

This invalidates the *response-basis* ingredient, not the learned-directional-map hypothesis.

## Votes (`decisions/consult_1025_identifiability/`)
- **Astra** (`astra.md`): **B — recipient-control-state-conditioned low-rank map.** Rank-4 bilinear map on
  recipient control-PC coordinates x pooled response-PC coordinates. Learn response PCs and map parameters
  from the seven training families only; fit a **common control-PC coordinate system on training controls**
  and project the held-out recipient's legal controls into it (avoids treating independently rotated PCs as
  comparable). Keep incumbent amplitude/pooling unchanged. Kill if any held-out treated data enter fitting,
  or if mean composite fails to improve across the two incumbent-comparable folds, either fold regresses,
  mean direction_fidelity does not improve, or mean Jaccard declines.
- **Kimi** (`kimi.md`): **B.** Control-state PCs from unperturbed wells; bilinear W on training-family
  treated responses; score on CD4/H1 only. Kill if direction_fidelity gain <5% absolute.
- **Fable** (`fable.md`): **B, with A nested inside as the mandatory ablation**, plus a **B-shuffled** arm
  (control state permuted across families) to separate real recipient-state signal from added capacity.
  Three arms, identical code/seeds. Kill B if the B-real minus A-shared direction-cosine CI includes zero, or
  if B-shuffled reaches within 50% of B-real's gain. Kill the whole map if neither beats the incumbent
  composite on **both** CD4 and H1.

**Unanimous: option B.** Fable's nesting is a strict improvement — it converts "is recipient state worth it?"
into a paired comparison and controls for capacity, at near-zero extra cost. I adopt it.

## Decision
Relaunch as **`dirmap2_1025`** with three arms and no ambiguity left to interpretation:

1. **A-shared**: one map + one response-PC basis from the seven training families; conditioning term zeroed.
2. **B-real**: same map, plus a rank-r bilinear term driven by the held-out family's own **unperturbed
   control** PCs, projected into a training-family control-PC coordinate system.
3. **B-shuffled**: B-real with control states permuted across families (capacity control).

Frozen rules: incumbent amplitude and pooling untouched (the amplitude scalar stays dead); held-out family
contributes **controls only**; rank/ridge chosen by inner leave-one-family-out on training families;
all 8 folds for map-level diagnostics (per-gene cosine, Jaccard proxy); official composite only on CD4 and H1
(the only folds with eligible frozen checkpoints); paired bootstrap over targets; preregister before scoring.

## Honest expectation, stated up front
Only 2 folds can bear the official-composite gate, so even a pass is weak evidence (n=2). If B does not beat
A on the map-level cosine with a CI excluding zero, the whole learned-map family is likely dead on this
panel, and the next move is representation-level (sequence/STRING features) rather than more linear maps.

## Constraints
Frozen entry `5NeeoWn96lfm23hgExcm`, `preserve_1025/` untouched. No submission, no incumbent retrain.
Write only in `dirmap2_1025/`. CPU only, `nice -n 10`. Copy (never launch) other jobs' scripts.
