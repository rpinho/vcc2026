**Verdict:** Run H4; quarantine H5 until its norm’s provenance is clear. I would not approve a final-panel build from a small local win, and I would amend the scaled-MSE gate **now, before seeing results**.

### 1. H4: coherent, but not implied by v9
Independent selection of call support does not imply independence of DE log2FC, signs, confidence ordering, or the amplitude/budget solution.

**Cheapest falsification:** on a few targets, freeze v9 and amplitude, then diff the selected 885 calls and the scorer’s DE inputs across incumbent/`eta0`/`dense`. If they differ, compute only the DE components. A DE advantage for gating/deflation falsifies H4; identical DE inputs mean DE scores cannot move in that frozen comparison.

`eta0` is the right intermediate **if deflation is the main suspect**: incumbent→`eta0` estimates deflation’s effect with the gate on, while `eta0`→`dense` estimates the gate’s effect without deflation. It is an ordered ablation, not a full separation; deflation-on/gate-off is still needed to assess interaction.

### 2. H5: only conditionally well-posed
Training in log-bulk space is insufficient to establish commensurable norms. Both norms must describe **raw response vectors**, before per-target normalization, using compatible genes, centering, and aggregation—and must not include the control baseline merely because the MLP predicts absolute log-bulk. A norm computed after normalization contains no learned magnitude.

**Cleanest alternative:** blend raw responses in native log-bulk space, apply one out-of-fold scalar shrinkage, and let the blended vector determine its own norm.

### 3. Local-build bar and MSE
Treat **0.074 as a pragmatic decision-error scale, not a confidence interval**. My minimum screen would be:

- overall gain **>0.074** over the pre-specified same-run comparator;
- paired target-family stability, with an interval excluding zero;
- no gate failure or pool1-like fidelity/Jaccard collapse;
- subsequently, one pre-specified official-validation transfer check.

Even that would not certify transfer. **pool1’s +0.0017 is a clear no-build.**

`scaled MSE > 0` is not a sound general necessity gate: a capped/floored component can remain zero despite genuine raw-MSE improvement. Record raw/uncapped MSE and per-target changes. Positive scaled MSE is useful evidence of escaping the scorer floor, especially for these expression-motivated arms, but its absence is not proof of no improvement.

### 4. If both fail: ranked remaining experiments
1. **Official-validation transfer test of the exact incumbent versus v5b control**, adding pool1 only if cheap. Their local overall gap is about 0.057, with substantial fidelity/Jaccard disagreement; this is the result most likely to change the package decision and also tests whether the harness has any decision value. If no transfer check is possible, stop surrogate optimization and retain the incumbent.
2. **Fixed-bulk, fixed-amplitude emitter ablation:** v9’s generic 885-call layout versus v5b calls, with at most one target-family-cross-validated call-budget alternative. This attacks the largest observed weakness rather than the saturated MSE floor.
3. **A fully gate-legal boundary amplitude check**, nominally requested ratio 1.05 versus 1.00. It follows the only visible dose-response, but the expected gain is small and remains subject to the bar above.

### 5. Sign-off
I would refuse promotion based on best-of-several-arms selection and “overall strictly up” alone, or any claim that a pre-registered local pass implies blind transfer. Same-run controls remove run drift, not harness misranking. Freeze comparators and non-ablated settings, report every component with paired uncertainty, and document the MSE amendment before unblinding. An overall gain accompanied by DE degradation may still be interesting, but it does **not** support H4’s mechanism.
