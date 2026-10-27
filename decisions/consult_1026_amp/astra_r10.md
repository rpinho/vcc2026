**Yes—give it a 48-hour feasibility gate, not all nine days.** R9 closes the tested call-side routes. v8a→v8b provides an official amplitude signal: PDS+MSE contributed approximately **+0.018 overall**, before the NMAE loss.

The main risk is **transferability**: increasing B does not guarantee better MSE. Spread must remain truth-aligned, the claim changes with the emitted layout, and the official slope comes from another baseline. The four-cell kernel’s tolerance failure makes representability the first hurdle.

**First experiment:** freeze B0’s per-cell CPM vectors and change only their relative library depths. This preserves CPM means and Wilcoxon inputs mathematically while changing depth-weighted pseudobulks. Solve nonnegative, bounded depth weights to approximate amplified B0 pseudobulk responses. Start with baseline and two modest amplitudes; approach B/claim≈0.95 only if feasible. With fixed claim, moving 0.50→0.95 suggests √(0.95/0.50)≈1.38×, but recompute after emission.

**Gate: DM48 — emitted-channel invariance and budget recovery.**

- Verify the serialized, reloaded counts: NMAE unchanged within a predeclared numerical tolerance; DE calls, signs and ordering unchanged. Mathematical invariance before rounding is insufficient.
- Require increased PDS and report per-context B, predicted jackknife claim, depth extremes and pseudobulk reconstruction error.
- Treat control-resampled **B/claim as a label-free proxy**, not the hidden true ρ.
- Require a central projected gain above +0.01 using the official amplitude slope; report half-credit and zero-credit MSE scenarios separately. Local clamped MSE cannot validate that projection.

Stop at 48 hours if feasibility or invariance fails. Passing earns further validation, **not promotion**: the official projection remains a hypothesis. B0 stays frozen at **0.16438**; any submission remains Ricardo’s decision.