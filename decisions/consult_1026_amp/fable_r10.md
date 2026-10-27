**1. Yes, conditionally. Nine days is enough for a cheap kill test plus one emitter build, and nothing else has an official datum above +0.01.** Main risk: the two channels are not as separable as the brief assumes. DE reads Wilcoxon ranks, not means. Depth-correlated placement holds the CPM mean but widens within-group CPM dispersion, so borderline calls shift and reach/fid move. The 10-22 kernel failure on the 4-cell pooled template says the usable depth variance is small, so the achievable pseudobulk gain may be well under 1.5x. Second risk: the one datum (v8b) is confounded, and the credit model already mispredicted candA by 4x.

**2. First experiment (half day, label-free): compute B, Σ min(Ĉ,C) and ρ for the emitted B0 package exactly as the scorer does, with C from jackknifed control resamples at n_p.** Kill rule: if ρ ≥ 0.9, the budget does not bind, the lever is dead, B0 final. If ρ ≈ 0.5 as expected, build the decoupled emitter on the v9 layout at pseudobulk gains 1.25x and 1.5x, borrowing counts from uncalled genes in deep cells so group library size holds.

Gate name: **DM-R10**. Pass requires all of:

| readout | threshold |
|---|---|
| DE call set per target | identical to B0 |
| local fid, jac, reach | within seed floor (±0.003) |
| local NMAE | ≥ −0.002 scaled |
| local PDS | ≥ B0 |
| emitted B/claim | ≥ 0.9 |
| MSE projection | +0.079 × (ΔB/claim ÷ 0.5), then halved for the candA miss |

Submit only if projected overall ≥ +0.008 net of the DE/NMAE deltas. That is half of v8b's PDS+MSE and the margin for MSE being locally unreadable. Timeline: ρ readout 10-27, emitter 10-30, gate 11-01, decision by 11-02.

**3. Not applicable, but for the record:** the only other datum above +0.01 is reach via call ordering (v8: +0.085 scaled), already banked in billmA. If DM-R10 fails any row, B0 final.
