# dm_1026 — decoupled-amplitude feasibility (consult_1026_amp R10), 2026-10-26 ~20:15 UTC
B0's cells were re-emitted for 5 H1 targets (counts_sha256 matches render_B0.json on all 5). Label-free throughout.

**Astra's mechanism (per-cell integer depth multipliers k_c; CPM vectors fixed, so DE/NMAE-invariant) is
INFEASIBLE.** Pooled log1p(5e4) delta along B0's own delta, relative to B0 (amp), and orthogonal error (orth),
both in units of |delta|:
- Threshold weights (probe.py; K in {2,3,4,6}, top 15-50% of cells by projection): best amp 1.03 at orth 0.20;
  max amp 1.12 at orth 0.77 (cosine 0.83). 5/5 targets alike.
- NNLS continuous weights aimed at s=1.2 / 1.38 (probe_nnls.py; an upper bound for any reweighting): amp
  1.00-1.01, orth 0.21-0.25 / 0.45-0.49 (ACAT2, GSK3B, SRC).
Reason: the emitted cells are draws around B0's predicted profile. Their spread is sampling noise, not variation
along the delta, so no reweighting of whole cells moves the pooled profile along the delta.

**Remaining mechanism (Opus R10 / Fable R10): depth-tilted placement.** Per gene, move counts toward high-depth
cells with per-cell zero-sum compensation. The pooled sum rises while the per-cell CPM mean is held. Efficiency
per moved count is 1 - d_lo/d_hi. Emitted depths: p10/p25/p50/p75/p90 = 34.9k/44.6k/55.6k/69.3k/85.1k (efficiency
0.36 IQR, 0.59 at 10/90). Arithmetically feasible, but it changes within-group CPM dispersion, so the Wilcoxon
ranks (DE calls, reach, fidelity) may move: the risk both named. Not built.
