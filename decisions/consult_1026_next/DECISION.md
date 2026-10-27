# DECISION — 2026-10-26T01:30Z

Advice on record: `kimi.md`, `fable.md`, `astra.md`. Measurement overrides advice where they conflict;
each override is stated.

## The night's real result (produced by Fable's objection)

Fable noticed PDS 0.706 (official) and cosine 0.0045 (local) cannot both describe one model. I verified
both and found the resolution: **CD4 is a cell line we hold out entirely** (`public_data/` has only
HepG2, Jurkat, H1, RPE1, K562, K562_GWPS), and the shipped package predicts contexts D/E/F whose
targets overlap our local 79 by only 4.

| | local CD4 fold | official panel D/E/F |
|---|---|---|
| PDS raw | 0.497 (none) | **0.836** (strong) |
| MSE scaled | worse than control | **0.061** (real skill) |

Both true: the model generalises to an unseen line poorly, but performs on the scored panel. **Every
"dead" verdict from conclusions 2–5 is therefore a verdict about cross-line generalisation, not about
the competition.** Full write-up: `diagnostics/THE_REAL_DIAGNOSIS_2026-10-26.md`.

## Advisor verdicts — what was accepted, what was overridden

### Kimi — REJECTED on measurement
Proposed adaptive budget `B_i = c·Σ|β|`. Tested her exact rule: `Σ|β|` gives R² 0.006, Spearman 0.054;
calibrated to match the median, its per-target predictions correlate with true `n_conf` at **0.10**.
She assumed the predictor works; it does not. Rejected.

### Astra — three corrections, all ACCEPTED
1. **My `budget_predictability.py` CV was wrong** (`groups` built then unused; it partitioned
   consecutive rows). Re-ran with **leave-one-family-out**: R² = **−0.108**, worse than a constant
   baseline (−0.022). `n_conf` is definitively unpredictable. This *strengthens* the negative result.
2. **The "MSE 0.2325 ⇒ cosine 0.875" inference is invalid.** `R = 1 + a² − 2aρ = (a−ρ)² + 1 − ρ²`, so
   `√(1−R)` is only a lower bound at optimal amplitude, not an identified cosine. **Retracted.**
3. **The "leader has replicates" claim is not evidence.** Correct — it was an inference, now withdrawn.
Her fidelity→Jaccard→NMAE→reach priority ranking is adopted, and her pass thresholds are adopted
(fidelity +0.05 raw, composite +0.02, CI excluding zero, ≥5 of 8 families).

### Fable — the objection that mattered
Her contradiction catch produced the night's key finding. Her recommendation (six-component
reproduction harness) is subsumed: we do not need to reproduce the hidden panel — we established the
regime difference directly by showing CD4 is an unseen line.

## DECISION — what runs, and the stop rule

**Running now:** H1 Arm A, an **in-distribution** fold (H1 is trained on; 400
cells/target, reliability 0.528 vs CD4's 0.317). This is the first evaluation tonight on a fold that
matches the scored panel.

**Pass (Astra's thresholds, preregistered):** mean raw fidelity +0.05, composite +0.02 on validated
anchors, paired CI excluding zero, improvement in ≥5 of 8 families.

**Stop rule:** if H1 does not pass, stop the signals-and-renderer programme and write up. Astra's and
Kimi's independent conclusion is the same: without expression-direction signal on the scored
distribution, DE-set engineering moves us ~20 places, not to the board.

**Not doing:** building an 08:00 submission package (submission is frozen and unavailable — Astra
agrees), amplitude/shape/transfer variants (all tested on the wrong fold), global LFC inflation
(measured null).

Frozen entry `5NeeoWn96lfm23hgExcm` untouched. Nothing submitted, nothing pushed.
