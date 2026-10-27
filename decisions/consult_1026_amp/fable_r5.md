**Ruling, round 5.**

One fact from the scorer settles Q1. The layout step is constrained to preserve per-gene column sums and per-cell depth (score_run.py line 71-72). So budget and confidence cannot move PDS or raw MSE at all. They touch only reach, fidelity and jaccard, where the paired CI half-width is about ±0.005 and the largest effect ever observed is +0.0008. The +0.054 bonus lives entirely in the emission transform, not the 885 calls.

**1. Experiment.** Run Astra's uniform-confidence ablation, and only as a diagnostic, not a candidate. It is the single-factor test of whether the layout does anything. Fable's 2x budget doubles a knob with no evidence the first dose works. Pre-register the reading: composite within ±0.002 of B0 means the layout is inert and the renderer line closes.

**2. Legitimacy.** Changing budget or confidence is legitimate if the frozen spec lists them as parameters, but legitimacy is not a reason. To replace B0 a renderer change would need: paired composite CI excluding zero on the 79 targets, PDS and raw MSE not worse than B0, the same sign on a second disjoint target split, and a bit-exact re-render matching the submitted hash. Nothing within reach of a layout change can clear the first bar.

**3. Stop.** Yes, stop. Nine levers failed, the remaining lever is provably inert on the two metrics that matter, and the final phase punishes a late swap with hidden scores. Finish the one diagnostic arm, then declare B0 final and spend the week on submission hygiene: hash check of the uploaded file, a clean re-render from the frozen spec, nothing else.
