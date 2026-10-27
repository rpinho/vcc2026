# Decisions index

The main calls we made between 10-21 and 10-26, and what the evidence was. The AI
advisors were Fable and Opus (Anthropic), Astra (OpenAI) and Kimi (Moonshot). An
advisor vote is not a measurement; several rulings below were later corrected by data.

| Date (UTC) | Question | Decision | Votes |
|---|---|---|---|
| 10-21 | billmA scores 0.16438, our best. What next? | Scorer-aware MSE work before the validation phase closes. | Fable, Astra |
| 10-22 | Validation scoring closes; final phase opens on new contexts D/E/F. | Submit billmA as the final entry. A layout variant was also submitted, then replaced on 10-23 by billmA after it lost 0.008 locally. | — |
| 10-25 | Can any local candidate justify replacing billmA? | No. Freeze the entry; the local harness cannot rank candidates reliably. | Astra, Kimi, Fable 3–0 |
| 10-25 | Research continues (owner's call), without touching the entry. | Test amplitude calibration and a low-rank donor→target map conditioned on recipient controls. | Astra, Fable, Kimi |
| 10-26 | Head-to-head on 79 held-out H1 targets (fit on 168). | STRING/manifold and donor maps improved raw MSE but lost the composite. Closed. | — (measured) |
| 10-26 | Magnitude (P1a) and direction (N1) calibration. | Both lost rendered PDS and the composite. Closed. | — (measured) |
| 10-26 | Renderer: seeds, uniform confidence, 2× calls (D1770). | Seeds move the composite ≤0.0008. Uniform confidence lost. D1770 gained +0.034 locally but calibrates to about seed noise; not promoted. | Fable, Opus, Astra 3–0 |
| 10-26 | Is the local harness trustworthy? | It inflates fidelity 9–12× and Jaccard 7–9×, and gets the composite sign right on 2 of 4 checkable contrasts ([calib3](calib3_1026/REPORT.md)). | — (measured) |
| 10-26 | Call-side screen: 250-call budget, reach head, detectability selection. | All three neutral or worse ([R9](h2h2_1026/results/r9_gate.md)). Closed. | — (measured) |
| 10-26 | Decoupled amplitude (raise pooled spread, keep per-cell values). | Whole-cell reweighting can't do it ([dm_1026](dm_1026/REPORT.md)); depth-tilt not built. | Fable, Opus stop; Astra build. 2–1 |
| 10-26 | Better biology? | Large gains need measured knockdowns in the same cell type, which D/E/F don't have; what transfers is tiny ([bio_gap_1026](bio_gap_1026/REPORT.md)). Closed for this entry. | Fable, Opus, Astra 3–0 |

The final entry is billmA. Final-phase scores are hidden until after 11-05.
