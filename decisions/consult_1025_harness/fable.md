**(A) VERDICT: Keep the finalized version unchanged. Do bounded diagnostic work only — with a hard rule that no new version is submitted before the deadline, whatever the diagnostics show.**

The decisive fact is not in your list of four: the scoring server is closed. Any new version is a blind swap of a validated 0.164 for an unmeasured number, ranked by a harness that is (i) fold-dependent, (ii) possibly off by 2x on the key quantity, and (iii) recently wrong on every attempt it endorsed. Your base rate for "harness says better → actually better" is currently 0 for n>0. Replacing a known score under those conditions is not a bet, it's a coin flip with your own stake and no payout visible until after you've lost the option to revert. Diagnostics, by contrast, cost nothing to the finalized version and set up the next cycle.

**(B) WHAT TO DO, ranked (≈20 operator-hours total)**

1. **Freeze (1 h).** Verify the published artifact is the latest version, make it read-only, remove submission credentials from the working machine, and write the no-submit rule into the project log. This is the one action that must happen today.

2. **Resolve the 2x discrepancy (5–6 h).** This is the only diagnostic whose outcome changes future decisions. Test the three hypotheses in order of cheapness: (a) recompute the anchor conversion by hand with the anchor values and spreads written out; (b) run the exact published artifact through the harness and confirm the emitter/baseline path is byte-identical to what was scored; (c) compare per-fold difficulty against the 10%-overlap fold. Stop when one hypothesis explains ≥0.1 of the gap. Can run unattended overnight in a sandbox directory.

3. **Fix the harness fold design (3–4 h).** Pool the folds, report rank accuracy with a confidence interval rather than three point estimates, and construct one fold whose target overlap exceeds 10%. Do not pick the flattering 18/21 fold as primary.

4. **Rerun the aborted experiment (2 h setup, overnight compute).** Fix the disk guard, run in the sandbox, harness-only. Void results deserve one rerun; nothing more.

5. **Spec the four untried ideas (2 h writing; optional 6–8 h implementing one).** Write each as a one-page plan with a falsifiable harness test. If time remains after 1–4, implement feature-based regression in harness-only mode. It is not submitted. It exists to be the first thing measured when the server reopens.

Stop the moment the operator budget is spent. Hour 21 is not worth more than hour 20 here.

**(C) Did facts 1–4 change the advice?**

There is no earlier advice from me in this thread, so I'll answer the question as asked: no, they would not have changed it, and they shouldn't change yours. Facts 1, 2, and 4 are all about how well candidates can be *ranked*; they are irrelevant when nothing can be *scored*. Fact 3 is the tempting one — "the right ideas were never tried" — but "never tried" is zero evidence of payoff, and theory saying scaling knobs can't raise the ceiling is not theory saying these four will. The closed server dominates every one of these facts.

**(D) RISK**

- *Primary risk: drift.* Diagnostics produce an encouraging harness number, someone argues "just this once," and the validated score is swapped for an unmeasured one. Mitigated by step 1, done first, done physically.
- *Accidental overwrite* of the finalized artifact by an unattended job. Mitigated by sandbox + read-only + separate credentials.
- *Sunk diagnostics:* the 2x gap may remain unexplained after 6 hours. Acceptable; cap and record.
- *Reputational:* 0.164 vs median 0.315 is poor and stays poor. That is already true; a blind submission cannot fix it and can only make it worse or unknown.