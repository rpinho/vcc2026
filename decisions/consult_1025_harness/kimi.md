## (A) VERDICT: HYBRID

**STOP AND SHIP for the November 05 deadline**, while permitting **one strictly time-boxed, learning-only diagnostic of Fact 2**. Keep `5NeeoWn96lfm23hgExcm` as the last FINAL. Do **not** launch an 11-day rho-model sprint or replace the incumbent with an unvalidated candidate.

## (B) WHAT TO DO, ranked

1. **Freeze the submission now — ≤1 hour.**
   Verify the entry ID, SHA `ce78c81c...9dc4`, FINAL status, and artifact integrity. Confirm no automated process can publish or replace it. The agent submits nothing. “Ship” now means **leave this entry untouched through the deadline**.

2. **Run the Fact-2 diagnostic — hard cap 4 hours, attended, read-only.**
   Use retained predictions and identical inputs to:
   - compare packaged-emitter versus local-replay outputs and baselines;
   - recompute raw \(u\), MSE, PDS, and the anchor conversion from raw arrays;
   - check target sets, preprocessing, scaling, and fold stratification;
   - determine whether an emitter/baseline mismatch or arithmetic error can be excluded.

   Treat officially implied rho only as an aggregate sanity check, not per-target truth. If the necessary retained artifacts are unavailable, stop rather than recreate the pipeline. No new scorer, folds, anchor fitting, training, or overnight jobs. Even finding a bug would **not** authorize a replacement FINAL.

3. **Close the record — ≤1 hour afterward.**
   Mark `learnednorm` **void**, not negative; do not rerun it. Archive the diagnostic and defer ideas #6/#7/#11/#12 until a trustworthy validation setting exists.

## (C) DID THE NEW FACTS CHANGE ANYTHING?

**They change the research conclusions, but not the shipping decision.**

- **Fact 1:** The earlier “harness failure” was too uniformly stated. CD4 and pooled-macro results show the harness is not globally useless. They still do not establish official-panel transfer: H1_rich has only 10.1% target overlap, the pooled test still fails on v4a/v4b, and there are zero exact local/official recipe pairs from which to learn a transfer rule.
- **Fact 2:** Yes, this is worth a **bounded diagnostic**, because it offers concrete implementation hypotheses that existing artifacts may falsify. That is different from rebuilding the harness or attempting to infer official truth from anchors. It cannot validate a new candidate.
- **Fact 3:** The rho-raising search space was **not exhausted**; any earlier implication otherwise should be withdrawn. But “unimplemented and mechanistically relevant” is not evidence of improvement. Eleven CPU days without a yardstick would produce an unselectable artifact. Idea #11 additionally creates imputed predictions for uncovered targets, not independent evidence.
- **Fact 4:** This removes a purported result; it supplies no positive evidence for changing the final.

Thus: **(i) yes, diagnose briefly for learning; (ii) no, do not spend the remaining deadline window pursuing those ideas; (iii) STOP AND SHIP still governs the submission.**

## (D) RISK

The dominant risk is **displacing a live FINAL whose blind score is unknown with a candidate that cannot be officially scored or credibly selected**. A local gain, even after a successful diagnostic, would not establish an official gain. The secondary risk is false confidence from treating fold-dependent anchors or implied rho as transferable ground truth. Unbounded jobs also recreate the disk/CPU failure mode just observed. The four-hour diagnostic cap limits the cost of learning without reopening the submission decision.
