## (A) VERDICT: what should run next

**Run one bounded, read-only research job on #7, preceded by a provenance check. Then stop unless its artifacts justify another experiment.** Do not schedule eleven days of work merely because the machine is available.

The submission decision remains unchanged: no replacement, no new submission candidate. These experiments can produce future research evidence—not evidence that anything would improve the frozen entry.

The representativeness question remains **unmeasured**. Missing response truth is not repaired by complete feature coverage.

## (B) WHAT TO DO, ranked, first action + hours

These are **elapsed-time caps**, not runtime predictions.

1. **Establish a trustworthy execution record — first action, cap 2 hours.**
   Inspect the instrument report against its referenced files. Create a read-only input manifest, an executable run specification, and a separate research output directory. Record unavailable dependencies as unavailable.

   Resolve the inconsistent entry identifiers in this request using the existing authoritative receipt or platform record. Do not silently select one. If unresolved, record the discrepancy and leave submission-related operations disabled.

   Deliverable: a human-inspectable run bundle with exact inputs, hashes, commands, environment, planned outputs, and stopping rules.

2. **#7: response-reliability source weighting — cap 6 additional hours.**
   First establish whether the existing summaries support a leakage-safe comparison against an identifiable incumbent baseline. Specify the weighting rule, comparison, and falsification criterion **before** examining results.

   If evaluable, run the smallest paired baseline-versus-weighting test. Report target-level results, coverage, uncertainty, and sensitivity to a few influential sources.

   If only descriptive analysis is possible, label it descriptive and stop. A favorable result on the known problematic harness is only a local signal—not validation of ranking reliability or submission benefit.

   Deliverable: a reproducible baseline/delta table, underlying observations, and explicit applicability limits. A documented failure of feasibility is an acceptable result.

3. **#6: only the untested STRING-space/ESM2/PIE-gene-text combination — optional, cap 12 additional hours.**
   Start only after a human reviews the preceding bundle and confirms an evaluable question, frozen comparison protocol, and usable baseline. Prefer a small ablation over architecture search. Do not repeat the rejected STRING-only bilinear core.

   Feature alignment alone does not satisfy this gate. If evaluation is not defensible, skip the job.

4. **#11: residual-value feasibility only — optional, cap 2 hours.**
   Inventory whether representative response truth exists for a meaningful residual-value test. Do not substitute feature availability or performance on source-covered targets for evidence about the 161 no-source targets.

**After these caps, default to idle.** Further work requires a new, inspectable evidence-based decision—not an orchestrator summary.

## (C) WHAT SHOULD *NOT* RUN

- Submission rebuilding, replacement, or mutation.
- The infeasible representativeness replay, unless its missing artifacts actually become available.
- Repeated optimization against the unreliable ordering yardstick.
- Broad hyperparameter sweeps, speculative long CPU jobs, or #12 without a defined falsifiable test.
- Synthetic “truth” presented as representative response evidence.
- Repeated ceremonial audits or invented monitoring infrastructure.

Idle time is preferable to producing misleading confidence.

## (D) HOW SHOULD PROGRESS BE VERIFIED, given the failure mode below

Each job must produce a **self-contained evidence bundle**:

- Immutable input manifest and hashes.
- Exact executable specification, code version, environment, and seeds.
- Captured stdout/stderr, timestamps, exit status, and timeout record.
- Raw observations, computed tables, and explicit missing-data/failure records.
- An independent checker that recomputes reported comparisons from raw outputs and verifies input identities.

Run jobs without submission credentials or write access to frozen-package files. Use bounded automatic termination; do not assume overnight supervision.

Human checkpoints should inspect artifacts directly and rerun at least one calculation. Status should be generated from execution records; every numerical claim must cite its source file.

**Exit 0 proves completion, not scientific validity. Hashes prove identity, not truth.** Neither an orchestrator narrative nor its own `verification.json` is sufficient evidence by itself.