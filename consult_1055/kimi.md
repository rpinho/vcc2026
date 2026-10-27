1. **All four listed superseded packages may be deleted**, including `final_billmA_repro.vcc`: its byte identity is already recorded. The two-package cap is a maximum, not a requirement to retain two. Never delete `model_billmB/out/final_billmA.vcc`.

   Preserve hashes, sizes, original receipts, G1 failures, rebuild sources/configuration, the 900/900 comparison and verification evidence, and `final_receipt.json`. Retain AUDIT26, its manifest, frozen sources and test/probe results. Record the deletions separately; preserve historical receipts unchanged and distinguish historical verification from present file availability.

2. **Defer until sustained storage headroom with competing writes controlled.** If growth is unbounded or the margin cannot be maintained, defer. AUDIT26 clears audited code; it does not establish a production result.

3. **Minimum launch gate:**
   **Available space must cover the reserve, snapshot, peak decode/temp footprint, bounded competing growth and allowance.**

   Maintain the reserve throughout, not merely at launch; a momentary reading is insufficient.

   First preserve evidence and log cleanup; arrange a quiet window with the other job’s owner, without deleting its files. Recheck safety-package identity and all 17 live hashes against the audited manifest. Use the audited launcher with fresh run identity/timestamps, monitor space, and stop safely before breaching the floor. Require the actual production qualification result before claiming PASS; clean only that run’s temporary files.

**Recommendation: cleanup, then a controlled storage window; otherwise defer tonight. No submission.**