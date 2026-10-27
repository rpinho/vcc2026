**Delete eligibility: yes. Defer the decode until sustained storage headroom and controlled competing writes are established.**

1. **Those four files are eligible for deletion**, provided no active job still needs them. Preserve:
   - Each file’s exact path, size, SHA-256, and receipt establishing failed G1, supersession, or byte identity.
   - The rebuild’s `final_receipt.json`, 900/900 comparison and independent verification evidence, attempt records, logs, and retained source inputs needed to reproduce it.
   - AUDIT26, revision 24’s manifest, all 17 frozen files, and supporting test/probe results.
   - **`model_billmB/out/final_billmA.vcc`**, its identity receipts, and entry `5NeeoWn96lfm23hgExcm`.

   The duplicate rebuild package can go; its proof must remain. Its deletion means future verification of those particular bytes relies on the retained identical package and historical evidence.

2. **Defer until storage capacity is dependable.** AUDIT26 establishes audited behavior, not production capacity or a completed decode. A prior successful capacity observation is not a proven peak requirement; decoded size and competing write rate are not supplied. Actual reclaimed space must be measured.

3. **Minimum launch gate:**
   - Live chain still matches the audited 17-file manifest; required receipts and inputs remain valid.
   - Safety package matches the pinned size and hash.
   - Sustained capacity on the actual output filesystem, with competing writes paused or bounded sufficiently for the run. Use consistent units: the bridge and independent decoder implement different GB/GiB floor checks.
   - Run through the audited launcher; accept only a current-run bridge PASS qualified by a fresh independent decode.

**Next concrete action:** identify the competing writer and establish an exclusive or bounded storage window. Then perform separately authorized cleanup and measure headroom before launch. If that window cannot be secured tonight, defer.

No files changed.