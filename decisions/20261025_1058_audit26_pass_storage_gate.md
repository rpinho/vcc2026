# 2026-10-25 10:58Z — AUDIT26 PASS; bridge authorized, production decode held by storage gate

## Audit result

AUDIT26 is **PASS** on revision 24. It reviewed exactly the 16 bridge
files and one rebuild support file named by the manifest, and all 17 live and
17 frozen copies matched before and after testing. The report's 17-row hash
table was rechecked against the live bytes; every row matches.

- Manifest: `package_bridge_1025/audit26_scratch/revision24_manifest.json`
- SHA-256: `3d4c833dbf6b1487d4094f8b1f7abd437ed85e3517b54af10a8a6b488d68b4c0`
- Frozen suites: **252/252** cases, including repeated receipt and timestamp runs
- Independent/supplemental probes: **115/115**
- Calendar counterexamples: all three return 15 with no decoder entry through
  both direct orchestration and the isolated launcher tail; valid controls
  accept

AUDIT24 and AUDIT25 remain FAIL and are not used as authorization. The bridge
code is now audited, but that is not a production decode result.

## Storage constraint

The decode requires a snapshot and temporary decoded output. **Decision: defer until sustained storage headroom and bounded competing writes can support the full run.** Recheck the 17-file identity and safety-package hash before launch. Transient machine-specific capacity readings are omitted.

## Cleanup

Four superseded packages are eligible for deletion. Each was rehashed and
matches its recorded receipt:

- `model_v7/v7a.vcc` — failed gate G1
- `model_v7/v7b.vcc` — failed gate G1
- `model_billm/billmA.vcc` — superseded by `final_billmA`
- `model_billmB/out/final_billmA_repro.vcc` — byte-identical rebuild proof

Their sizes, hashes, gate failures, `final_receipt.json`, 900/900 comparison,
independent monitor and audit evidence remain. **The safety package
`model_billmB/out/final_billmA.vcc` is not touched.** Deletion is limited to
these four files and is recorded in `package_bridge_1025/CLEANUP_1025.md`.

## Consultation

**Kimi (`orq best`, `consult_1055/kimi.md`):** delete all four; preserve hashes,
sizes, receipts, G1 failures, rebuild evidence and AUDIT26. Defer until the reserve, snapshot, peak decode, competing growth and allowance can all be supported throughout the run.

**Astra (`codex exec -m gpt-6-astra`, `consult_1055/astra.md`):** same deletion
eligibility and same launch conclusion. It adds that a prior successful capacity observation is not a proven peak, and that the competing writer must be
identified and bounded before launch. It also flags that the bridge and
independent decoder use different GB/GiB floor checks.

## Safety

Nothing submitted, uploaded or replaced. Safety entry
`5NeeoWn96lfm23hgExcm`; package 2,595,584,000 bytes with SHA-256
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`.
No production decode has run.
