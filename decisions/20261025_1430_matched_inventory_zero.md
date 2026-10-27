# 2026-10-25 14:30Z — zero exact-matched local/official records

## What ran

Astra's ranked option 2: a bounded, read-only inventory of retained local/official matched
records (`_codex/vcc_matched_inventory_1025.md`, job `vcc_matched_inventory_1025`). It was
explicitly forbidden from fitting any promotion rule or promoting approximations.

## Result: 0 exact matched contrasts

67 records inventoried (63 distinct local full-six-component emissions from 80 retained
rows/aliases, plus 4 official-only identities):

| bucket | records |
|---|---:|
| EXACT-MATCHED | **0** |
| APPROXIMATE | 0 |
| UNVERIFIED | 23 |
| NO-OFFICIAL | 44 |

Zero in APPROXIMATE is not a promotion of the known surrogates — they stay UNVERIFIED with
`documented_approximation: true`, because their official emitter/selection bindings cannot be
resolved from retained evidence.

## Why each family fails to match

- **v1–v4b:** local recipe/prediction/implementation/reference/selection digests are retained, but
  historical official emitter/source and complete adaptive-selection bindings are unresolved. The
  local `exact_target_matched: true` flag means target/batch matching *inside the local task*, not
  local/official recipe identity.
- **v5a/v5b:** score→spec→audit chain checks, but the frozen specs describe 79-target
  available-input ports, not the hidden 300-target submission, and carry
  `external_official_label: null`. The v5b spec itself says 0.1528 remains unpaired.
- **Recent four arms:** v9a_885 is a generic 885-budget surrogate for v9a; v9a_694 the analogous
  694-budget surrogate for v8a; v8b_like uses 1.5 amplitude at fixed 885. The v5b arm repeats the
  historical H1 emission and adds no observation.
- **Native billm:** model_billm and model_billmB have different H1/CD4 prediction digests; both are
  native-emitter LOCO sanity checks and neither proves equality to the officially scored billmA v9
  emission. The 885-call incumbent is a third adaptation.
- **Family-quality H1:** closure receipt PASS, but the scientific decision is NOT PROMOTED —
  requested-ratio difference 0.051027 exceeds the 0.02 tolerance. The CD4 continuation was never run.
- **Sidechain:** replay equality is numerical repeatability only, `transfer_noise_bound: null`; its
  `audit.json` is a copy of the historical comparison, not a new cohort.
- **Safety package:** blind scores are hidden; the two submission IDs are aliases of one package, so
  they add zero independent observations.

## Consequence

Astra's prior conclusion is now confirmed by direct inventory rather than inference: there is **no**
retained pairing that could serve as an independent transfer observation. No acceptance rule is
proposed, and none can be derived from these records.

Minimum missing observation, stated by the job: one frozen recipe-faithful candidate-minus-incumbent
contrast in the intended modification family, with exact head/source/amplitude/emitter/per-target
call-budget identities and selection history for both arms, linked to independent official component
outcomes on a versioned panel/anchor bundle. Even that would permit a first transfer-error check, not
a reusable acceptance bound — independent prospective errors and a validation-to-final bridge remain
missing. The job names historical billmA vs v9a as a concrete recovery target but explicitly does not
authorize a replay or submission.

## Integrity

- 6,983 input paths hashed before and after; identical path sets, no changed entries.
- All 18 original audit input hashes match current files.
- No model, scorer, or verifier module imported; no fitting, emission, or scoring.
- Safety package 2,595,584,000 bytes
  `ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`, unchanged.
## Decision

Close the transfer-calibration question for this deadline: within the bounded retained inventory there
is no certified exact pairing on which to build a promotion rule, and manufacturing one would require
new official-scored contrasts that this deadline does not accommodate. Retain the safety package. Do
not reopen the score search on the basis of this inventory.

## Correction after independent review

The independent reproducibility check (`closeout_repro_check_1025/REPRO_CHECK.md`) correctly flagged
that "the evidence base does not exist and cannot be manufactured" overstates a bounded inventory.
What the data support is narrower and still decisive for the decision:

- **Established:** 0 of 67 records pass the four-field exact-identity test, and the inventory
  disclaims any full-repository absence proof.
- **Judgement, not measurement:** that no promotion rule should therefore be built. This remains my
  recommendation, but it rests on the certification result plus the deadline, not on the number alone.
- **Loose-reading sensitivity:** under a procedural definition the v5a/v5b LOCO ports are arguably
  recipe-matched. They still fail the purpose, because they are 79-target available-input ports with
  `external_official_label: null` — no official counterpart exists to pair them with.

The decision to retain the safety package is unchanged.

No submission, upload, or package replacement. Safety entry `5NeeoWn96lfm23hgExcm` unchanged.
