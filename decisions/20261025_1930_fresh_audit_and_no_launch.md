# Decision — fresh competitor audit and new-benchmark vote

**Date:** 2026-10-25 19:30Z
**Deciders:** Astra (`openai/gpt-6-astra`), Fable (`anthropic/claude-fable-5.1`), Kimi
(`moonshotai/kimi-k3.1-pro`). The agent did not vote.
**Question:** with the final entry frozen, public validation apparently closed, and a new Oct-25
leave-one-cell-line-out benchmark available, should the team spend the remaining 11 days building
and testing on it?

Raw answers:
- `decisions/consult_1025_fresh/astra.md`
- `decisions/consult_1025_fresh/fable.md`
- `decisions/consult_1025_fresh/kimi.md`

## Votes

| advisor | verdict | first experiment |
|---|---|---|
| Astra | Do not commit 11 days; permit a bounded pilot contingent on unanimous approval | one paired frozen-recipe vs source-reliability-weighting run on a preselected K562 holdout, all five anchors, leakage-audited |
| Fable | **No** — not worth 11 days; close out, stand down | write a one-page close-out; no benchmark download or training |
| Kimi | **No** — no actionable official upside; archive and write postmortem | close the record; revisit only if scoring reopens or a separate research objective is approved |

**Result: 2–1 against launching.** Because Astra conditioned her permission on unanimous approval,
there is **no authorized benchmark run**. Nothing was downloaded, trained, repackaged, or submitted.

## Points of agreement

1. Do not replace, resubmit, or alter the frozen entry.
2. The new benchmark is not D/E/F data and cannot produce an official score.
3. A benchmark gain must never be described as an official improvement.
- Any future work requires either a changed validation/selection landscape or an explicit new objective.

## Why the earlier questions had different answers

- **Trash recovery:** the deleted `.vcc` files are not recoverable through git (never committed) and
  were not in a backup. The vault holds a byte-identical copy of the safety package. The superseded
  packages are reconstructible from source and receipts, not byte-identical.
- **Git:** a dedicated local repo now exists at `vcc2026/.git` with two commits; packages, datasets,
  caches, and nested repositories are excluded. No remote was added and nothing was pushed.
- **Audit recency:** `competitors/audit_1025_live/REPORT.md` records the live Oct-25 checks. No
  post-Oct-22 validation entries were found. Remaining gaps are listed in that report.
- **Model progress:** the recorded validation sequence improved from v4a 0.1438 to billmA 0.16438,
  about +0.0206, with no official feedback available since Oct 21. The leader is 0.47203.

## Follow-through

- Verified coverage is 298/300, with only `ANG` and `DPY30` on generic fallback. (The independent receipt is not included in this source export.)
- Frozen package SHA is unchanged:
  `ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`.
- No compute job launched. The guard and leaderboard monitor remain active and read-only.
