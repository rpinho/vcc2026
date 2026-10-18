# Publication check — 2026-10-18

## Verdict

**CONDITIONAL GO for a public release of our original code and methods, with the
included MIT dependencies and notices. No blanket during-competition publication
embargo was found in the current rules, FAQ or website terms.** This is an
interpretation of the published text, not organizer approval of a particular
collaboration. Preserve team eligibility, do not coordinate across teams to
undermine the competition, and distribute no restricted data. The export remains
local, with no commits, remote, remote repository creation or push.

## Official evidence: code and method publication

The [FAQ, Overview](https://virtualcellchallenge.org/faq#overview), question about
researchers/students publishing results, states:

> We encourage publishing and sharing research with the scientific community, if that is allowed by your employer or academic institution.

This supports open scientific publication. There is no stated requirement here to
wait until the competition ends. Existing personal-capacity work is the scope of
this export; publication does not override obligations to an institution.

The [Rules, Responsibility to Others](https://virtualcellchallenge.org/rules#responsibility-to-others)
prohibit:

> collude with other teams to share models, results or strategies in a manner that undermines the integrity of the Challenge

Public availability is not a blanket waiver of this integrity provision. Publishing
for everyone and inviting general improvements does not establish permission for
private reciprocal arrangements or evasion of team/submission limits.

The [Rules, Finalists, Winners and Prizes](https://virtualcellchallenge.org/rules#finalists-winners-prizes)
require finalists seeking prizes to provide a public high-level method description;
code and weights disclosure is optional. This is a finalist obligation, not a
publication deadline for other entrants. The rules take precedence over conflicting
guidance. [Official rules](https://virtualcellchallenge.org/rules).

## Team boundaries

The [Rules, Teams](https://virtualcellchallenge.org/rules#teams) specify teams of:

> one to eight

Each individual may participate on only one team, and all members must be eligible.
Employment/university and personal capacities cannot be mixed in one team; employee
teams have an organization restriction with a university exception. Institutional
representation requires authorization. The FAQ directs team-merger requests to
organizers; a public fork does not itself add someone to our registered team.
[Team FAQ](https://virtualcellchallenge.org/faq#teams),
[participation capacity](https://virtualcellchallenge.org/rules#participating-as-employee).

## Data redistribution and rights

The [Terms of Use](https://virtualcellchallenge.org/terms), effective August 4, 2026,
section 1, Using our Content, states:

> The Content is protected by copyright and other intellectual property laws

The terms define Content to include data and explain that downloading it does not
remove those protections. Section 2 restricts reproduction/distribution and similar
uses of the Services:

> except as expressly permitted by these Terms

The terms also protect other users' privacy and require rights in material posted
through the Services. These provisions are not a general redistribution license.
The rules separately require rights to use third-party code/data.
[Website terms](https://virtualcellchallenge.org/terms),
[rights requirement](https://virtualcellchallenge.org/rules#responsibility-to-others).

No separate, explicit permission to republish organizer datasets, hidden responses,
or derived per-cell data was established by this review. Nor was a standalone
universal data-redistribution ban located in those pages. Therefore this release
redistributes **no data at all**: readers obtain authorized challenge inputs from
the organizers and public-source datasets under each provider's terms. Our MIT
license applies to software; it does not relicense any organizer or external data.

## Retrieval and evidence scope

Checked the live public home, /rules, /faq and /terms on 2026-10-18. The page HTML
initially displays a loading placeholder; the actual text was read from its public
Next.js content chunks, without login, account access or browser automation:

| Page | Served content chunk | SHA-256 of fetched JavaScript |
|---|---|---|
| Rules | [2622-4f80a691fca6a841.js](https://virtualcellchallenge.org/_next/static/chunks/2622-4f80a691fca6a841.js) | 4a4ed038b116b95e2a94411a2883f5ce65617e683b50a00bfd22a0d1da7ae5af |
| FAQ | [334-291236033da9ac0d.js](https://virtualcellchallenge.org/_next/static/chunks/334-291236033da9ac0d.js) | 8ab6ab697c9efc18fa4ba8981a0b1912ee7683823c4932aa26f9e100f697deab |
| Terms | [5584-2d1dd1b066f1342c.js](https://virtualcellchallenge.org/_next/static/chunks/5584-2d1dd1b066f1342c.js) | 910fec2f64395b6fcbe8ef6ccc1168f03a5c7171ea8fa49cffb8a375df465fac |

Quoted excerpts above are brief; read the linked full provisions in context.
No private organizer correspondence or account-specific download agreement was
inspected. Recheck rules and any download-specific terms before an eventual release.
This review concerns publication; it is not certification of any model's eligibility
or score. The current rules also distinguish model training from directly inserting
experimental results into a submitted prediction.

## Export scope and license audit

Included: model_v4, model_v5 and model_v6 Python pipeline/preparation/verification
code; proxy scoring/calibration code; model_de experiments and independent count
checks; competitors/monitor.py; RSI acquisition, baseline, experiment-launch and
verification source; selected MIT vendor Python modules; README with our official
score table, requirements, MIT LICENSE and original third-party notices. Only our
aggregate official score table is copied; no raw leaderboard payload is included.

Explicitly excluded:

- Credential stores, tokens, keyring contents, environment files and SSH keys.
- All challenge/public datasets, hidden/restricted data, raw or derived per-cell
  matrices, selected rows/metadata, source banks, embeddings and pretrained weights.
- `.vcc` packages, artifacts, caches, logs, ledger, journals and historical results
  bundles. RSI `record_*.py` journal-writing scripts are omitted.
- `_codex/`, agent instructions, private agent/launchd supervisor and submission helper.
- Raw leaderboard responses, other teams' member lists and upstream site dumps.
- All stage1-lab source; it has no license in the inspected checkout. Ideas only
  are credited. Arc PIE noncommercial/share-alike source and model assets are omitted.
- Two numba-mwu-derived kernels with unspecified upstream license:
  proxy/presorted_mwu.py and model_de/mwu_kernel.py. The export retains the standard
  installed scorer; the optional model_de fast path explicitly refuses activation.

[THIRD_PARTY.md](THIRD_PARTY.md) records MIT notices and reviewed revisions. AtlasShift,
Sidechain and zero-shot-perturbation notices are preserved verbatim. Arc cell-eval2
is an external installed dependency with its MIT notice included. Dependencies
retain their own terms. The original working model/RSI sources were not altered.

Publication-only edits normalize local/example paths, discover the participant CLI
on PATH, remove the v4 implicit cache-deletion hook, disable machine-specific
workload hooks, and use the standard scorer instead of the unlicensed optimization.
The public full-adapter fixture now creates its scratch directory on first use.

## Validation and scans

Final measured results and scan output follow below. Tests run in disposable copies
outside the export, with synthetic data only. Their scratch matrices, receipts,
logs and paths are not included in this release. Full challenge-score reproduction
was not attempted; original source/coverage choices and historical caches are absent.

## Local repository and eventual publication

Directory: `../vcc2026-public/` relative to the original research directory.
It has its own empty Git repository on branch main, no commit and no remote.
An initial commit and publication remain actions for the owner. Nothing was staged,
committed, uploaded or submitted by this preparation.

The following owner-run commands would create and publish the repository; they
have **not** been executed:

```sh
cd ../vcc2026-public
git add .
git commit -m "Publish VCC 2026 pipeline and evaluation harness"
gh repo create rpinho/vcc2026 --public --source=. --remote=origin --push
```

Before running them, inspect `git status --short` and rerun the scans below.

## Final verification results

| Check | Result |
|---|---|
| All exported Python sources parse | PASS (360 files) |
| Exported shell scripts parse with their interpreter | PASS |
| v6 pipeline, bank builder and prep CLI startup | PASS |
| Leaderboard monitor and RSI launcher CLI startup | PASS |
| Released AtlasShift numerical fixture | PASS; exact integer totals and independent rational bounds |
| Bounded DE versus official full scorer/anchors | PASS bitwise equality |
| Bounded scorer versus official dispatcher | PASS bitwise equality |
| DE strategy/count-margin/metric fixtures | PASS; official 1/2-thread and dispatcher parity |
| Proxy full adapter with standard MWU and exact caches | PASS bitwise official full scorer, 400-cell baseline and two-split anchors |
| Fresh local H1-bank → v5b build (synthetic) | PASS; 3 targets × 400 cells × 80 genes; integer counts/structure independently checked |
| Original vendor MIT notices | PASS; all 8 AtlasShift/Sidechain notices byte-identical |
| Nontrivial exact Python-file matches against stage1-lab | 0 (supplemental provenance check, not proof against every possible derivative) |

One clean-checkout fixture initially lacked scratch/output directories. The export
now creates both and its full parity rerun passed. No private precomputed fixture
was copied into the release to make these checks pass.

## Secret and privacy scan result

Gitleaks **8.30.1**, built-in default rules, redacted output, no per-file allow
comments honored. Both installed scanner tools were detected; Gitleaks was used.
The scan is local and performs no credential verification against online services.

Run from the export root:

```sh
env -u GITLEAKS_CONFIG -u GITLEAKS_CONFIG_TOML gitleaks dir . \
  --no-banner --no-color --redact --ignore-gitleaks-allow
python tools/check_export.py
```

Observed Gitleaks result on the completed export (rechecked after embedding this result):

```text
Gitleaks 8.30.1
scanned ~3023415 bytes (3.02 MB) in 5.43s
no leaks found
exit code: 0
findings: []
```

```text
Files scanned: 389; Python sources parsed: 360
Privacy/structure findings: 0
PASS
```

The additional offline check scans every export file outside Git metadata for
emails, personal absolute paths, credential literals and private keys; it checks
Python syntax, unexpected/binary file types, symlinks and runtime/data directories.
The Sidechain Python package named `data` is a source-code directory and is allowed;
no dataset files are allowed. Initial broad patterns mistook that package and an
upstream `/home/` URL fragment for data/a local path; the check was narrowed to
recognize those cases. No credential or personal-data finding was suppressed.

Final result: **PASS — zero privacy/structure findings, zero secret findings**.
No real token store, keychain, environment file or hidden challenge dataset was
opened or copied. Regex/Gitleaks checks are heuristic, not a guarantee that every
possible secret representation can be detected.

SOURCE_MANIFEST.sha256 fingerprints the export other than this audit and the
manifest itself, so the source snapshot can be independently checked. Tests ran
outside this tree; the public export contains no generated test data.
