# Prior art and publication scope — refreshed 2026-10-23

The 2026-10-23 refresh began by reading the original prior-art review, official
score record and the later public-code census before exporting additional code.
Our best recorded validation result is now **billmA 0.16438**. This is our
reconstruction of a publicly described response MLP, not the original author's
implementation or a reproduction of that author's score.

Fresh public reads include GitHub repository searches sorted by update time for
"virtual cell challenge", "vcc2026", "cell-eval", "perturbation prediction 2026"
and "State model arc"; the official challenge website; Arc's State commit history;
and the cell-eval2 changelog. Sidechain's public README now reports SER-29 at
0.2226, with broader source coverage and context-dependent bulk strength. Its
exact scored settings and assets have not been reproduced here.

- [Sidechain source and reported results](https://github.com/saberhq/sidechain)
- [Official final-phase dates and hidden-score policy](https://virtualcellchallenge.org/)
- [cell-eval2 changelog](https://github.com/ArcInstitute/cell-eval2/blob/main/CHANGELOG.md):
  0.20.0 adds evaluation settings/results replay. The historical experiments in
  this release used 0.18.0; no silent scorer upgrade or new score claim is made.
- [Arc's PIE release](https://arcinstitute.org/news/pie) predates this refresh;
  source/weights remain excluded under their own terms.

The dated 2026-10-18 observations below are retained as history. This refresh
publishes existing experiments, including failures; it makes no new performance
or record claim. No complete hidden official score was reproduced locally.

## Original publication review — 2026-10-18

This export packages existing work; it introduces no new model or performance claim.
The original project prior-art review (2026-10-18) was read before assembling this directory.

- Official leaderboard: https://virtualcellchallenge.org/api/leaderboard?get_final=false . Earlier project audit reported Illumina iCell 0.4708712865; this is a dated observation, not a current record claim.
- Our official baseline: v5b 0.1528, from the project's OFFICIAL_SCORES.md. Local tests do not establish official performance.
- Public baseline: Sidechain SER19 (https://github.com/saberhq/sidechain), previously observed 0.1684131816; AtlasShift (https://github.com/kaipengm2/Virtual-Cell-Challenge-2026), published baseline 0.1546. Full identical reproduction of either official score is not claimed. AtlasShift numerical kernels were reproduced in the original work.
- Scoring: https://github.com/ArcInstitute/cell-eval2 . Use the official package; original work used 0.18.0.
- Stage1-lab ideas: https://github.com/danielchen26/vcc2026-stage1-lab . No license found in the inspected source; no source code is included here.
- Existing review covered GitHub sorted by update time, Arc issues/PRs, the broad AI-lab/results repository sweep, recent primary papers, public community posts and publicly released methods. That research is inherited, not represented as a fresh exhaustive search for this export.
- Representative primary papers: https://arxiv.org/abs/2511.16954 and https://arxiv.org/abs/2604.13986 . No new scientific or leaderboard record is claimed.

Publication checks and current rules review are recorded in PUBLISH_CHECK.md.

Current publication-time aggregate recheck: {"checked_utc": "2026-10-18T14:48:00.188987+00:00", "entries": 1575, "leader_team": "Illumina AI", "score": 0.4708712865331793, "source": "https://virtualcellchallenge.org/api/leaderboard?get_final=false"}. Raw response was processed only in memory; no member information exported.

Public-code heads rechecked read-only at 2026-10-18T14:55:22.994999+00:00:

- https://github.com/kaipengm2/Virtual-Cell-Challenge-2026: d24ce4fdae0cd7cba3ba29546cd8737094eae98a
- https://github.com/saberhq/sidechain: 0a09ee1d29dcf861384c5689d34866a85d786503
- https://github.com/ArcInstitute/cell-eval2: a2630dcbda95fc0b819bf841e56d9561d9c1e1da
- https://github.com/danielchen26/vcc2026-stage1-lab: 3a5ea5e9f640a6dd652516fd1492711618829870

Released AtlasShift numerical fixture was reproduced successfully in a synthetic disposable copy during publication preparation. No new model or record claim is made. The exhaustive research sweep was not repeated because this task exports existing code rather than searching for a new method.
