# Independent VCC consult, 2026-10-21

You are one of two independent senior advisers (the other is a different lab's model; you will not see each other's answer).
The team requested a plan to improve its Virtual Cell Challenge 2026 result.

Read first (all under the repository root unless noted): PUSH_1021.md (today's honest audit: what was tried vs not),
OFFICIAL_SCORES.md (every official score with metric components), CONSTRAINTS.md, RULES.md if present, model_v7/NOTES.md,
model_billm/PACK_REPORT.md, the unpublished top-500 competitor study (competitor code hunt),
the unpublished adaptation study, rsi/NEXT_SUBMISSIONS.md, and the official scorer/metric code or docs you can
find in the repo (search for the metric definitions: pds, mse noise credit / claim, nmae, fid, reach, jac, scaling).

Facts: best official 0.16438 (billmA, rank 463 of ~1,580): pds .706, mse .061, nmae .116, fid -.066, reach .186, jac -.016.
Ranks 10-50 median scaled MSE ~.44; top-10 median overall ~.32; leader .47 (proprietary data). Two validation submissions
per day (window resets 00:00 UTC). v7a/v7b (scorer-aware emission, never scored) are being submitted 2026-10-22 00:01 UTC.
A Codex job (model_billmB/) is retraining billmA on many more donor rows. The final phase is blind: new contexts D/E/F,
new 300-target panel, only the last entry counts, due 2026-11-05; final controls appear ~Oct 22. Local proxy scores have
proven unreliable (inverted ordering on 10-21), so the official scorer is the only judge. Compute: one Mac Mini (CPU/MPS,
shared execution) plus GPU rentals up to ~$20 each.

Your task (READ-ONLY on everything except your one output file; do not run long jobs, do not download data, do not submit):
1. Diagnose from first principles why our scaled MSE is ~.06-.095 when ranks 10-50 reach ~.44. Work out exactly how the
   official scaled MSE (with its noise credit) rewards a prediction, and what an emitter must do to score well on it
   without wrecking nmae/fid. Be concrete and quantitative; cite the files/lines you relied on.
2. Say what you think the top teams are doing that we are not, using the competitor descriptions/notes in the reports.
3. Give a day-by-day plan to Nov 5 for the two daily validation slots (each a single-factor test with a stated
   hypothesis and the metric component it should move), plus the final-phase plan. Rank levers by expected official gain
   per day of work, and say plainly which ideas to drop.
4. Point out any errors, wasted effort or wrong assumptions in our existing work.
Write your answer (max ~1,500 words) to the output file named in your prompt and stop.
