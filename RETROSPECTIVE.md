# What we tried in the Virtual Cell Challenge 2026

Our best validation score was 0.164. When validation closed on October 22, that put us around 463rd of about
1,580. The top ten were around 0.32 and the leader was at 0.47, so we were nowhere near the top. That entry,
billmA, is our final submission. Final scores come out after November 5.

I ran most of this with AI agents and spent a lot of tokens on it. Here's what happened.

## Submissions

| Entry | Score | What changed |
|---|---:|---|
| v1 | 0.0877 | Public signatures for the same target, added to the control mean |
| v3b | -0.2482 | The same average shift on every target, plus shrunk signatures. Bad idea. |
| v4a | 0.1438 | Started from AtlasShift, a public method |
| v5b | 0.1528 | Matched the spread to how the scorer credits MSE |
| v9a | 0.1613 | DE calls ordered by confidence |
| billmA | 0.1644 | Our version of a public response-MLP idea |

Most of the gain came from v1 to v4a, when I stopped building on my own first model and started from a public
one. I should have done that on day one. After that it was small steps. The full list is in
[OFFICIAL_SCORES.md](OFFICIAL_SCORES.md).

## After October 22

Once the server stopped scoring, the only way to compare ideas was to run the official scorer locally on
training targets we held out. We assumed that would be close enough. It wasn't. The local scorer said the v9
layout change was a clear win (+0.066). The official board had already scored the same change as a small loss
(-0.009). We had four changes scored both ways, and the local scorer got the direction right on two of them.
That's a coin flip. ([calibration](calib3_1026/REPORT.md))

We kept trying things anyway: gene networks, magnitude and direction calibration, more DE calls, fewer DE
calls, renderer changes, amplitude changes. Some looked good locally. None held up once we corrected for how
far off the local scorer was. So billmA stayed.

## The biology

I think this is the interesting part. On 168 training targets, our predicted response barely pointed the right
way. The mean cosine with the true response was 0.084. Two halves of the same real data agree at about 0.65, so
this isn't noise. We were at roughly a tenth of what's achievable ([noise ceiling](bio_gap_1026/ceiling/REPORT.md)). It wasn't missing data, because the targets we had
source data for carried 99.5% of the error. And most fixes that lowered the error did it by shrinking the
predictions toward zero, which the MSE part of the score doesn't reward.

The one thing that helped a lot was using other knockdowns measured in the same cell type. On the scorer's
alignment measure it went from 0.115 to 0.42. The problem is that the final test is three new cell lines with
only control cells, so that data doesn't exist there. When we held out whole cell types to mimic that, it went
from 0.115 to 0.118. ([biology diagnostics](bio_gap_1026/REPORT.md))

My guess is that some of the gap to the leaders on the validation board comes from same-cell-type data that
nobody will have in the final phase. I could be wrong. We'll see when the final scores come out.

## The AI advisors

Decisions went to a panel of models (Fable, Opus, Astra, Kimi). They caught real mistakes, some of them mine.
They also agreed with each other on things that turned out to be wrong, more than once. After a while I stopped
treating a unanimous vote as evidence and only trusted the test results. The main decisions are in
[DECISIONS_INDEX.md](DECISIONS_INDEX.md).

## What I'd do differently

Start from the best public baseline. Before trusting a local scorer, check that it ranks your past submissions
the way the official board did. Ours wouldn't have passed. And keep one cell type untouched until the end.

The code and reports are in this repo. There's no data or model weights.
