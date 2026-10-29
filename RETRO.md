# Retro

The [Virtual Cell Challenge 2026](https://virtualcellchallenge.org/) asks one thing. Arc Institute knocks genes
down with CRISPRi and measures every gene in single cells. Can you predict what a few hundred knockdowns do in a
cell type you have not seen? The final test is three new cell lines. You only get control cells. Scores come out
after November 5.

## What we did

We took measured responses from public screens and moved them onto the new cell type's control cells. My own first
model scored 0.088. Starting from a public method,
[AtlasShift](https://github.com/kaipengm2/Virtual-Cell-Challenge-2026), got 0.144. Then three small steps: spread
calibration, gene calls gated on agreement between screens, and a small response network. That got us to 0.164,
about 463rd of 1,580 on the validation board. The top ten were near 0.32. The leader, Illumina, was at 0.47.

## What broke

We could tell knockdowns apart. We could not say what they do.

On 168 training targets our alignment with the real response was 0.084 mean cosine. Two halves of the same real
cells agree at 0.65. So we were at about a tenth of what the data allows. On those same targets our discrimination
score was 0.86. A prediction can be good enough to pick the right knockdown out of a lineup and still barely point
the right way.

It was not missing rows. Targets with source data carried 99.5% of the error.

What helps is in-context data. Other knockdowns measured in the same cell type lifted alignment from 0.115 to 0.42.
The final phase withholds exactly that. We held out whole cell types to mimic it, and 0.115 became 0.118.

After validation closed we were blind. Our local copy of the official scorer got the direction of a change right on 2 of 4 comparisons. It inflated one component 9 to 12 times. We tried size and direction
calibration, more calls and fewer calls, and a small test of network features. None of it held up.

## What I would do differently

Build the local scorer first, and calibrate it against the real board while the board still answers.

Spend the time on expression error. That was our biggest gap, about 0.06 against roughly 0.44 for entries ranked 10
to 50. We did work on it, mostly calibrating size and spread. But by then we had no scorer we could trust to tell
us whether it helped.

And treat the close of validation as the end. Everything after it was untestable.

## Where I think this is going

With public data, predicting a new cell line is mostly a data problem. It is not a modelling problem.
[Ahlmann-Eltze, Huber and Anders](https://doi.org/10.1038/s41592-025-02772-6) showed deep models do not yet beat
simple linear baselines on perturbation prediction. [Arc's own wrap-up](https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up)
said much the same. The signal that helps is in-context perturbation data. The public does not have enough of it across enough cell types. The leader has a private
[atlas of more than 200 cell lines](https://www.illumina.com/company/news-center/press-releases/2026/fda84c92-b4b3-4691-a402-35555abe8605.html).

I could be wrong, and the final scores should show it. If a top entry used no same-cell-type data and still beat
0.3, this is a modelling problem and I had the bottleneck wrong.

## What would make it easier

- The same knockdowns measured across many more cell lines.
- Matched controls and replicates, so everyone knows the noise ceiling per target.
- Measured knockdown strength per target.
- A public held-out panel that ranks entries the way the official board does.

## Still open

Was the leaders' edge same-cell-type data? How should response size transfer into a cell type where you only have
controls? Which public screens add independent signal?

## Thanks

Arc Institute for the challenge and [cell-eval2](https://github.com/ArcInstitute/cell-eval2). kaipengm2 for
AtlasShift, our first real baseline. Saber HQ for [Sidechain](https://github.com/saberhq/sidechain). Bill M, whose
described response network became our best entry. Matin Gerami for
[zero-shot-perturbation](https://github.com/MMatinGerami/zero-shot-perturbation),
[danielchen26](https://github.com/danielchen26/vcc2026-stage1-lab) and the Rasmussen Lab for
[SLIM](https://github.com/RasmussenLab/SLIM).
