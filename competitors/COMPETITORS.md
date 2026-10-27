# VCC 2026 competitor methods sweep (2026-10-17)

Top 50 of the validation leaderboard, from the snapshot `leaderboard-2026-10-17.txt`. Live data comes from the public endpoint
`https://virtualcellchallenge.org/api/leaderboard?get_final=false`, which needs no login and returns organisation, model
name, a free-text **description**, and all six raw and scaled components for 1,564 entries. The local copy is `raw/lb-2026-10-17.json`,
every non-empty description is in `raw/leaderboard-descriptions-2026-10-17.txt`, and raw member rosters are not exported.

Search method: four parallel sweeps covering ranks 1-11, 12-31 and 32-50, plus public repos, the 2025 winners and method families.
Each team was searched by its public identifiers on the web, GitHub repos, GitHub code, bioRxiv/arXiv/Europe PMC, Hugging Face,
and public research sites. X and LinkedIn need a login and were not read. Zhihu returned 403. GitHub code search was rate-limited for
part of the run. Per-team submission history (`/api/leaderboard/team/<id>`) requires authentication, so it was not read. A "nothing
found" below means nothing turned up in those sources; it does not mean nothing exists.

Confidence levels:
- **high**: the team's own write-up or code.
- **med**: strong circumstantial evidence, such as the team's own leaderboard description or a reconstruction by another team.
- **low**: a guess from a name or a score signature.
- **none**: nothing found.

## Headline findings

1. **No top-50 team has published code or a write-up (0/50 high confidence).**
   - 7/50 have method-level evidence at medium or low-medium confidence: #1 Illumina, #13 pop, #16 YaochunL, #17 C-team, #31 William, #36 isee and #37 11111.
   - 1 more has a weak hint (#25 null_hypothesis).
   - Other leads identified an organisation without establishing a method.
   - The remaining ~33 have nothing. em2.bio (#35) says "Full recipe at em2.bio", but the page is still "Coming soon".
2. **The 0.30 tier wins mostly on MSE and reach, not PDS.**

   Median scaled components:

   | Group | pds | mse | nmae | fid | reach | jac |
   |---|---:|---:|---:|---:|---:|---:|
   | Ranks 11-50 | .807 | .442 | .238 | .043 | .258 | .011 |
   | AtlasShift copy | .682 | .047 | .137 | -.017 | .080 | -.001 |
   | Our v4a | .678 | **0** | .129 | -.019 | .076 | -.001 |

   About half of the gap from v4a (0.144) to 0.30 is MSE: +0.44 scaled, about +0.07 overall.
3. **The MSE gap is mostly the scorer's sampling-correction budget, and it is label-free.** cell-eval2's `expr_mse_unbiased_capped`
   subtracts a sampling correction: `r * min(jk_pred, jk_real)`, where `r = min(1, B_pred/claim)`. `B_pred` is the submission's own
   across-perturbation centred sum of squares (pseudobulk spread).
   - A submission that **pins each perturbation's per-gene totals** to its predicted pseudobulk, keeps realistic cell-to-cell
     dispersion, and has B/claim >= 1 collects the full correction.
   - Arc measured that pinning totals on an honest control-paste arm moved scaled MSE from **0.000 to 0.903**.
   - Arc documents the non-binding regime (B >= claim) as a residual it chose not to close: cell-eval2 `CHANGELOG.md`, issue #348
     section, and the `src/cell_eval2/metrics/delta.py` constants at commit a2630dc, verified locally.
   - Top teams' descriptions say they tune exactly this. Bill M #127: "across-target spread calibrated to B/claim 1.00", MSE .349.
     Hosein Fooladi #228: "0.97x the MSE sampling claim". binbin #261: "refund-ratio target rho 0.95". The reconstruction of William's #31
     recipe uses "noise/claim=0.99". Jason Jowu #328 gives a label-free rule: "predicted across-target pseudobulk spread = 0.95 x
     claim".
   - The rule depends only on our own prediction and the public controls, so it **carries to the Oct 22 new-cell-line test**.
4. **Most top-50 entries come in tight score clusters**, with near-identical six-component vectors submitted within hours of each
   other. These are observed score similarities; they do not establish shared code, common ownership or a public recipe:
   - **A**: #23 SQ W, #26 Ava Brown, #28 Anna Xie, #29 HY H, #30 Vivian Xu, all at pds .831 / mse .442 / nmae .220.
   - **B**: #21 BioFarmer and #24 Jianhua Wang, with #32 ITHACA ZERO and #35 em2.bio nearby.
   - **C**: #14 CYH and #19 nailong king, with #12 WCK nearby.
   - **D**: #13 pop, #16 YaochunL and #20 CellWM, the "M014" lineage.
   - **E**: #33 TOBE, #34 Layne Wong, #40 Lei Huang, #42 Lingkuan Meng and #46 黄龙. Jaccard is exactly .0115 for four of them.
   - **F**: #41 Patheia and #47 LiuLab are identical, with #48 Gru a near-twin.
   - **G**: #36 isee and #37 11111.
5. **Illumina's lead (0.467) comes from magnitude accuracy, not discrimination.**

   | Raw component | Illumina | Ranks 2-11 |
   |---|---:|---:|
   | MSE | 0.235 | 0.46-0.57 |
   | NMAE | 0.741 | 0.83-0.87 |
   | Reach | 0.483 | 0.31-0.41 |
   | PDS cosine | 0.898 | 0.85-0.88 |

   The most likely explanation is the proprietary **Illumina Billion Cell Atlas**: genome-wide CRISPR Perturb-seq in more than 200
   cell lines, explicitly "to train virtual cell models". The rules allow private data. We cannot reproduce that data advantage.
6. **Leaderboard probing is common** (85-112 submissions; per-context weight probes; "order-only probe"). It will not carry to a new
   300-target, new-cell-line panel. Validation ranks overstate some teams' generalization.

## Per-team table

Public team/model names and scores are from the snapshot plus the public API; individual member rosters and identity guesses are omitted. "Live" is the live rank at about 10:00 EDT, which moves within hours.

| # | Team / org / model | Score | What was found (links) | Method (1-3 lines) | Conf. |
|---|---|---:|---|---|---|
| 1 | Illumina AI / Illumina / iCell  | 0.4670 | Public Illumina Billion Cell Atlas announcements: https://www.illumina.com/areas-of-interest/genomics-in-drug-development/cell-atlas.html ; no executable scored recipe established. | Model trained on measured same-gene knockdown responses across 200+ proprietary lines (some likely close to the hidden lines), conditioned on target-line controls, with 3'-to-Flex chemistry calibration. The lead is effect magnitude (MSE/NMAE/reach) | med (data) / low (arch) |
| 2 | X-team / ABC / X-91  | 0.3381 | No verified public implementation linked to this leaderboard entry. | Method unconfirmed. | none |
| 3 | Hhx Hhx / Myrtle  | 0.3328 | Withdrawn from the live board. Nothing found | - | none |
| 4 | cqawesome / m0998t  | 0.3312 | No verified public implementation found; 112 submissions in this snapshot. | - | none |
| 5 | realcell / test  | 0.3274 | Nothing found (GitHub user realcell is empty) | - | none |
| 6 | JJJ / GG-33  | 0.3263 | Nothing found | - | none |
| 7 | TAO / UT Austin / v67  | 0.3232 | Nothing found | - | none |
| 8 | hya / pp  | 0.3180 | Nothing found | - | none |
| 9 | roxs / halcyon43  | 0.3146 | Nothing found | - | none |
| 10 | Testh / test  | 0.3134 | Nothing found. Top 3 in late Sept / early Oct | - | none |
| 11 | Aginglab.com / AgingLab / VirtualOmics  | 0.3131 | Public team organisation: AgingLab. No verified VCC method found; later dropped off the live board. | - | none (method) |
| 12 | WCK / Rice / genoguessr-v16  | 0.3105 | Nothing found ("genoguessr" code and web search empty). Score cluster C | - | none |
| 13 | pop / M014-q10  | 0.3093 | Own description: "M014-q8 with the delta slot from B171 (B157 recipe + MSE weight x100)". No code (M014/B157/psig code search empty) | Slot-based ensemble of recipes. One "delta slot" comes from a recipe trained with a 100x MSE loss weight. Same lineage as #16 | low-med |
| 14 | CYH / Westlake / Cyan  | 0.3090 | Nothing found. Cluster C twin of #19 | - | none |
| 15 | FlowState / flowstate  | 0.3082 | No verified VCC implementation found. | Method unconfirmed. | none |
| 16 | YaochunL / M014-q8g  | 0.3080 | Own description: "M014-q8 with a weaker psig delta shrink (gamma 1.5 -> 0.75) to test the pds lever". Submitted 15 min after #13 | Same M014 code. The predicted delta is shrunk by a per-target significance score ("psig"), trading MSE against PDS | low-med |
| 17 | C-team / NNF CBMR, Univ. Copenhagen / The_Way45  | 0.3076 | Public SLIM method: https://github.com/RasmussenLab/SLIM (MIT); no scored implementation established for this team. | SLIM is closed-form bilinear ridge on gene embeddings; connection to this scored recipe is unverified. | med (tools) / low (recipe) |
| 18 | chill / Excuse  | 0.3074 | Nothing found. Highest Jaccard in the band (.035) | - | none |
| 19 | nailong king / nailong  | 0.3071 | Nothing found. Cluster C | - | none |
| 20 | CellWM / CellWM-94  | 0.3053 | Nothing found. Scores next to the M014 pair | - | none |
| 21 | BioFarmer / SJTU_Agri  | 0.3047 | Nothing found. Cluster B | - | none |
| 22 | Hopkin Jiang / bbbbbbbbulk  | 0.3045 | Nothing found | Name suggests pseudobulk-centred | low |
| 23 | SQ W / sq-1017b  | 0.3037 | Nothing found. Cluster A | - | none |
| 24 | Jianhua Wang / Milk ball  | 0.3036 | Nothing found. Cluster B | - | none |
| 25 | null_hypothesis / e62  | 0.3032 | Public CD4 buffering analysis: https://github.com/urwahnawaz/cd4_buffering . Team implementation not established. | Possibly a responder/non-responder gate built from paralog, degree and constraint features (affects reach, fid, jac). Unconfirmed | low |
| 26 | Ava Brown / cellpred-0.9.64  | 0.3032 | Nothing found. Cluster A | - | none |
| 27 | Junru Jin / dld05skip  | 0.3027 | No verified public VCC implementation found. | - | none |
| 28 | Anna Xie / M69  | 0.3027 | Nothing found. Cluster A | - | none |
| 29 | HY H / Insight / hh-v74  | 0.3024 | Nothing found. Cluster A | - | none |
| 30 | Vivian Xu / run_2026_68  | 0.3013 | Nothing found. Cluster A, loosely | - | none |
| 31 | William / transfer-protocol-v50c  | 0.3008 | Recipe reconstructed in **another team's** description (gugugaga #560, "Public-information reconstruction of transfer-protocol v41a", which scored 0.140 because three sources were missing) | Sources: K562 GWPS, Orion HCT116 and HEK293T, H1 CRISPRi, CD4 T CRISPRi, HipSci iPSC, KOLF2.1J, VIPerturb-seq, LINCS L1000. Source-centred, target-specific transfer in two branches (DE; depth-weighted log-of-mean pseudobulk), adapted to each context from its controls, per-context bulk amplitude at **noise/claim = 0.99**, dual-moment raw counts | med |
| 32 | ITHACA ZERO / SJTU_AB / Viking_Shiiiip  | 0.3007 | Nothing found. Best reach in its group (.271) | - | none |
| 33 | TOBE / 261017_v2  | 0.3003 | Nothing found. Cluster E | - | none |
| 34 | Layne Wong / C20  | 0.3001 | Nothing found. Cluster E (C18/C20 numbering is shared with #40) | - | none |
| 35 | em2.bio / THATWASWEIRD  | 0.2992 | Public description says full recipe at https://em2.bio ; VCC project still marked Coming soon in this historical scan. | Unpublished | none |
| 36 | isee / Larch-17  | 0.2977 | Own description: "cleaner borrowed responses". Cluster G | Atlas-style transfer of public responses, with source denoising | low-med |
| 37 | 11111 / v1023  | 0.2976 | Own description: "atlas". Cluster G twin of #36 | Atlas transfer (AtlasShift family) far past 0.155, presumably with the emission and calibration levers above | low-med |
| 38 | Frankbest18 / Candidate-027  | 0.2974 | No verified public VCC implementation found. | Method unconfirmed. | none |
| 39 | TJU-SynBio / Tianjin Univ / v38  | 0.2962 | Nothing found. Best MSE in its group (.412) and reach (.287) | - | none |
| 40 | Lei Huang / C18  | 0.2962 | Nothing found. Cluster E | - | none |
| 41 | Patheia / Eonark / Cell Papa  | 0.2962 | Organisation: https://eonark.ai (Shenzhen longevity virtual-cell startup, $14M seed https://thesaasnews.com/news/eonark-raises-14m-seed). Cluster F, identical to #47 | - | none (method) |
| 42 | Lingkuan Meng / test  | 0.2960 | Nothing found. Cluster E | - | none |
| 43 | Yifan SONG / E_primary_v82  | 0.2958 | Nothing found | - | none |
| 44 | Advaita Labs / Advaita_V9  | 0.2956 | https://advaita.xyz ("Aegis": auto research + virtual cell), https://github.com/AdvaitaLabs (no VCC code) | Probably produced by their autoresearch system | low |
| 45 | TerminatorJ / kkk  | 0.2953 | No verified public VCC implementation found. | Method unconfirmed. | none |
| 46 | 黄龙 / toy model  | 0.2953 | Nothing found. Cluster E | - | none |
| 47 | LiuLab / Larkspur  | 0.2952 | Nothing confirmed. github.com/Liuqi3372/vcc2026 (Chinese AtlasShift repackage at 0.156) is not linked to it. Cluster F | - | none |
| 48 | Gru / Wintergreen  | 0.2932 | Nothing found. Cluster F near-twin | - | none |
| 49 | Q Zhao / P2  | 0.2920 | Nothing found | - | none |
| 50 | Lei Wang / simula-v11  | 0.2912 | Nothing found | - | none |

## Recipes disclosed by other teams' leaderboard descriptions (ranks 50-600)

These are the competitors' own words, from `raw/leaderboard-descriptions-2026-10-17.txt`. Each item is technique: who used it.

- **Pooled amplitude calibrated to the MSE claim**:
  - Bill M #127 (0.248, MSE .349): "B/claim 1.00".
  - Hosein Fooladi #228: "0.97x claim".
  - binbin #261: "refund ratio rho 0.95; adding a 5th source changed spread and undershot".
  - Jason Jowu #328: label-free "spread = 0.95 x claim", amplitudes A/B/C 0.515/0.543/0.546.
  - Mengran Zhang #383: "auto B/claim >= 1".
- **Pinned group totals / dual-moment emission**:
  - Jae #339: "pseudobulk totals constrained to expected values, generated cells are dependent ... the iid-based noise deduction can improve the reported MSE without changing the predicted mean".
  - RAISE #529: "exact integer group gene totals".
  - Zhu Lab #1369: "exact gene totals round(400*mean)".
  - qikun zhao #169: "dual-moment fitting".
- **Control-axis deflation** (project the top-k control PCs of the context's own NTC cells out of the predicted bulk response):
  - Jason Jowu #328, Hao Tian #376 and Hualuo Pang #381 (top-50/100, eta 1).
  - Becker Lab variant: 0.1669 vs plain AtlasShift 0.1546.
  - Chaos #273 (10 PCs, eta 1).
- **Explicit DE call selection**:
  - Chen Zhang #82 (0.273): PIE p_de rank to choose padding DE calls.
  - SIVC #167: "select 300 non-self DE calls".
  - eye-contact #96: Wilcoxon-constrained count exchanges.
  - Suyash Dhoot #358: 3000-4500 DE calls per context.
- **Context identity**: context A is "T-lineage / Jurkat by DepMap identity" (ASK SATHVIK #356, Sulgi Kim #357, FoRJu #362). These teams weight CD4 T-cell or Nadig Jurkat sources up for A.
- **Extra sources**: CD4 primary T (Zhu 2026, 10x Flex, the same chemistry as the hidden data), HipSci iPSC (Feng), KOLF2.1J, VIPerturb-seq, LINCS L1000, Pan hESC GSE295214, Tian 2021 CRISPRi, Nadig Jurkat/HepG2. One proprietary set, "ELEM-1200RS" from ELEM Biotech, is used by entries scoring about 0.22 and was not findable.
- **Sign / direction models**: Kaipeng Zheng #288 "Atlas-SymmetricRank-v4" (0.199) trains a binary sign model on K562 GWPS directions with PIE/atlas features; it changes only the LFC direction confidence.

## Public code (all below 0.20 claimed)

| Repo | Claimed official | Method | License |
|---|---:|---|---|
| https://github.com/saberhq/sidechain | 0.1684 (SER-19) | Multi-source transfer + post-pooling adaptive shrinkage; full ablation log at https://saberhq.com/sidechain/models/ | MIT |
| https://github.com/kaipengm2/Virtual-Cell-Challenge-2026 (AtlasShift d24ce4f) | 0.1546 | 5-source weighted transfer, promoter correction, dual-summary emission | MIT |
| https://github.com/ferry-ai/virtual-cell-challenge-2026 | 0.149 | CRISPRi transfer + KO votes | none |
| https://github.com/MMatinGerami/zero-shot-perturbation | 0.137 | Equal-weight, norm-restored transfer across 8 screens (+0.038 vs 6 similarity-weighted sources); noise-corrected cross-line correlation of the same knockdown is about 0.47 | MIT |
| https://github.com/danielchen26/vcc2026-stage1-lab , https://github.com/danielchen26/vcc2026-dossier | n/a | Scorer analysis: the DE metrics depend only on the predicted per-gene mean and the mean rank of predicted cells against the published controls. Fast psi-rank DE (41x on an M1 Pro, 16 GB); target-of-200-375 DE genes estimate | none (read-only) |
| https://github.com/forrestsheldon/vcc2026-h1-benchmark | 0.098 | Local H1 harness | MIT |
| https://github.com/RasmussenLab/SLIM | n/a | Closed-form bilinear ridge on gene embeddings (C-team's lab) | MIT |
| https://github.com/ArcInstitute/pie | n/a | DE-gene probability in unseen contexts (1.2-3.2x AUPRC over baselines) | see repo |

## 2025 winners (prior art)

Arc wrap-up: https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up

- **#1 BioMap xTrimoSCPerturb.** scFoundation plus PLM target embeddings, trained on pseudobulk with a PDS+DES loss and DEG-frequency / mean-expression features. No code.
- **#2 XLearning (Sichuan).** MLP on mean control expression, an ESM-2 target embedding and depth, predicting the delta.
- **#3 Outlier TransPert.** Purely statistical cross-line pseudobulk/Wilcoxon transfer, similarity-weighted, predicts which genes are significant, one global scale. Related papers:
  - https://arxiv.org/abs/2511.16954: direction transfers across lines, magnitude does not.
  - https://www.biorxiv.org/content/10.64898/2026.10.05.756839v1: an essentiality-tied global response plus a transferable specific part.
- **Generalist prize, Altos PRiMeFlow.** https://arxiv.org/abs/2604.13986 , https://github.com/altoslabs/primeflow : flow matching pretrained on 14 lines. GPU-heavy.
- **Turbine.** Ridge on pseudobulk with a cells-per-perturbation feature: https://blog.turbine.ai/p/how-did-we-get-a-regression-model

## Evidence on unseen-cell-line generalization (the Oct 22 test)

- **Per-knockdown mean across source lines beats deep models on an unseen line.** AdaPert Table 13 (https://arxiv.org/abs/2602.18885), train RPE1/HepG2/Jurkat, test K562, Pearson-delta:

  | Method | Pearson-delta |
  |---|---:|
  | Perturb Mean | 0.650 |
  | Linear on GenePT | 0.515 |
  | AdaPert | 0.460 |
  | TxPert | 0.405 |
  | STATE | 0.309 |

- **The response decomposes** (Molina & Zhang, https://www.biorxiv.org/content/10.64898/2026.07.24.740459v1): line template 28%, shared effect 29%, line-by-knockdown interaction 24%, noise 19%. No method recovers the interaction part, and every extra source line helps.
- **Linear and embedding models transfer between lines.** Ahlmann-Eltze et al., Nat Methods 2025: https://pmc.ncbi.nlm.nih.gov/articles/PMC12328236/ .
- **Generic models struggle.** STATE matches only the perturbation mean on genetic data, and a public STATE fine-tune (CaiZongyuan) scored below control cells.
- **The levers that should transfer** are the label-free ones: scorer-aware emission and calibration, wider equal-weight source banks, control-axis deflation, PIE DE ranking, and responder gating. The ones that should not are per-context weights probed on the leaderboard and target-specific tuning on validation.

## Synthesis

### Techniques that recur among the top teams (inferred from descriptions, clusters and component profiles)

1. **Scorer-aware count emission**: exact per-(target, gene) totals, dual-moment cells, per-context amplitude tuned to B/claim of about 0.95-1.0. This is the MSE lever, worth about +0.4 scaled in the top 50 versus about 0 for us.
2. **Multi-source same-target transfer**, with 6-9 public screens, source centring, and separate DE and pseudobulk branches.
3. **Explicit DE call-set control**: choosing which genes are called (about 250-300 per target), with what sign and magnitude, using PIE or cross-source agreement. This is the reach/fid lever.
4. **Context-identity routing**, e.g. A = T-lineage, so CD4 or Jurkat sources are weighted up.
5. **Control-axis deflation**.
6. **Shrinkage of the delta by significance** (M014 "psig").
7. **Heavy leaderboard iteration** (50-112 submissions).

### Reproducible in under 2 days on the 16 GB M4 (CPU/MPS)

- **(a) B/claim amplitude calibration plus pinned totals.** Uses cell-eval2's own `correction_for` / `delta.py` functions; CPU. 0.5-1 day.
- **(b) Control-axis deflation.** A few lines of numpy. Hours.
- **(c) Sidechain reproduction.** MIT, CPU. Already partly done.
- **(d) SLIM ridge for coverage gaps.** MIT, seconds.
- **(e) PIE p_de ranking for DE calls.** Check the model size first.
- **(f) psi-rank fast DE harness.** Re-implement the idea; the danielchen26 repo has no license.
- **(g) Extra sources.**
  - HipSci, KOLF2.1J and VIPerturb-seq need new downloads; check the 15 GB budget and the per-reader ceilings in `rsi/REPORT.md`.
  - LINCS L1000 is a different modality and lower priority.

### Ranked "adopt next" list

Expected gains are official-scale, overall, and UNVERIFIED until scored locally.

| Rank | Adopt | Expected gain | Evidence | Generalizes to new lines? |
|---|---|---|---|---|
| 1 | **Scorer-aware emission.** Pin per- integer totals to the predicted pseudobulk), keep within-target cell dispersion at about reference level, and choose per-context pooled amplitude label-free so the predicted across-target centred SS is B = 0.95-1.0 x claim  | MSE 0 -> 0.3-0.45 scaled, so **+0.05 to +0.075** | Public Illumina Billion Cell Atlas announcements: https://www.illumina.com/areas-of-interest/genomics-in-drug-development/cell-atlas.html ; no executable scored recipe established. | Yes. Uses only our prediction and the public controls |
| 2 | **Explicit DE call-set control.** About 250-300 calls per target ranked by cross-source agreement + PIE p_de; consensus sign; make called genes Wilcoxon-detectable against the published controls; leave the rest uncalled  | reach +0.1-0.18, fid +0.05, nmae +0.05, so **+0.03 to +0.06** | No verified public implementation linked to this leaderboard entry. | Method unconfirmed. | none
| 3 | **Wider equal-weight source bank, William-style.** Add CD4, HipSci iPSC, KOLF2.1J and VIPerturb-seq to K562/Orion/H1; source-centre; norm-restore; separate DE and depth-weighted log-mean pseudobulk branches  | pds +0.05-0.10, so **+0.01 to +0.03** (plus it feeds 1 and 2) | William #31 at 0.3008 vs the 0.140 reconstruction missing 3 sources; AdaPert Table 13 (Perturb Mean best on an unseen line); Molina & Zhang (each extra line helps); MMatinGerami +0.038 | **Yes, the strongest biology lever for unseen lines** |
| 4 | **Control-axis deflation**: project out the top 50 PCs of each context's own NTC cells from the predicted bulk response  | **+0.01** (pds .768 vs .682) | No verified public implementation found; 112 submissions in this snapshot. | Yes. Label-free |
| 5 | **Responder gating + coverage fill**: predict silent knockdowns and shrink them toward zero; SLIM ridge on STRING embeddings for targets with no source measurement; M014-style significance shrink of the delta  | **+0.005 to +0.02** | null_hypothesis cd4_buffering (p<1e-14 predictors); C-team / SLIM; YaochunL "psig delta shrink"; essentiality-tied global response preprint | Yes |

**Do not** spend rounds on:
- per-context weights tuned on validation probes (they won't carry to the Oct 22 lines);
- STATE fine-tunes (below control cells);
- deep generative emitters (Sidechain PHE-2 scored -0.98);
- control-cell shift emitters (Sidechain "r" collapsed fid);
- another global amplitude sweep without the B/claim rule.

A realistic target from 1+2+4 on the v4a base is about 0.24-0.28; that is an estimate, not a measurement. Illumina's 0.467 needs proprietary
multi-line data we don't have.

## Re-check list before Oct 22

- em2.bio recipe page.
- The leaderboard descriptions: re-pull the API and diff, since new descriptions leak recipes.
- X-Atlas/Pisces data release (HepG2/iPSC/Jurkat).
- Illumina bioRxiv / blog posts on iCell.
- New commits to Sidechain or AtlasShift.
