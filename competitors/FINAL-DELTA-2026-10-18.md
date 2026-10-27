# Leaderboard changes after the initial audit — 2026-10-18

Read-only refresh 2026-10-18T17:07:19.974091+00:00: 1577 entries, 604 descriptions. 105 new IDs/changed descriptions versus the initial1,569-entry snapshot; 40 described changes below. All other current descriptions are covered by DESCRIPTION-AUDIT-2026-10-18.md. Historical removed entries remain audited there. The initial whole-board cluster tables and medians remain explicitly pinned to that earlier snapshot.

All component vectors below are scaled PDS / MSE / NMAE / fidelity / reach / Jaccard. Participant descriptions are claims, not verified causal ablations. No new candidate was chosen using these vectors or the failed proxy. Full records and descriptions: raw/final-leaderboard-delta.json.

## #33 11111 — 0.300320
Model: `v1026`. Vector: 0.8185 / 0.4329 / 0.2495 / 0.0416 / 0.2436 / 0.0158.
Quote: “registry…”
Concrete lever: No concrete reproducible rule beyond model/status name.
Component: Unidentified from description.
Label-free / final transfer: Unknown.
Expected gain / cost: No defensible gain estimate or implementation cost.

## #86 tenglong yuan — 0.268613
Model: `mamahaha-v195-linear-scgrn-tftg-collectri-a010-20261018`. Vector: 0.8004 / 0.3921 / 0.1522 / 0.0362 / 0.2317 / -0.0008.
Quote: “V195 baseline + one-hop linear GRN correction alpha=0.1 before PCA. CollecTRI signed edges plus scGRN blood T-cell and…”
Concrete lever: One-hop linear GRN residual alpha=.1 before PCA; CollecTRI signed edges, scGRN blood-T and TFTG CD4/HepG2 promoter evidence; context control correlations weight edges; GRAND unavailable403.
Component: Direction / network prior; V195 parent not released, no isolated contribution established.
Label-free / final transfer: Potentially label-free via each context controls; fixed public graph choices need validation.
Expected gain / cost: Unidentified; hours for graph integration, original parent unavailable.

## #138 HANYOUNG JEONG — 0.246824
Model: `VC-Nova-v33`. Vector: 0.7787 / 0.2974 / 0.1534 / 0.0297 / 0.2191 / 0.0026.
Quote: “VC-Nova-v33…”
Concrete lever: No concrete reproducible rule beyond model/status name.
Component: Unidentified from description.
Label-free / final transfer: Unknown.
Expected gain / cost: No defensible gain estimate or implementation cost.

## #203 Lu Zhao — 0.223746
Model: `86a4b3c7e29546`. Vector: 0.7045 / 0.2292 / 0.1820 / 0.0074 / 0.2177 / 0.0017.
Quote: “This work used the ELEM Biotech Virtual Cell Perturbation Transcriptomics Dataset, Release [ELEM-1200RS V1.0], provided under the ELEM…”
Concrete lever: ELEM priority-access data, no emission recipe.
Component: Biological source coverage.
Label-free / final transfer: Unavailable proprietary/priority-access asset.
Expected gain / cost: Unknown; unavailable data.

## #210 Jiabao Zhu — 0.220962
Model: `cee4247cee9032`. Vector: 0.7010 / 0.2145 / 0.1818 / 0.0078 / 0.2190 / 0.0017.
Quote: “This work used the ELEM Biotech Virtual Cell Perturbation Transcriptomics Dataset, Release [ELEM-1200RS V1.0], provided under the ELEM…”
Concrete lever: ELEM priority-access data, no emission recipe.
Component: Biological source coverage.
Label-free / final transfer: Unavailable proprietary/priority-access asset.
Expected gain / cost: Unknown; unavailable data.

## #237 Javier Villela Castrejon — 0.212889
Model: `JaviCell`. Vector: 0.7434 / 0.1529 / 0.1561 / 0.0005 / 0.2173 / 0.0071.
Quote: “JaviCell replicate…”
Concrete lever: No concrete reproducible rule beyond model/status name.
Component: Unidentified from description.
Label-free / final transfer: Unknown.
Expected gain / cost: No defensible gain estimate or implementation cost.

## #317 xxx yan — 0.190474
Model: `restore best artefact (amp13_pk2)`. Vector: 0.7500 / 0.1033 / 0.1408 / -0.0091 / 0.1562 / 0.0017.
Quote: “Restoring our best-scoring artefact. Content is bit-identical to the entry that scored 0.19047403659226078 (sha256 0b727e0095e3022838eaf0f638b2f2aa6385cfacfb1aace3dc4924f…”
Concrete lever: Restored pool2 artefact; author reports pool2→1 reach+.046449 but fidelity−.092812 and Jaccard−.014609; manual A-source-weight probe also negative.
Component: DE dispersion tradeoff; author ablation, not independently reproduced.
Label-free / final transfer: Pooling rule yes; hand-picked A weights excluded.
Expected gain / cost: Negative all-gene-scatter evidence; emitter hours.

## #327 Haibaio Yao — 0.188360
Model: `atlas-p2-elem-bulkonly-v1`. Vector: 0.7632 / 0.1010 / 0.1000 / 0.0006 / 0.1600 / 0.0054.
Quote: “Experimental canonical P2 + ELEM bulk-only0.5, mean0; fixed five-source BASE recipes. 5 direct targets; H132 MSE+0.519%, NMAE slightly…”
Concrete lever: ELEM bulk-only .5, mean0 on five-source P2; only5 directly covered targets; reported H132 MSE+.519% crosses zero.
Component: Bulk direction; small/uncertain source addition.
Label-free / final transfer: Priority-access ELEM unavailable.
Expected gain / cost: Unidentified; unavailable data.

## #338 ants cell — 0.185438
Model: `p67 zero callhead`. Vector: 0.6588 / 0.1613 / 0.0959 / 0.0210 / 0.1759 / -0.0002.
Quote: “Restore the published P67 zero-callhead package. No fold change.…”
Concrete lever: Restore zero-callhead P67; no fold-change modification specified.
Component: Bulk-only output can outperform some DE heads; causal gain unknown.
Label-free / final transfer: No executable zero-callhead specification.
Expected gain / cost: Unknown; insufficient recipe.

## #349 Jiaxin Han — 0.182008
Model: `HJX-noise-credit-reliability-v1`. Vector: 0.6297 / 0.1830 / 0.1466 / -0.0051 / 0.1465 / -0.0088.
Quote: “Noise credit calibrated transfer with gene level transfer reliability. Same output design as v2, with each gene's centred…”
Concrete lever: Gene effect reliability learned from3 public proxy cell lines, clipped .5–1.5; same noise-credit renderer.
Component: PDS/NMAE and DE ranking; no isolated official attribution.
Label-free / final transfer: External-source-only training can transfer; proxy validity must be re-established.
Expected gain / cost: Unknown; days for model/data work.

## #361 qin han — 0.179209
Model: `kd70+cv35 (bulk-exact, CV=35)`. Vector: 0.6838 / 0.1113 / 0.1410 / -0.0036 / 0.1290 / 0.0138.
Quote: “Exact group-bulk-preserving dispersion sweep. Same kd70 base, per-cell gamma CV varies, integerisation preserves every group column sum exactly…”
Concrete lever: Gamma dispersion CV35; every group column sum preserved; kd70 parent.
Component: DE ranks and JK claim; fixed bulk guarantees fixed PDS, NOT universally fixed MSE because JK changes.
Label-free / final transfer: Universal CV/column-total rule can transfer.
Expected gain / cost: Unknown; hours; no causal improvement over our model shown.

## #362 Yixuan Hu — 0.179039
Model: `kd70+cv25 (bulk-exact, CV=0.25)`. Vector: 0.6838 / 0.1113 / 0.1437 / -0.0020 / 0.1270 / 0.0104.
Quote: “Exact group-bulk-preserving dispersion sweep, CV=0.25 -> ~4300 DE calls, targeting the peak found between cv35 (3044 calls, DE-block…”
Concrete lever: Exact group-column-total gamma dispersion CV=.25; author targets about4300 DE calls between cv35/3044 calls and anchor/6243 calls.
Component: DE call-count tradeoff; fixed bulk fixes PDS but does not generally fix JK-adjusted MSE.
Label-free / final transfer: Universal dispersion rule can transfer; hidden true-call optimum unavailable.
Expected gain / cost: Unidentified gain for us; emitter hours; .17904 is not a verified improvement over our renderer.

## #365 VCC@ZJU — 0.177859
Model: `Joyee-69`. Vector: 0.6939 / 0.1011 / 0.1460 / -0.0115 / 0.1359 / 0.0019.
Quote: “AtlasShift Joyee-69; 按预批准规则择优精确拼接ABC；重复已有配置则仅B bulk=0.55。 raw2:1:1:2/CD4/GTF; native-first; float L1<=0.001; projection<=0.03; seed20260910; 300xABCx400x18533.…”
Concrete lever: Explicit context-specific mean/bulk/CD4 coefficients, raw source2:1:1:2, native-first, seed20260910, projection≤.03 and floatL1≤.001; full recipe in description.
Component: Independent moments / context coefficients.
Label-free / final transfer: NO as written: ABC splicing or manually specified A/B/C coefficients excluded.
Expected gain / cost: Unknown; hours; transferable part is independent moments, already tested.

## #384 Yiyang Li — 0.174152
Model: `Joyee-62`. Vector: 0.7051 / 0.0850 / 0.1209 / -0.0114 / 0.1427 / 0.0026.
Quote: “AtlasShift Joyee-62; 相对26仅B：mean=1.5、bulk=0.5、CD4 mean/bulk=0.5/2；其他两context counts不变。 raw2:1:1:2/CD4/GTF; native-first; float L1<=0.001; projection<=0.03; seed20260910; …”
Concrete lever: Explicit context-specific mean/bulk/CD4 coefficients, raw source2:1:1:2, native-first, seed20260910, projection≤.03 and floatL1≤.001; full recipe in description.
Component: Independent moments / context coefficients.
Label-free / final transfer: NO as written: ABC splicing or manually specified A/B/C coefficients excluded.
Expected gain / cost: Unknown; hours; transferable part is independent moments, already tested.

## #385 lvrongcen lv — 0.174100
Model: `kd70+cv75 (bulk-exact, CV=75)`. Vector: 0.6838 / 0.1113 / 0.1332 / -0.0462 / 0.1436 / 0.0188.
Quote: “Exact group-bulk-preserving dispersion sweep. Same kd70 base, per-cell gamma CV varies, integerisation preserves every group column sum exactly…”
Concrete lever: Gamma dispersion CV75 with the same described kd70 parent and exact group column totals.
Component: DE ranks and JK claim; MSE-invariance claim is not general.
Label-free / final transfer: Universal rule can transfer.
Expected gain / cost: Unknown; hours; worse total than CV35 observationally.

## #386 Ruan Jiafa — 0.174040
Model: `Joyee-60`. Vector: 0.7001 / 0.0792 / 0.1328 / -0.0121 / 0.1419 / 0.0024.
Quote: “AtlasShift Joyee-60; 相对26仅B：mean=1.25、bulk=0.5、CD4 mean/bulk=0.5/1；其他两context counts不变。 raw2:1:1:2/CD4/GTF; native-first; float L1<=0.001; projection<=0.03; seed20260910;…”
Concrete lever: Explicit context-specific mean/bulk/CD4 coefficients, raw source2:1:1:2, native-first, seed20260910, projection≤.03 and floatL1≤.001; full recipe in description.
Component: Independent moments / context coefficients.
Label-free / final transfer: NO as written: ABC splicing or manually specified A/B/C coefficients excluded.
Expected gain / cost: Unknown; hours; transferable part is independent moments, already tested.

## #392 gogogo go — 0.173256
Model: `Joyee-58`. Vector: 0.7045 / 0.0821 / 0.1214 / -0.0122 / 0.1411 / 0.0027.
Quote: “AtlasShift Joyee-58; 相对26仅B：mean=1.5、bulk=0.5、CD4 mean/bulk=0.5/1；其他两context counts不变。 raw2:1:1:2/CD4/GTF; native-first; float L1<=0.001; projection<=0.03; seed20260910; …”
Concrete lever: Explicit context-specific mean/bulk/CD4 coefficients, raw source2:1:1:2, native-first, seed20260910, projection≤.03 and floatL1≤.001; full recipe in description.
Component: Independent moments / context coefficients.
Label-free / final transfer: NO as written: ABC splicing or manually specified A/B/C coefficients excluded.
Expected gain / cost: Unknown; hours; transferable part is independent moments, already tested.

## #394 Ethan Ding — 0.172938
Model: `Joyee-64`. Vector: 0.7064 / 0.0716 / 0.1251 / -0.0108 / 0.1427 / 0.0027.
Quote: “AtlasShift Joyee-64; 相对26仅A：mean=1.5、bulk=0.4、CD4 mean/bulk=1/1；其他两context counts不变。 raw2:1:1:2/CD4/GTF; native-first; float L1<=0.001; projection<=0.03; seed20260910; 30…”
Concrete lever: Explicit context-specific mean/bulk/CD4 coefficients, raw source2:1:1:2, native-first, seed20260910, projection≤.03 and floatL1≤.001; full recipe in description.
Component: Independent moments / context coefficients.
Label-free / final transfer: NO as written: ABC splicing or manually specified A/B/C coefficients excluded.
Expected gain / cost: Unknown; hours; transferable part is independent moments, already tested.

## #400 Li Xi — 0.171981
Model: `Joyee-66`. Vector: 0.7073 / 0.0656 / 0.1245 / -0.0112 / 0.1429 / 0.0027.
Quote: “AtlasShift Joyee-66; 相对26仅A：mean=1.5、bulk=0.3、CD4 mean/bulk=1/2；其他两context counts不变。 raw2:1:1:2/CD4/GTF; native-first; float L1<=0.001; projection<=0.03; seed20260910; 30…”
Concrete lever: Explicit context-specific mean/bulk/CD4 coefficients, raw source2:1:1:2, native-first, seed20260910, projection≤.03 and floatL1≤.001; full recipe in description.
Component: Independent moments / context coefficients.
Label-free / final transfer: NO as written: ABC splicing or manually specified A/B/C coefficients excluded.
Expected gain / cost: Unknown; hours; transferable part is independent moments, already tested.

## #404 Akui Ice — 0.171069
Model: `Joyee-68`. Vector: 0.7063 / 0.0627 / 0.1247 / -0.0109 / 0.1409 / 0.0026.
Quote: “AtlasShift Joyee-68; 相对26仅C：mean=1.5、bulk=0.3、CD4 mean/bulk=1/0.5；其他两context counts不变。 raw2:1:1:2/CD4/GTF; native-first; float L1<=0.001; projection<=0.03; seed20260910; …”
Concrete lever: Explicit context-specific mean/bulk/CD4 coefficients, raw source2:1:1:2, native-first, seed20260910, projection≤.03 and floatL1≤.001; full recipe in description.
Component: Independent moments / context coefficients.
Label-free / final transfer: NO as written: ABC splicing or manually specified A/B/C coefficients excluded.
Expected gain / cost: Unknown; hours; transferable part is independent moments, already tested.

## #409 Sulgi Kim — 0.170353
Model: `ctx-lineage-rs07`. Vector: 0.6042 / 0.1392 / 0.1277 / 0.0068 / 0.1284 / 0.0159.
Quote: “Per-context lineage-weighted sources: A=CD4T x3+K562; B/C=HCT116 x3+HEK293T x2+K562+CD4T x0.5; hESC sources dropped. Real controls, exact mean, rank shift…”
Concrete lever: A CD4×3; BC HCT116×3, HEK293T×2, CD4×.5; drop hESC; rank shift .07.
Component: Context-specific source choice / DE rank.
Label-free / final transfer: NO as written: hard-coded ABC lineage recipes excluded.
Expected gain / cost: Not carried forward; hours.

## #443 Pan xy — 0.162121
Model: `h16rc06`. Vector: 0.6678 / 0.1251 / 0.1318 / -0.0612 / 0.1150 / -0.0057.
Quote: “pc5_plus_official_space_residual_correction_a0.6…”
Concrete lever: Model string identifies PC5 and official-space residual alpha .6; residual estimator unspecified.
Component: Potential pooled geometry.
Label-free / final transfer: Unknown; insufficient residual definition.
Expected gain / cost: Unknown; insufficient code.

## #458 Giovanni Marco Dall'Olio — 0.159225
Model: `XFER040 sampling calibration strength 1.25`. Vector: 0.5651 / 0.0000 / 0.1254 / 0.0191 / 0.2350 / 0.0108.
Quote: “XFER039 with the control-only sampling-calibration strength raised 1.00 to 1.25, single factor, same SPLICE006 parent, model, 239 targets,…”
Concrete lever: Control-only sampling calibration strength1→1.25,239 targets, same parent/renderer/seeds; author says local checks negative and rolls back.
Component: Sampling geometry; source reports a one-factor probe.
Label-free / final transfer: Yes in principle; exact renderer not released.
Expected gain / cost: Unidentified; hours if renderer available.

## #541 SLee — 0.144183
Model: `B1 fold-change transfer x1.75, iid with replacement`. Vector: 0.5592 / 0.0000 / 0.1099 / 0.0089 / 0.1845 / 0.0026.
Quote: “Plain B1 source-mean fold-change transfer (epsilon 0.1, scale 1.75) on each context control profile. No PC removal, no…”
Concrete lever: Source mean FC epsilon .1, amplitude1.75; IID with-replacement control seeds; no PC removal, fixed totals or sharpening; seed0.
Component: Direction/magnitude plus unpinned dispersion.
Label-free / final transfer: Yes; no challenge treated labels needed.
Expected gain / cost: Unknown; cheap but lower score than our baseline, not an isolated comparison.

## #547 Teddy Wasserman — 0.141509
Model: `summary-transfer v18 cd4t05bc`. Vector: 0.5489 / 0.0000 / 0.1044 / 0.0005 / 0.1896 / 0.0056.
Quote: “v14 recipe (a_spec 2.5) with CD4+ T refs at 1/2 per condition in contexts B/C…”
Concrete lever: A-specific2.5 parent; halve CD4 references in BC.
Component: Source weights.
Label-free / final transfer: NO as written: hand-chosen ABC recipe.
Expected gain / cost: Excluded from unseen-context rule; hours.

## #607 vDesigner — 0.132149
Model: `vCell`. Vector: 0.6080 / 0.0000 / 0.1116 / -0.0388 / 0.1127 / -0.0005.
Quote: “B3r125: K562+HCT116 same-target borrow (B3), per-target LFC rescaled to 1.25x H1 Training median knockdown norm, NB-rendered 400 cells/target,…”
Concrete lever: K562/HCT116 same-target borrow; norm1.25×H1 training median knockdown; negative-binomial400-cell renderer.
Component: Magnitude and dispersion.
Label-free / final transfer: Yes if norm uses only external H1 training.
Expected gain / cost: Unknown; hours; lower score is not a causal ablation.

## #630 Intercontinental — 0.126832
Model: `AtlasShift + co-response residual`. Vector: 0.5924 / 0.0073 / 0.1119 / -0.0149 / 0.0646 / -0.0003.
Quote: “EXP152: context-specific AtlasShift co-response residual, strength 0.25, using public RPE1 proxy perturbation deltas with official target perturbations withheld.…”
Concrete lever: AtlasShift co-response residual strength.25 from public RPE1 proxy perturbation deltas; official targets withheld.
Component: Direction / source residual transfer.
Label-free / final transfer: External-only data eligible; exact context-specific residual estimator unspecified.
Expected gain / cost: Unknown; hours for integration; .12683 is negative observational evidence.

## #676 John Carlson — 0.114400
Model: `235-tailclip-q99-s075`. Vector: 0.5086 / 0.0000 / 0.0557 / -0.0227 / 0.1457 / -0.0008.
Quote: “PROBES 235: 226 fill (same 1,729,475-entry mask, same median gain 19.5953) with |log r| winsorized at K562 own…”
Concrete lever: Same1,729,475-entry fill mask and median gain19.5953; winsorize log-ratio at K562 q99.
Component: Tail magnitude / coverage.
Label-free / final transfer: Source-only q99 can transfer; exact mask absent.
Expected gain / cost: Unknown; hours; no evidence to prioritize.

## #691 JOOHEON LEE — 0.109664
Model: `V3dmAMcFG_am0.25_ap4.0`. Vector: 0.5148 / 0.0000 / 0.0758 / -0.0203 / 0.0858 / 0.0019.
Quote: “Ridge on control co-expression embeddings (trained on Nadig Jurkat/HepG2 CRISPRi) + composite lookup transfer of same-target responses from…”
Concrete lever: Nadig-control embeddings ridge plus five-family lookup; missing Nadig genes get lookup; mean .25 plus target-specific component4; dual moments.
Component: Coverage / direction / amplitude.
Label-free / final transfer: External training rule can transfer; coefficients universal.
Expected gain / cost: Unknown; days for retraining, no positive official evidence here.

## #772 Yiquan Wang — 0.088137
Model: `wyqmath`. Vector: 0.4197 / 0.0701 / 0.0983 / -0.1687 / 0.1042 / 0.0052.
Quote: “validation iter…”
Concrete lever: No concrete reproducible rule beyond model/status name.
Component: Unidentified from description.
Label-free / final transfer: Unknown.
Expected gain / cost: No defensible gain estimate or implementation cost.

## #793 HW — 0.082440
Model: `v2_pie_unified_strat_wilcoxon`. Vector: 0.3462 / 0.0000 / 0.0460 / -0.0147 / 0.1017 / 0.0154.
Quote: “7-dataset PIE_vcc26_unified (80.5M) ensemble + full-pool template baseline alignment + library-stratified symmetric mean-invariant integer Wilcoxon rank amplifier…”
Concrete lever: 7-dataset80.5M PIE ensemble; full-pool alignment; library-stratified symmetric mean-invariant rank amplification.
Component: DE rank engineering; does not supply accurate biological direction.
Label-free / final transfer: Can be label-free if trained only on public sources; exact rule absent.
Expected gain / cost: Unknown; large model/engineering cost; .0824 is not a success.

## #847 Marcus Edwards — 0.070203
Model: `ridge-v11bsoft-c1`. Vector: 0.3335 / 0.0000 / 0.0343 / -0.0290 / 0.0932 / -0.0109.
Quote: “Rules-clean v11b ridge + public Marson CD4+ T-cell CRISPRi (Stim48hr) training line at row weight 0.5: trained on…”
Concrete lever: Ridge on HepG2/Jurkat/RPE1/CD4Stim48hr weight.5; soft inner LOLO top3 scale grid; no global shrink.
Component: External training/scale selection.
Label-free / final transfer: Yes if held-out sources genuinely excluded; our proxy does not validate their ranking.
Expected gain / cost: Unknown; days; not evidence that more training closes our gap.

## #943 dancell — 0.034641
Model: `BioFlow-R165-v12-PIE-only`. Vector: 0.1825 / 0.0000 / 0.0454 / -0.0399 / 0.0159 / 0.0039.
Quote: “PIE official A/B/C exploratory candidate; delta_p/lfc source with conservative mean and bulk amplitudes; submitted to measure cross-context generalization.…”
Concrete lever: PIE-only delta_p/LFC with unspecified conservative mean/bulk coefficients.
Component: PIE response/magnitude.
Label-free / final transfer: Potentially; exact coefficients absent.
Expected gain / cost: Unknown; poor official result, not adopted.

## #1053 VITAL-DGU — -0.003126
Model: `gl_pert_v4`. Vector: 0.0126 / 0.0000 / 0.0026 / -0.0363 / 0.0234 / -0.0210.
Quote: “Full-vector response decoder; expanded public CRISPRi data; context-balanced training; 2000 steps; seed 0; 6242 modeled genes; mean-preserving count…”
Concrete lever: Context-balanced decoder,2000steps,6242genes, mean-preserving counts.
Component: Learned direction/coverage.
Label-free / final transfer: Training provenance insufficiently specified.
Expected gain / cost: Unknown; training cost; negative official result.

## #1096 IntrInterestExtrAversion — -0.016758
Model: `VCC2026-Validation-v14`. Vector: -0.0193 / 0.0000 / -0.0136 / -0.0631 / -0.0088 / 0.0042.
Quote: “VCC2026 validation submission.…”
Concrete lever: No concrete reproducible rule beyond model/status name.
Component: Unidentified from description.
Label-free / final transfer: Unknown.
Expected gain / cost: No defensible gain estimate or implementation cost.

## #1108 kairui hu — -0.019724
Model: `S05 COV recovered111 fixed TRAIN response heads`. Vector: 0.0814 / 0.0000 / -0.0651 / -0.1082 / 0.0098 / -0.0362.
Quote: “Original COV response retained; all111 previously missing official measurement outputs restored from fixed same-scope TRAIN residual heads, with…”
Concrete lever: Recover111 missing outputs with fixed TRAIN residual heads and exact measurement projection; NTC-only counts; original6177 normalization, amplitude1.
Component: Measurement completeness versus biological recovery.
Label-free / final transfer: Controls-only destination rule stated; training scope needs verification.
Expected gain / cost: Unknown; no positive biological result.

## #1134 Daniel T. — -0.027366
Model: `VCC26 experimental 2026-10-18 pooled`. Vector: 0.2845 / 0.0000 / 0.0827 / -0.5997 / 0.1072 / -0.0388.
Quote: “Experimental baseline submission.…”
Concrete lever: No concrete reproducible rule beyond model/status name.
Component: Unidentified from description.
Label-free / final transfer: Unknown.
Expected gain / cost: No defensible gain estimate or implementation cost.

## #1323 HaveFun — -0.202879
Model: `GCF-PIE`. Vector: 0.1112 / 0.0000 / 0.0198 / -1.3093 / 0.0207 / -0.0598.
Quote: “Per-pool CellFM control sequences with a cross-attention perturbation state. Gene text, ESM2, DepMap, STRING, and NCBI enter as…”
Concrete lever: CellFM cross-attention with gene/protein/DepMap/STRING priors; call set≤256, omit shared bulk offset.
Component: Call-set truncation and lost bulk geometry.
Label-free / final transfer: Context text only when known; unseen-context rule unspecified.
Expected gain / cost: Unknown; expensive; negative result warns against this cap/zero-bulk combination.

## #1360 WANG JIZHE — -0.267921
Model: `PRiMeFlow e42 H1 g20 n8 ESM-proxy14 v2`. Vector: 0.0531 / 0.0000 / 0.0040 / -1.5853 / 0.0008 / -0.0801.
Quote: “Paired comparison with entry By7G7CJ9nh3QDdzkGQBW. Only 14 absent target conditions use fixed ESM2 cosine top1 encoded proxies. All…”
Concrete lever: Replace14 missing targets with ESM2 cosine top1 proxies;286 known outputs exact reuse; PRiMeFlow epoch42, CFG3, Euler20, H1n8,alpha.5,535-gene fill.
Component: Missing-target imputation; no isolated paired score available here.
Label-free / final transfer: External protein proxies can transfer; full source/training contract required.
Expected gain / cost: Unknown; large generative model; negative result.

## #1554 CytoAlphaThoughts — -0.997200
Model: `state40k-r2-step18000-blend-a0475`. Vector: 0.0006 / 0.0000 / -6.0000 / -0.0511 / 0.0656 / 0.0016.
Quote: “State 40k continuation, step 18000, calibrated log1p blend alpha 0.475; generated from independent trajectory selection…”
Concrete lever: STATE40k continuation step18000; independent trajectory selection; log1p blend .475.
Component: Learned generation and scale.
Label-free / final transfer: Training/selection scope unspecified.
Expected gain / cost: Unknown; expensive; negative result.
