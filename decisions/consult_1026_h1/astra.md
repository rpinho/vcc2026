VERDICT: Finish the manifold gate; reject unconditional control-mean replacement. Freeze submissions.
EXPERIMENT: If manifold fails, run one cross-fitted sign-correction classifier with fixed DE call sets.
PROTOCOL: Hold out entire targets and their replicates; train on permitted training features and labels; evaluate all six metrics.
KILL TEST: Kill sign correction unless held-out direction fidelity and composite both beat the incumbent without changing call counts.
CEILING: 0.20 as a seven-day planning estimate, not a measured upper bound.

The control replacement has **no demonstrated official gain**. H1’s raw MSE improvement would be \(1.4097-1=0.4097\), but official MSE would worsen from \(0.973\) to \(1\). These are different evaluation populations.

Let \(s\) be the official MSE scaling function. With other components unchanged, replacement gains \([s(1)-0.061]/6\). **If control scores zero**, that is **−0.01017 composite**, taking **0.16438 → 0.15421**. Expression replacement could also alter PDS and DE scores. Neither scaling nor head independence is supplied. Consequently, a numerical *expected* advantage for manifold is unidentifiable; claiming one would fabricate evidence. Its existing gate tests an actual improvement over control.

The fallback isolates sign errors: learn whether to retain or reverse each predicted sign using training-only features, retaining ranked genes and call counts. Precision **0.752/0.705** supports preserving selection; sign accuracy **0.636** leaves room, without proving learnability.

Replicates offer a conditional construction: recipient control mean plus the averaged treated-minus-matched-control response from **permitted training replicates of that same target and line**. Test on an independently acquired biological replicate, excluding shared cells, batches and label-derived preprocessing. If performance disappears under that separation, reject the claim. Same-line controls alone do not establish **0.8 response cosine**.

Reaching **0.20** requires **0.03562 composite**, or **0.21372 summed component improvement**. That is already ambitious given the failed mechanisms; leader-level recovery lacks supporting evidence.