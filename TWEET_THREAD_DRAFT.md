# Thread draft for Ricardo to edit and post

Draft only. No posts have been published by this job. Each numbered paragraph is
one post, including its number; each is at most 280 characters.

1/9 We opened our VCC 2026 research repo, including the failures. Official validation: 0.0877 → 0.16438 (billmA). Final scores are hidden. Code, results and a retrospective: https://github.com/rpinho/vcc2026

2/9 Our biggest early gain came from building on public methods: AtlasShift transfer/count generation, then spread calibration and DE layout. We should have reproduced stronger baselines earlier. Credits and license boundaries are in the repo.

3/9 The local scorer was mathematically correct and still misleading. Layout deltas inflated fidelity 9–12× and Jaccard 7–9×, while understating reach. Component signs could agree and the overall ranking could reverse.

4/9 Local MSE was often clamped to zero. Lifting the clamp exposed the wrong sign on a historical amplitude contrast. None of our tested replacement instruments passed calibration. Numerical parity did not validate candidate ranking.

5/9 Doubling DE calls gained +0.03449 locally. We did not promote it: the gain depended on strong targets, and calibrated projections were near seed noise. It remains indeterminate, not a demonstrated official win or loss.

6/9 Biology audit: mean response cosine 0.084 despite retrieval PDS ≈0.856. Already-covered targets carried 99.5% of error. Filling missing source rows alone could fix at most 0.5% on that panel; better donors are a separate question.

7/9 Large distance gains were mostly shrinkage. In-context truth raised aligned cross 0.115→0.423, but isn't available for final contexts. Leave-one-context-out transfer added only +0.0027, with PDS risk. No credible official-score conversion.

8/9 AI advisors helped catch errors, and also shared wrong premises. We corrected scale invariance, fixed-pseudobulk MSE invariance, and repeated suggestions to use a closed validation server. Written tests mattered more than agreement.

9/9 What we'd change: reproduce public baselines first, validate historical rankings, save every component, hold out contexts, and separate biology from rendering. Negative results are the point. No data or weights: https://github.com/rpinho/vcc2026
