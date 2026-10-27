**R12 verdict: B0 stays final. Close the biology line. Don't build an arm.**

**1. Is the arithmetic right?** The conclusion holds, and the real gain is probably smaller than +0.0004.
- **Unit mismatch.** Cross is a fraction of the true delta energy (Σ‖b‖²). The 0.811 rate converts raw MSE, which is measured per cell and includes cell noise. Treating 2×0.0027 as raw MSE assumes delta energy is the whole raw MSE. It is a fraction of it, so +0.0044 scaled is an upper bound. The ÷6 to the overall score and the PDS conversion (0.27–0.6) look consistent.
- **The local sign isn't trustworthy.** No local MSE instrument passes C4: with the clamp lifted it reads −0.0012 where official was +0.0791. We have no local check that +0.0027 cross would show up officially.
- **The gain is small and partial.** control_bounded only changes the 6,006 shared genes, and its H1 PDS interval reaches −0.0052. Even under the most generous conversion it is about 4× below the +0.0028 seed floor.
- Full maps gain +0.03 cross on RPE1 and HepG2, but they fail PDS in all five contexts and lose cross on K562, H1 and Jurkat. That doesn't transfer reliably.

No finding justifies building before 11-03.

**2. Any missed biology lever?** I don't see a credible one. The direction error (58%) can only be fixed with information about how a target responds in the new context. D/E/F give us only control cells, and control_bounded is the strongest control-only test we have: +0.0027. Coverage is capped at 0.5%, and shrinkage, amplitude and SLIM are already closed. Gene-network or GO priors would be another form of the family correction, and that turned out to be shrinkage (cross fell).

**3. Close it?** Yes. B0 = billmA 0.16438 is the final entry.

Finding 4 (kNN 0.423 and ridge 0.398 with in-context truth) is still worth a sentence in a write-up: better biology here depends on data the final phase doesn't provide.
