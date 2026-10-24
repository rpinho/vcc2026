"""Single-worker final public rebuild and enrolled-metric/rational verification."""
import runtime,gc
import public_build
R=runtime.ROOT
runtime.save(R/'artifacts/public_final_manifest.json',dict(source_sha256={p:runtime.sha(R/p) for p in ['layout.py','call_budget.py','public_build.py','public_projection.py','public_compare.py','public_verify.py','scorer_check.py']},scorer_sha256=runtime.scorer_hashes()))
for source in ['Jurkat','HepG2','RPE1']:
 public_build.main(source);gc.collect()
import public_projection
import public_compare
import public_verify
runtime.save(R/'artifacts/public_final_complete.json',dict(status='PASS final layout rebuilt on three held public panels; official producer/dispatch and independent raw-metric check complete'))
