"""Verify production uses original official functions, plus structural invariants."""
import runtime as rt
import official_scoring as s
import cell_eval2,cell_eval2.run as run,cell_eval2.anchor as anchor
assert s.compute_metrics is cell_eval2.compute_metrics
assert run._compute_de_side is s.ORIGINAL_DE
assert anchor.compute_metrics is cell_eval2.compute_metrics
assert not hasattr(s,'install_bounded')
rt.save('artifacts/scorer_identity.json',dict(status='PASS unmodified official implementation identity',production=s.provenance(),rejected_adapter='rejected_de_partition.json'));print('PASS unmodified scorer identity')
