# Official validation scores

Source: official submissions on vcc2026-val-1, recorded 2026-10-18. Values are rounded scaled metrics.

| entry | overall | pds | mse | nmae | fid | reach | jac | what it was |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| v1 transfer | 0.0877 | 0.517 | 0 | -0.006 | -0.067 | 0.099 | -0.017 | public same-target signatures on the control mean (model/) |
| v2b posterior transport | 0.0617 | 0.446 | 0 | 0.042 | -0.188 | 0.092 | -0.020 | 1x signatures, transport of real control cells (model_v2/) |
| v3b mean-shift + 0.3x signatures | **-0.2482** | 0.050 | 0 | 0.014 | **-1.539** | 0.028 | -0.042 | shared average-perturbation shift on every target + shrunk signatures (model_v3/) |
| v4a AtlasShift-anchor transfer | **0.1438** | 0.678 | 0 | 0.129 | -0.019 | 0.076 | -0.001 | AtlasShift anchors + transfer (model_v4/)|
| v4b dual-summary + ash shrinkage + guide | 0.1243 | 0.658 | 0 | 0.070 | -0.025 | 0.042 | -0.000 | dual-summary anchors + ash shrinkage + guide (model_v4/)|
| v5a spread-calibrated v4a | **0.1499** | 0.678 | 0.078 | 0.094 | -0.020 | 0.071 | -0.001 | v4a anchors, per-context spread scaled to 1.0x the MSE noise credit (model_v5/)|
| v5b DE-gated + spread-calibrated + control axes removed | **0.1528** | 0.680 | 0.082 | 0.109 | -0.031 | 0.077 | -0.000 | v5a + DE calls gated by cross-source sign agreement + top-50 control axes removed (model_v5/)|

v6 has no official score recorded. Ranks are omitted because they change.
