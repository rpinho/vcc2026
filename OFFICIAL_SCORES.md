# Official validation scores

Recorded on vcc2026-val-1 through the end of validation. Values are rounded
scaled metrics; missing components were not recorded.

| entry | overall | pds | mse | nmae | fid | reach | jac | what it was |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| v1 transfer | 0.0877 | 0.517 | 0 | -0.006 | -0.067 | 0.099 | -0.017 | public same-target signatures on the control mean (model/) |
| v2b posterior transport | 0.0617 | 0.446 | 0 | 0.042 | -0.188 | 0.092 | -0.020 | 1x signatures, transport of real control cells (model_v2/) |
| v3b mean-shift + 0.3x signatures | **-0.2482** | 0.050 | 0 | 0.014 | **-1.539** | 0.028 | -0.042 | shared average-perturbation shift on every target + shrunk signatures (model_v3/) |
| v4a AtlasShift-anchor transfer | **0.1438** | 0.678 | 0 | 0.129 | -0.019 | 0.076 | -0.001 | AtlasShift anchors + transfer (model_v4/)|
| v4b dual-summary + ash shrinkage + guide | 0.1243 | 0.658 | 0 | 0.070 | -0.025 | 0.042 | -0.000 | dual-summary anchors + ash shrinkage + guide (model_v4/)|
| v5a spread-calibrated v4a | **0.1499** | 0.678 | 0.078 | 0.094 | -0.020 | 0.071 | -0.001 | v4a anchors, per-context spread scaled to 1.0x the MSE noise credit (model_v5/)|
| v5b DE-gated + spread-calibrated + control axes removed | **0.1528** | 0.680 | 0.082 | 0.109 | -0.031 | 0.077 | -0.000 | v5a + DE calls gated by cross-source sign agreement + top-50 control axes removed (model_v5/)|
| v6 | — | — | — | — | — | — | — | arbitrary-panel builder and coverage fill; no official score recorded |
| v7a / v7b | — | — | — | — | — | — | — | native bulk / selective-scatter experiments; validation attempts refused after the phase changed |
| v8a layout-only | 0.1470 | 0.680 | 0.016 | 0.102 | -0.064 | 0.162 | -0.015 | layout changes on the exact v5b totals |
| v8b amplitude | 0.1412 | 0.708 | 0.095 | -0.059 | -0.060 | 0.178 | -0.015 | increased response amplitude |
| v9a | 0.16135 | — | — | — | — | — | — | v5b totals with confidence-ordered DE layout |
| v9b | 0.12049 | — | — | — | — | — | — | alternative v9 layout |
| slim_v5b | 0.08833 | 0.667 | 0 | -0.276 | -0.048 | 0.200 | -0.013 | real-control sampling with gain and excessive spread |
| **billmA** | **0.16438** | **0.706** | **0.061** | **0.116** | **-0.066** | **0.186** | **-0.016** | response MLP blend + v5b amplitude + v9 layout |

Best recorded validation is billmA 0.16438; historical rank 463, not a current
rank claim. Final-test scores are hidden. See README.md for the final entry
history and RESEARCH_UPDATE_20261026.md for separately labelled local diagnostics.
