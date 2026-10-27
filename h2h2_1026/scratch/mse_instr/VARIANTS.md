# Prespecified diagnostic variants

No parameter is fitted on either the 79 targets or the official contrasts. All formulations use the installed 50,000 group-sum log1p space; target-gene coordinates are omitted from squared distances and budget while scalar jackknives remain whole, exactly as installed. Use weights w_p=1/G_p and D=sum w_p(||real_p-real_ctrl||²-C_p-C_ctrl). In these units E=sum w_p||pred_p-real_p||²; J=sum w_p Chat_p; R=sum w_p C_p; Q=sum w_p min(Chat_p,C_p).

1. **Current clamp**: installed capped ratio, scaled with existing SEARCH baseline/replicate anchors and clipped to [0,1].
2. **Capped, score clamp lifted**: (E-R-min(B,Q))/D, affine scaling without [0,1] clamp.
3. **Uncapped unbiased**: (E-R-J)/D; neither trace cap nor budget.
4. **Trace cap lifted, budget retained**: (E-R-min(B,J))/D.
5. **Budget credit, labels**: (E-R-rho Q)/D, rho=min(1,B/Q). Must reproduce #2, an algebra identity.
6. **Budget credit, label-free rho**: same labeled residual and real subtraction, but Q0=sum w_p min(Chat_p,C0(n_p)), rho0=min(1,B/Q0), and prediction credit rho0 Q0. C0(n_p) is the average of 16 input-control bootstrap jackknives at real n_p, seed 20261026. Only rho/claim are label-free; reading MSE still requires labels. Also record the alternative rho0 applied to labeled Q to disambiguate that interpretation.
7. **Plain pseudobulk delta MSE**: sum w_p||(pred_p-input_ctrl)-(real_p-real_ctrl)||² / D; no variance credit. Also report direct plain expression error E/D.

Contrast scale: apply the existing local anchor span to every D-normalized variant, without refitting anchors, and divide the resulting difference by the supplied official scaled difference. This is a fixed common ruler, not a calibrated official transfer coefficient. Report raw differences too.

C1=v5b-v9a885; C3=v5b-v9a694; C4=v8b_like-v9a694; C5=billm_v9_885-v9a885. For C6, v9a694-v9a885 is explicitly a **budget surrogate**, not a genuine v8-v9 layout replay. A missing genuine C6 cannot pass the instrument gate. LOO uses the same fixed formulation holding out each contrast; no fitted mapping or selection on remaining official outcomes is allowed. Its signs therefore equal the direct signs, and shared legs make these checks dependent.

Qualification: C4/C3/C6 signs correct and C4 ratio 0.3–3. No gate is relaxed for missing evidence. Only a qualifying instrument triggers A125/A150/1.38x official predictions.
