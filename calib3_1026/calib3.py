"""calib3_1026: component-level transfer of the local H1 SEARCH harness to the official validation board.

Analysis only. Reads retained score files; runs no scoring, rendering or emission.
Run: nice -n 10 python -B calib3_1026/calib3.py
Writes calib3_1026/contrasts.json (and prints the tables used in REPORT.md).
"""
import csv, json, os
import numpy as np

from pathlib import Path
ROOT = str(Path(__file__).resolve().parents[1])
OUT = os.path.join(ROOT, "calib3_1026")
P = lambda *a: os.path.join(ROOT, *a)
C = ["pds_cosine", "expr_mse_unbiased_capped_norm", "de_wilcoxon_lfc_nmae",
     "de_wilcoxon_direction_fidelity_yield_raw", "de_wilcoxon_direction_reach_raw", "de_wilcoxon_sig_jaccard"]
SHORT = dict(zip(C, ["pds", "mse", "nmae", "fid", "reach", "jac"]))
DE = ["fid", "reach", "jac"]

# ---------------------------------------------------------------- official (validation vcc2026-val-1, anchors r4)
LB_KEYS = dict(pds="score_pds", mse="score_mse", nmae="score_nmae", fid="score_fid", reach="score_reach", jac="score_jac")
RAWK = dict(pds="pds_cosine", mse="expr_mse_unbiased_capped_norm", nmae="de_wilcoxon_lfc_nmae",
            fid="de_wilcoxon_direction_fidelity_yield_raw", reach="de_wilcoxon_direction_reach_raw", jac="de_wilcoxon_sig_jaccard")


def lb_entry(path, name_prefix):
    d = json.load(open(P(path)))
    for e in d["entries"]:
        if e.get("team_name") == "Ricardo Pinho" and e["model_name"].startswith(name_prefix):
            return dict(scaled={k: e[v] for k, v in LB_KEYS.items()}, raw={k: e[v] for k, v in RAWK.items()},
                        overall=e["score_avg"], source=path + " (leaderboard API row, exact)", exact=True)
    raise KeyError(name_prefix)


official = {
    "v5b": lb_entry("model_sidechain/audit/leaderboard.json", "v5b"),
    "v8a": lb_entry("model_v9/artifacts/leaderboard.json", "v8a"),
    "v9a": lb_entry("model_billm_pack/leaderboard.json", "v9a"),
    "v2b": lb_entry("research/sources/leaderboard_live.json", "v2b"),
    "slim_v5b": lb_entry("mse_1025/inputs/leaderboard.json", "SLIM"),
    # rounded 3 dp rows (OFFICIAL_SCORES.md / PUSH_1021.md); no exact API row retained
    "v8b": dict(scaled=dict(pds=0.708, mse=0.095, nmae=-0.059, fid=-0.060, reach=0.178, jac=-0.015), overall=0.1412,
                source="OFFICIAL_SCORES.md (rounded 3 dp)", exact=False),
    "billmA": dict(scaled=dict(pds=0.706, mse=0.061, nmae=0.116, fid=-0.066, reach=0.186, jac=-0.016), overall=0.16438,
                   source="OFFICIAL_SCORES.md / PUSH_1021.md (rounded 3 dp)", exact=False),
    "v5a": dict(scaled=dict(pds=0.678, mse=0.078, nmae=0.094, fid=-0.020, reach=0.071, jac=-0.001), overall=0.1499,
                source="OFFICIAL_SCORES.md (rounded)", exact=False),
    "v4a": dict(scaled=dict(pds=0.678, mse=0.0, nmae=0.129, fid=-0.019, reach=0.076, jac=-0.001), overall=0.1438,
                source="OFFICIAL_SCORES.md (rounded)", exact=False),
    "v4b": dict(scaled=dict(pds=0.658, mse=0.0, nmae=0.070, fid=-0.025, reach=0.042, jac=-0.000), overall=0.1243,
                source="OFFICIAL_SCORES.md (rounded)", exact=False),
    "v1": dict(scaled=dict(pds=0.517, mse=0.0, nmae=-0.006, fid=-0.067, reach=0.099, jac=-0.017), overall=0.0877,
               source="OFFICIAL_SCORES.md (rounded)", exact=False),
    "v3b": dict(scaled=dict(pds=0.050, mse=0.0, nmae=0.014, fid=-1.539, reach=0.028, jac=-0.042), overall=-0.2482,
                source="OFFICIAL_SCORES.md (rounded)", exact=False),
    "v9b": dict(scaled=None, overall=0.12049, source="OFFICIAL_SCORES.md (overall only)", exact=False),
}

# ---------------------------------------------------------------- local (H1_rich_search, 79 targets, local anchors)
anch = json.load(open(P("h2h2_1026/results/official_B0.json")))["anchors"]


def scale_local(raw):
    out = {}
    for c in C:
        b, r = anch["baseline"][c], anch["replicate"][c]
        s = (raw[c] - b) / (r - b)
        if c == "expr_mse_unbiased_capped_norm":
            s = max(s, 0.0)
        out[SHORT[c]] = s
    return out


def loc(path, label, note):
    d = json.load(open(P(path)))
    raw = d["raw"]
    sc = scale_local(raw)
    stored = d.get("scaled")
    if stored:  # verify our re-scaling reproduces the stored scaled values
        for c in C:
            assert abs(stored[c] - sc[SHORT[c]]) < 1e-9, (path, c, stored[c], sc[SHORT[c]])
    return dict(raw={SHORT[c]: raw[c] for c in C}, raw_mse_unbiased=raw["expr_mse_unbiased"], scaled=sc,
                overall=float(np.mean([sc[k] for k in SHORT.values()])), source=path, note=note, label=label)


FP = "final_push/artifacts/de_check/H1_rich_search/"
local = {
    "L_v5b": loc(FP + "v5b/score.json", "v5b", "v5b biology, native dense layout (~3.5-4.5k calls); same recipe as official v5b ported to 79-target LOCO (rsi v5b_available_loco_r33 digest)"),
    "L_v9a_885": loc(FP + "v9a_885/score.json", "v9a surrogate (885)", "v5b biology + v9 layout, fixed 885 budget; official v9a used v8A per-target budgets (actual val medians 557/619/652)"),
    "L_v9a_694": loc(FP + "v9a_694/score.json", "v9a/v8a surrogate (694)", "v5b biology + v9 layout, fixed 694; inventory maps it to v8a, but it uses the v9 (MSE-safe) layout, not v8's"),
    "L_v8b_like": loc(FP + "v8b_like/score.json", "v8b surrogate", "v5b biology x1.5 amplitude + v9 layout 885 (official v8b = B=claim amplitude on v8 layout)"),
    "L_billm_native": loc(FP + "incumbent/score.json", "billmA LOCO, native layout", "billmA LOCO generic incumbent, no v9 layout (= layout_cd4 H1 layout_native)"),
    "L_billm_v9_885": loc(FP + "levers_knob_none/score.json", "billmA LOCO, v9 885", "billmA LOCO generic incumbent + v9 layout 885 (= layout_cd4 H1 incumbent)"),
    "L_mainA2_450": loc(FP + "levers_budget_450/score.json", "main_a2 biology, v9 450", "budget lever"),
    "L_mainA2_650": loc(FP + "levers_budget_650/score.json", "main_a2 biology, v9 650", "budget lever"),
    "L_mainA2_885": loc(FP + "levers_base_885/score.json", "main_a2 biology, v9 885", "budget lever"),
}
for a in ["B0", "D1770", "U885", "A075", "A125", "A150", "B0s1", "B0s2", "N1", "P1aB0", "B3"]:
    local["H_" + a] = loc(f"h2h2_1026/results/official_{a}.json", a,
                          "h2h2 arm (pinned shipped billmA response, shipped v9 emitter, same 79 targets/anchors)")
for v in ["v1", "v2b", "v3b", "v4a", "v4b"]:
    d = json.load(open(P(f"rsi/artifacts/scores/H1_rich_search/{v}_s42.json")))
    local["R_" + v] = dict(scaled={SHORT[c]: d["scaled"][c] for c in C}, overall=d["overall"],
                           source=f"rsi/artifacts/scores/H1_rich_search/{v}_s42.json", note="historical RSI port", label=v)
d = json.load(open(P("rsi/artifacts/baseline_scores/v5a_available_loco_r30/H1_rich_search/s20260910.json")))
local["R_v5a"] = dict(scaled={SHORT[c]: d["scaled"][c] for c in C}, overall=d["overall"],
                      source="rsi/artifacts/baseline_scores/v5a_available_loco_r30/H1_rich_search/s20260910.json",
                      note="historical RSI port", label="v5a")

# panel identity check: h2h2 and final_push H1_rich_search
t_fp = sorted({r["perturbation"] for r in csv.DictReader(open(P(FP + "v5b/tidy.csv")))})
t_h2 = sorted(json.load(open(P("h2h2_1026/results/official_B0.json")))["targets"])
fp_v5b_raw = json.load(open(P(FP + "v5b/score.json")))["raw"]
panel_check = dict(same_79_targets=t_fp == t_h2, n=len(t_h2),
                   same_expr_distance=fp_v5b_raw["expr_distance_unbiased"] == anch["baseline"]["expr_distance_unbiased"],
                   stored_scaled_reproduced_with_h2h2_anchors=True)

# ---------------------------------------------------------------- pairs inventory
pairs = [
    ("v5b", "L_v5b", "same recipe ported to 79 LOCO targets; native layout; not the 300-target package"),
    ("v9a", "L_v9a_885", "same v5b biology + v9 layout code; fixed 885 vs official per-target v8A budgets"),
    ("v8a", "L_v9a_694", "v5b biology, 694 budget, but v9 layout code (official v8a used v8 layout, which lost 80% of MSE)"),
    ("v8b", "L_v8b_like", "x1.5 amplitude on v9 layout 885; official v8b amplitude B=claim on v8 layout"),
    ("billmA", "L_billm_v9_885", "billmA LOCO generic incumbent + v9 885 (biology construction differs from shipped)"),
    ("billmA", "H_B0", "pinned shipped billmA response + shipped v9 emitter; closest to the frozen entry"),
    ("v5a", "R_v5a", "RSI port"), ("v4a", "R_v4a", "RSI port"), ("v4b", "R_v4b", "RSI port"),
    ("v2b", "R_v2b", "RSI port"), ("v1", "R_v1", "RSI port"), ("v3b", "R_v3b", "RSI port"),
]
inventory = [dict(official=o, local=l, official_scaled=official[o]["scaled"], official_overall=official[o]["overall"],
                  official_exact=official[o]["exact"], local_scaled=local[l]["scaled"], local_overall=local[l]["overall"],
                  local_raw=local[l].get("raw"), exactness=n) for o, l, n in pairs]


# ---------------------------------------------------------------- contrasts
def dlt(a, b, src):
    return {k: a["scaled"][k] - b["scaled"][k] for k in SHORT.values()}


contrasts = {
    "C1_v5b-v9a @ 885": ("v5b", "v9a", "L_v5b", "L_v9a_885", "layout: native dense vs v9 (v5b biology)"),
    "C1b_v5b-v9a @ 694": ("v5b", "v9a", "L_v5b", "L_v9a_694", "same official contrast, 694 surrogate"),
    "C3_v5b-v8a": ("v5b", "v8a", "L_v5b", "L_v9a_694", "layout: native vs v8 (local uses v9 layout: MSE confound)"),
    "C4_v8b-v8a": ("v8b", "v8a", "L_v8b_like", "L_v9a_694", "amplitude x~1.5 at fixed layout"),
    "C5_billmA-v9a": ("billmA", "v9a", "L_billm_v9_885", "L_v9a_885", "biology: billmA vs v5b on v9 layout"),
    "C6_v8a-v9a": ("v8a", "v9a", None, None, "layout implementation v8 vs v9 at similar budgets; NO local analog"),
}
ctab = {}
for name, (oa, ob, la, lb, desc) in contrasts.items():
    od = {k: official[oa]["scaled"][k] - official[ob]["scaled"][k] for k in SHORT.values()}
    oo = official[oa]["overall"] - official[ob]["overall"]
    if la is None:
        ctab[name] = dict(desc=desc, official=od, official_overall=oo, local=None)
        continue
    ld = {k: local[la]["scaled"][k] - local[lb]["scaled"][k] for k in SHORT.values()}
    lo = local[la]["overall"] - local[lb]["overall"]
    ratio = {k: (ld[k] / od[k] if abs(od[k]) > 1e-12 else None) for k in SHORT.values()}
    sign = {k: (None if abs(od[k]) < 0.002 else bool(np.sign(ld[k]) == np.sign(od[k]))) for k in SHORT.values()}
    ctab[name] = dict(desc=desc, official=od, official_overall=oo, local=ld, local_overall=lo, ratio=ratio,
                      sign_agree=sign, overall_sign_agree=bool(np.sign(lo) == np.sign(oo)),
                      official_pair=[oa, ob], local_pair=[la, lb])

# ---------------------------------------------------------------- transfer fit + leave-one-contrast-out
INDEP = ["C1_v5b-v9a @ 885", "C3_v5b-v8a", "C4_v8b-v8a", "C5_billmA-v9a"]


def fit_rho(names, comps=("pds", "nmae", "fid", "reach", "jac")):
    rho = {}
    for k in comps:
        o = np.array([ctab[n]["official"][k] for n in names])
        l = np.array([ctab[n]["local"][k] for n in names])
        den = float(o @ o)
        rho[k] = float(l @ o / den) if den > 1e-8 else None  # local = rho * official, through origin
    return rho


def map_official(ld, rho, mse_term=0.0):
    comp = {}
    for k in ["pds", "nmae", "fid", "reach", "jac"]:
        comp[k] = 0.0 if (rho.get(k) in (None, 0.0)) else ld[k] / rho[k]
    comp["mse"] = mse_term
    return comp, float(np.mean(list(comp.values())))


loo = {}
for h in INDEP:
    train = [n for n in INDEP if n != h]
    rho = fit_rho(train)
    comp, pred = map_official(ctab[h]["local"], rho)
    off = ctab[h]["official"]
    off_nomse = float(np.mean([off[k] for k in ["pds", "nmae", "fid", "reach", "jac"]] + [0.0]))
    loo[h] = dict(train=train, rho=rho, mapped=comp, mapped_overall_ex_mse=pred,
                  official_overall=ctab[h]["official_overall"], official_overall_ex_mse=off_nomse,
                  sign_ok_full=bool(np.sign(pred) == np.sign(ctab[h]["official_overall"])),
                  sign_ok_ex_mse=bool(np.sign(pred) == np.sign(off_nomse)),
                  raw_local_sign_ok=bool(np.sign(ctab[h]["local_overall"]) == np.sign(ctab[h]["official_overall"])))

# DE-only LOO restricted to the layout contrasts (C1, C3) is degenerate (shared v5b; near-duplicate local legs).

# ---------------------------------------------------------------- official MSE slope (per raw unit) from exact rows
o5, o8, o9 = official["v5b"], official["v8a"], official["v9a"]
slope = (o5["scaled"]["mse"] - o8["scaled"]["mse"]) / (o8["raw"]["mse"] - o5["raw"]["mse"])  # scaled per unit raw decrease
billm_raw_est = o5["raw"]["mse"] + (o5["scaled"]["mse"] - official["billmA"]["scaled"]["mse"]) / slope
v9_vs_v5b_rel = o9["raw"]["mse"] / o5["raw"]["mse"] - 1


def mse_term(rel_change_local, frac):
    # official composite points if a fraction `frac` of the local relative raw-MSE change transfers to billmA's official raw
    d_raw = rel_change_local * billm_raw_est * frac
    return -d_raw * slope / 6.0


# ---------------------------------------------------------------- projections
RATIO_SETS = {
    "C1 @ 885": fit_rho(["C1_v5b-v9a @ 885"], comps=("fid", "reach", "jac", "nmae")),
    "C1 @ 694": fit_rho(["C1b_v5b-v9a @ 694"], comps=("fid", "reach", "jac", "nmae")),
    "C3": fit_rho(["C3_v5b-v8a"], comps=("fid", "reach", "jac", "nmae")),
    "pooled_C1C3C4C5": fit_rho(INDEP, comps=("fid", "reach", "jac", "nmae")),
}
# pooled nmae rho is dominated by C4 (sign-wrong); for layout arms nmae local deltas are ~1e-5 so it does not matter.
proj = {}
for arm in ["D1770", "U885", "A125", "A150", "A075", "B0s1", "B0s2"]:
    ld = {k: local["H_" + arm]["scaled"][k] - local["H_B0"]["scaled"][k] for k in SHORT.values()}
    rel = local["H_" + arm]["raw"]["mse"] / local["H_B0"]["raw"]["mse"] - 1
    res = {}
    for rs, rho in RATIO_SETS.items():
        comp = {k: (ld[k] / rho[k] if rho.get(k) not in (None, 0.0) else 0.0) for k in DE}
        comp["nmae"] = ld["nmae"] / rho["nmae"] if rho.get("nmae") and rho["nmae"] > 0 else ld["nmae"]
        comp["pds"] = ld["pds"]  # PDS local=official assumed (identical for layout arms; amplitude arms: see report)
        de_pts = sum(comp[k] for k in DE) / 6.0
        res[rs] = dict(components=comp, de_points=de_pts, nmae_pds_points=(comp["nmae"] + comp["pds"]) / 6.0)
    mse = dict(rel_raw_change_local=rel, full=mse_term(rel, 1.0), half=mse_term(rel, 0.5), zero=0.0)
    totals = {rs: res[rs]["de_points"] + res[rs]["nmae_pds_points"] + mse["half"] for rs in res}
    proj[arm] = dict(local_delta_scaled=ld, local_delta_composite=float(np.mean(list(ld.values()))),
                     by_ratio_set=res, mse_points=mse, total_central_by_ratio_set=totals,
                     range=[min(totals.values()) - abs(mse["half"]) - 0.0004 * (arm == "D1770"),
                            max(totals.values()) + abs(mse["half"])])

# ---------------------------------------------------------------- per-target pattern
cols = list(np.load(P("h2h2_1026/results/target_metrics_B0.npz"))["columns"])
def tm(a):
    z = np.load(P(f"h2h2_1026/results/target_metrics_{a}.npz"))
    return dict(zip(z["targets"], z["values"]))
B0t, Dt = tm("B0"), tm("D1770")
lc = {r["target"]: r for r in csv.DictReader(open(P("layout_cd4_1025/artifacts/per_target.csv"))) if r["recipient"] == "H1"}
def tidy(path):
    out = {}
    for r in csv.DictReader(open(P(path))):
        out.setdefault(r["perturbation"], {})[r["metric"]] = float(r["value"])
    return out
v5t, v9t = tidy(FP + "v5b/tidy.csv"), tidy(FP + "v9a_885/tidy.csv")
fi, ri, ji = (cols.index("de_wilcoxon_direction_fidelity_yield_raw"), cols.index("de_wilcoxon_direction_reach_raw"),
              cols.index("de_wilcoxon_sig_jaccard"))
T = t_h2
nconf = np.array([float(lc[t]["n_conf"]) for t in T])
strata = np.array([lc[t]["stratum"] for t in T])
def comp_arr(get_a, get_b, key):
    return np.array([get_a(t, key) - get_b(t, key) for t in T])
dD = {k: np.array([Dt[t][i] - B0t[t][i] for t in T]) for k, i in [("fid", fi), ("reach", ri), ("jac", ji)]}
dN = {k: np.array([float(lc[t][f"layout_native_{k}"]) - float(lc[t][f"incumbent_{k}"]) for t in T]) for k in ["fid", "reach", "jac"]}
dV = {k: np.array([v5t[t][m] - v9t[t][m] for t in T]) for k, m in
      [("fid", "de_wilcoxon_direction_fidelity_yield_raw"), ("reach", "de_wilcoxon_direction_reach_raw"), ("jac", "de_wilcoxon_sig_jaccard")]}
Dcalls_note = "D1770 calls per target min/median/max 378/1655/1770; B0 885 budget"

def pattern(d):
    f = d["fid"]
    tot = f.sum()
    order = np.argsort(-f)
    share5 = float(f[order[:5]].sum() / tot) if tot != 0 else None
    st = {}
    for s in ["<10", "10-999", ">=1000"]:
        m = strata == s
        st[s] = dict(n=int(m.sum()), **{k: float(d[k][m].mean()) for k in d},
                     share_of_fid_sum=float(d["fid"][m].sum() / tot) if tot else None,
                     share_of_reach_sum=float(d["reach"][m].sum() / d["reach"].sum()) if d["reach"].sum() else None,
                     share_of_jac_sum=float(d["jac"][m].sum() / d["jac"].sum()) if d["jac"].sum() else None)
    return dict(better=int((f > 1e-12).sum()), worse=int((f < -1e-12).sum()), mean={k: float(d[k].mean()) for k in d},
                top5_share_fid=share5, top5_targets=[T[i] for i in order[:5]], strata=st)

def spearman(a, b):
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])

per_target = dict(
    D1770_minus_B0=pattern(dD), billm_native_minus_v9=pattern(dN), v5b_minus_v9a885=pattern(dV),
    spearman_fid_D_vs_native=spearman(dD["fid"], dN["fid"]), spearman_fid_D_vs_v5b=spearman(dD["fid"], dV["fid"]),
    spearman_jac_D_vs_native=spearman(dD["jac"], dN["jac"]), spearman_reach_D_vs_native=spearman(dD["reach"], dN["reach"]),
    note=Dcalls_note)

out = dict(
    question="Does the component-calibrated local harness project D1770 - B0 to an official positive?",
    panel_check=panel_check, official=official,
    local={k: {kk: vv for kk, vv in v.items()} for k, v in local.items()},
    inventory=inventory, contrasts=ctab, loo=loo, ratio_sets=RATIO_SETS,
    official_mse=dict(slope_scaled_per_raw_unit=slope, billmA_raw_est=billm_raw_est,
                      v9a_vs_v5b_relative_raw_change=v9_vs_v5b_rel,
                      local_same_contrast_relative=local["L_v9a_885"]["raw"]["mse"] / local["L_v5b"]["raw"]["mse"] - 1),
    projections=proj, per_target=per_target)
json.dump(out, open(os.path.join(OUT, "contrasts.json"), "w"), indent=1, default=float)

# ---------------------------------------------------------------- print
f = lambda x: "  n/a " if x is None else f"{x:+.4f}"
print("panel", panel_check)
print("\nCONTRASTS (scaled deltas; official | local | ratio local/official)")
for n, c in ctab.items():
    print(n, c["desc"])
    print("   official", {k: round(v, 4) for k, v in c["official"].items()}, "overall %+.4f" % c["official_overall"])
    if c["local"]:
        print("   local   ", {k: round(v, 4) for k, v in c["local"].items()}, "overall %+.4f" % c["local_overall"])
        print("   ratio   ", {k: (None if v is None else round(v, 2)) for k, v in c["ratio"].items()}, "sign", c["sign_agree"])
print("\nLOO")
for h, r in loo.items():
    print(h, "rho", {k: (None if v is None else round(v, 2)) for k, v in r["rho"].items()},
          "mapped(ex-mse) %+.4f official %+.4f official_ex_mse %+.4f raw_local %+.4f" % (
              r["mapped_overall_ex_mse"], r["official_overall"], r["official_overall_ex_mse"], ctab[h]["local_overall"]),
          "sign full", r["sign_ok_full"], "ex-mse", r["sign_ok_ex_mse"], "rawlocal", r["raw_local_sign_ok"])
print("\nratio sets", {k: {kk: (None if vv is None else round(vv, 3)) for kk, vv in v.items()} for k, v in RATIO_SETS.items()})
print("official mse slope %.4f billmA raw est %.4f v9a/v5b rel %.5f" % (slope, billm_raw_est, v9_vs_v5b_rel))
print("\nPROJECTIONS (official composite points)")
for a, p in proj.items():
    print(a, "local d %+.4f" % p["local_delta_composite"], "local scaled", {k: round(v, 4) for k, v in p["local_delta_scaled"].items()})
    for rs, r in p["by_ratio_set"].items():
        print("    ", rs.ljust(16), "DE %+.5f nmae+pds %+.5f" % (r["de_points"], r["nmae_pds_points"]),
              {k: round(v, 4) for k, v in r["components"].items()})
    print("     mse rel %+.4f -> full %+.5f half %+.5f" % (p["mse_points"]["rel_raw_change_local"], p["mse_points"]["full"], p["mse_points"]["half"]),
          "totals", {k: round(v, 5) for k, v in p["total_central_by_ratio_set"].items()}, "range", [round(x, 5) for x in p["range"]])
print("\nPER TARGET")
for k in ["D1770_minus_B0", "billm_native_minus_v9", "v5b_minus_v9a885"]:
    p = per_target[k]
    print(k, "better/worse", p["better"], p["worse"], "top5 share %.2f" % p["top5_share_fid"], p["top5_targets"],
          "mean", {kk: round(v, 4) for kk, v in p["mean"].items()})
    for s, v in p["strata"].items():
        print("     ", s, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items()})
print({k: round(v, 3) for k, v in per_target.items() if k.startswith("spearman")})

# ---------------------------------------------------------------- mechanism model: n_conf-stratum reweighting
# Assumption: per-stratum local effects transfer; the official val panel differs only in its stratum mix w.
# Fit w on the one official layout contrast with exact raw components (v5b - v9a), using local v5b - v9a_885
# per-target deltas; then predict D1770 - B0 official raw deltas and scale with official slopes.
OFF_RAW = {k: official["v5b"]["raw"][k] - official["v9a"]["raw"][k] for k in DE}
OFF_SLOPE = {k: (official["v5b"]["scaled"][k] - official["v9a"]["scaled"][k]) / OFF_RAW[k] for k in DE}
v8a_check = {k: official["v9a"]["scaled"][k] + OFF_SLOPE[k] * (official["v8a"]["raw"][k] - official["v9a"]["raw"][k])
             - official["v8a"]["scaled"][k] for k in DE}
SN = ["<10", "10-999", ">=1000"]
grid = [(a, b, round(1 - a - b, 4)) for a in np.arange(0, 1.0001, 0.005) for b in np.arange(0, 1.0001 - a, 0.005)]
grid = np.array(grid)

def fit_w(dloc, idx):
    m = {k: np.array([dloc[k][idx][strata[idx] == s].mean() for s in SN]) for k in DE}
    err = sum(((grid @ m[k] - OFF_RAW[k]) * OFF_SLOPE[k]) ** 2 for k in DE)
    j = int(np.argmin(err))
    return grid[j], float(np.sqrt(err[j] / 3)), m

def predict(w, dtar, idx):
    m = {k: np.array([dtar[k][idx][strata[idx] == s].mean() for s in SN]) for k in DE}
    raw = {k: float(w @ m[k]) for k in DE}
    sc = {k: raw[k] * OFF_SLOPE[k] for k in DE}
    return raw, sc

all_idx = np.arange(len(T))
w0, rmse0, m0 = fit_w(dV, all_idx)
fit_resid = {k: float((w0 @ m0[k] - OFF_RAW[k]) * OFF_SLOPE[k]) for k in DE}
strat = {}
for arm, dtar in [("D1770", dD)] + [(a, None) for a in ["U885", "A125", "A150", "B0s1", "B0s2"]]:
    if dtar is None:
        za, zb = tm(a := arm), B0t
        dtar = {k: np.array([za[t][i] - zb[t][i] for t in T]) for k, i in [("fid", fi), ("reach", ri), ("jac", ji)]}
    raw, sc = predict(w0, dtar, all_idx)
    rng = np.random.default_rng(20261026)
    draws = []
    for _ in range(4000):
        idx = np.concatenate([rng.choice(np.where(strata == s)[0], (strata == s).sum()) for s in SN])
        w, _, _ = fit_w(dV, idx)
        _, s2 = predict(w, dtar, idx)
        draws.append([s2[k] for k in DE] + [sum(s2.values()) / 6.0])
    draws = np.array(draws)
    strat[arm] = dict(raw=raw, scaled=sc, de_points=sum(sc.values()) / 6.0,
                      de_points_ex_reach=(sc["fid"] + sc["jac"]) / 6.0,
                      boot_de_points_q05_q50_q95=np.quantile(draws[:, 3], [0.05, 0.5, 0.95]).tolist(),
                      boot_fid_scaled_q05_q95=np.quantile(draws[:, 0], [0.05, 0.95]).tolist(),
                      boot_reach_scaled_q05_q95=np.quantile(draws[:, 1], [0.05, 0.95]).tolist(),
                      boot_jac_scaled_q05_q95=np.quantile(draws[:, 2], [0.05, 0.95]).tolist(),
                      p_de_points_le_0=float((draws[:, 3] <= 0).mean()))
# check: does the same w predict official v5b - v8a from local v5b - v9a_694 (v8-layout confound)?
v694t = tidy(FP + "v9a_694/tidy.csv")
d694 = {k: np.array([v5t[t][m] - v694t[t][m] for t in T]) for k, m in
        [("fid", "de_wilcoxon_direction_fidelity_yield_raw"), ("reach", "de_wilcoxon_direction_reach_raw"), ("jac", "de_wilcoxon_sig_jaccard")]}
_, c3pred = predict(w0, d694, all_idx)
c3_check = dict(predicted_scaled=c3pred, official_scaled={k: ctab["C3_v5b-v8a"]["official"][k] for k in DE})
strat_model = dict(official_raw_v5b_minus_v9a=OFF_RAW, official_slope_scaled_per_raw=OFF_SLOPE,
                   slope_check_on_v8a_residual=v8a_check, local_strata_counts={s: int((strata == s).sum()) for s in SN},
                   fitted_w=dict(zip(SN, w0.tolist())), fit_residual_scaled=fit_resid, arms=strat, c3_check=c3_check)
out["stratum_model"] = strat_model
json.dump(out, open(os.path.join(OUT, "contrasts.json"), "w"), indent=1, default=float)
print("\nSTRATUM MODEL")
print("official raw", {k: round(v, 5) for k, v in OFF_RAW.items()}, "slope", {k: round(v, 3) for k, v in OFF_SLOPE.items()},
      "v8a slope check", {k: round(v, 5) for k, v in v8a_check.items()})
print("w", dict(zip(SN, w0.tolist())), "resid", {k: round(v, 4) for k, v in fit_resid.items()})
for a, s in strat.items():
    print(a, "scaled", {k: round(v, 4) for k, v in s["scaled"].items()}, "DE pts %+.5f ex-reach %+.5f" % (s["de_points"], s["de_points_ex_reach"]),
          "boot q05/50/95", [round(x, 5) for x in s["boot_de_points_q05_q50_q95"]], "P<=0 %.3f" % s["p_de_points_le_0"],
          "fid", [round(x, 4) for x in s["boot_fid_scaled_q05_q95"]], "reach", [round(x, 4) for x in s["boot_reach_scaled_q05_q95"]],
          "jac", [round(x, 4) for x in s["boot_jac_scaled_q05_q95"]])
print("C3 check", {k: round(v, 4) for k, v in c3pred.items()}, c3_check["official_scaled"])
