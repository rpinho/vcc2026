"""Expected per-metric payoff of the silent tier, from the official metric mechanics plus the measured n_real
distribution (artifacts/de_count_*.csv) and the strength predictor's measured accuracy (artifacts/strength_lofo.json).

Mechanics that make a correctly identified empty-reference target valuable (docs/vcc2026-metrics-brief.md):
  jaccard : "An empty union is defined as 1" -> reference empty AND we call nothing => J_p = 1 (vs ~0 today, since
            every entry we have submitted calls ~600 genes per target).
  fidelity: reference empty and we call nothing => 0/0, the perturbation is DROPPED from the mean (neutral, not a
            gain); reference empty and we call genes => the bare fraction correct.
  nmae    : scored only where the reference gate holds >= 10 genes, so small-n_real targets are often dropped anyway.
  reach   : the pool is the reference's significant set; empty pool => the perturbation is excluded from the mean.
  pds     : a prediction whose effect is zero scores EXACTLY 0.5 for that perturbation, which is BELOW our measured
            per-target pds (~0.65-0.73 on public folds), so a zeroed row COSTS pds.
  mse     : a zeroed row contributes its own |b_p - b_ctrl|^2 as pure numerator; the baseline does the same, so the
            scaled effect is ~0 for that target (it is what the baseline does) rather than a gain.

Scaled-score arithmetic uses the published anchors (section 8 of the brief): s = (u - b) / (r - b).
"""
import json
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
ART = HERE / 'artifacts'
# anchors (brief section 8); midpoints of the quoted per-context ranges
ANCHOR = {'jaccard': (0.029, 0.399), 'fidelity': (0.5135, 0.8135), 'reach': (0.072, 0.968), 'nmae': (1.0013, 0.400),
          'pds': (0.500, 0.9555), 'mse': (0.989, 0.0365)}


def scaled(metric, u):
    b, r = ANCHOR[metric]; return (u - b) / (r - b)


def main():
    files = sorted(ART.glob('de_count_*.csv'))
    if not files:
        print('no de_count_*.csv yet'); return
    de = pd.concat([pd.read_csv(f) for f in files])
    rep = {}
    for task, d in de.groupby('task'):
        n = d.n_real.to_numpy()
        rep[task] = dict(targets=int(len(n)), cells=int(d.cells.iloc[0]), zero=float((n == 0).mean()),
                         lt10=float((n < 10).mean()), lt50=float((n < 50).mean()), median=float(np.median(n)),
                         p90=float(np.percentile(n, 90)))
        print(task, rep[task], flush=True)
    # use the 400-cell task (the official cell count) for the headline estimate
    main_task = 'H1_rich_search' if 'H1_rich_search' in rep else list(rep)[0]
    f0 = rep[main_task]['zero']; f10 = rep[main_task]['lt10']
    lofo = json.loads((ART / 'strength_lofo.json').read_text()) if (ART / 'strength_lofo.json').exists() else None
    acc = np.nanmean([r['silent_auc'] for r in lofo['lofo']]) if lofo else np.nan
    recall = np.nanmean([r['bottom_third_recall'] for r in lofo['lofo']]) if lofo else np.nan
    print(f'\nempty-reference fraction (n_real == 0) at 400 cells: {f0:.3f}; n_real < 10: {f10:.3f}')
    print(f'predictor: silent AUC {acc:.3f}, bottom-third recall {recall:.3f}')
    rows = []
    for prec in [1.0, 0.75, 0.5, recall if np.isfinite(recall) else 0.5]:
        # we silence a set of size f0/prec targets to capture f0*prec... model it directly:
        # correctly silenced share of the panel = f0 * recall_of_empty_targets (assume = prec)
        good = f0 * prec
        # wrongly silenced share: targets that DO respond but we zero. Assume we silence the predicted bottom third,
        # so wrongly silenced = (1/3 - good) at most; take the pessimistic full remainder of the silenced set.
        bad = max(0.0, 1 / 3 - good)
        dj = good * (scaled('jaccard', 1.0) - scaled('jaccard', 0.0))
        # pds: every silenced row (right or wrong) drops to exactly 0.5 => scaled 0 for that row
        per_target_pds_now = 0.70
        dpds = -(good + bad) * (scaled('pds', per_target_pds_now) - scaled('pds', 0.5))
        # fidelity: correctly silenced rows are dropped (neutral); wrongly silenced rows score 0 (reference found
        # something, we called nothing => 0)
        per_target_fid_now = 0.53
        dfid = -bad * (scaled('fidelity', per_target_fid_now) - scaled('fidelity', 0.0))
        rows.append(dict(empty_recall=round(prec, 3), silenced_correct=round(good, 4), silenced_wrong=round(bad, 4),
                         d_jaccard=round(dj, 4), d_pds=round(dpds, 4), d_fidelity=round(dfid, 4),
                         d_overall=round((dj + dpds + dfid) / 6, 5)))
    out = pd.DataFrame(rows)
    print('\nSilenced-tier payoff (scaled component deltas, and overall = sum/6):')
    print(out.to_string(index=False))
    print('\nSame, but zeroing ONLY the DE calls and keeping the bulk response (no pds cost):')
    rows2 = []
    for prec in [1.0, 0.75, 0.5, recall if np.isfinite(recall) else 0.5]:
        good = f0 * prec; bad = max(0.0, 1 / 3 - good)
        dj = good * (scaled('jaccard', 1.0) - scaled('jaccard', 0.0))
        dfid = -bad * (scaled('fidelity', 0.53) - scaled('fidelity', 0.0))
        rows2.append(dict(empty_recall=round(prec, 3), d_jaccard=round(dj, 4), d_pds=0.0, d_fidelity=round(dfid, 4),
                          d_overall=round((dj + dfid) / 6, 5)))
    print(pd.DataFrame(rows2).to_string(index=False))
    (ART / 'metric_payoff.json').write_text(json.dumps(dict(de=rep, silent_auc=acc, bottom_third_recall=recall,
                                                            full_silence=rows, calls_only=rows2, anchors=ANCHOR), indent=2, default=float))


if __name__ == '__main__':
    main()
