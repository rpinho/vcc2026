"""Summarize measured raw ablations. No fit or candidate selection on hidden data.
Exact Fraction differences below certify arithmetic on recorded binary64 outputs,
not the scientific uncertainty or exactness of the statistical estimator.
"""
import pathlib,json,math,statistics,csv
from fractions import Fraction
P=pathlib.Path(__file__).resolve().parent
NAMES={'pds':'pds_cosine','mse':'expr_mse_unbiased_capped_norm','nmae':'de_wilcoxon_lfc_nmae',
       'fid':'de_wilcoxon_direction_fidelity_yield_raw','reach':'de_wilcoxon_direction_reach_raw','jac':'de_wilcoxon_sig_jaccard'}
LOWER={'mse','nmae'}
ARM_ORDER=['pool4','pool1','select100','select300','select1000','selectadaptive','selectmagnitude','compress300','amplitude09','shaped300','shapedadaptive','quant1000','quant1000_exchange']
def main():
    lines=['# Local DE ablations — raw official metrics','',
           'Higher PDS/fidelity/reach/Jaccard and lower MSE/NMAE are better. These are',
           'small-panel public-line diagnostics, not scaled validation forecasts.','']
    exact={}
    for folder in sorted((P/'artifacts').iterdir()):
        if not folder.is_dir() or not (folder/'pool4.json').exists():continue
        base=json.loads((folder/'pool4.json').read_text());design=json.loads((folder/'design.json').read_text())
        targets=design['task']['targets'];exact[folder.name]={}
        cohorts={}
        for arm in ARM_ORDER:
            cp=folder/(arm+'_per_target.csv')
            if not cp.exists():continue
            with cp.open() as f:rs=list(csv.DictReader(f))
            cohorts[arm]={k:sorted(r['perturbation'] for r in rs if r['metric']==NAMES[k]
                         and r.get('value') and math.isfinite(float(r['value']))) for k in ('nmae','reach')}
        if 'pool4' in cohorts:
            for arm,cohort in cohorts.items():
                assert cohort==cohorts['pool4'],('reference-defined cohort changed',folder.name,arm,cohort,cohorts['pool4'])
        lines.extend([f'## {folder.name}', '',f'{len(targets)} targets; real cells per target '+
                      f"{min(design['task']['real_counts'].values())}–{max(design['task']['real_counts'].values())}; "+
                      f"{design['task']['controls']} controls. Destination family excluded from sources.",'',
                      '| Arm | PDS | MSE ↓ | NMAE ↓ | Fidelity | Reach | Jaccard | Calls* | Seconds |',
                      '|---|---:|---:|---:|---:|---:|---:|---:|---:|'])
        for arm in ARM_ORDER:
            path=folder/(arm+'.json')
            if not path.exists():continue
            r=json.loads(path.read_text());raw=r['raw'];vals={k:raw[n] for k,n in NAMES.items()}
            calls=[x['predicted_calls'] for x in r['de_counts'] if x['target'] in targets]
            label='amplitude09 (bulk ratios ×1.5)' if arm=='amplitude09' else arm
            row=[label]+[f'{v:.6f}' for v in vals.values()]+[f'{statistics.mean(calls):.1f}',f"{r['seconds']:.1f}"]
            lines.append('| '+' | '.join(row)+' |')
            ds={k:Fraction.from_float(float(v))-Fraction.from_float(float(base['raw'][NAMES[k]])) for k,v in vals.items()
                if math.isfinite(v) and math.isfinite(base['raw'][NAMES[k]])}
            exact[folder.name][arm]=dict(delta_binary64_rational={k:str(v) for k,v in ds.items()},
                improves_raw_components=[k for k,v in ds.items() if (v<0 if k in LOWER else v>0)],
                worsens_raw_components=[k for k,v in ds.items() if (v>0 if k in LOWER else v<0)],
                seconds=r['seconds'])
        lines.extend(['','*Mean actual q<.05 calls before own-target exclusion; only treated-target rows counted.',
                      'Official metric columns already apply own-target exclusion. Support mask size is not actual call count.',''])
        if 'pool4' in cohorts:
            lines.extend([f"Reference-defined scored cohorts: NMAE {len(cohorts['pool4']['nmae'])}, reach {len(cohorts['pool4']['reach'])}; identical across completed arms.",''])
    (P/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    (P/'artifacts/delta_certificate.json').write_text(json.dumps(exact,indent=2))
    print('\n'.join(lines))
if __name__=='__main__':main()
