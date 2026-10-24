import common as job
import io_guard;io_guard.install()
from fractions import Fraction as F

def main():
    scores=job.json.loads((job.ROOT/'artifacts/local_scores.json').read_text())
    prior=job.json.loads((job.REPO/'model_billm/artifacts/local_scores.json').read_text())
    exact=job.json.loads((job.ROOT/'local_verification.json').read_text())
    assert exact['status']=='PASS both SEARCH receipts'
    assert len(exact['records'])==2
    rows=[];passed=True
    for task in ['H1_rich_search','CD4_depth192_search']:
        base=job.json.loads((job.REPO/f'rsi/artifacts/baseline_scores/v5b_available_loco_r33/{task}/s20260910.json').read_text())
        current=scores['tasks'][task];old=prior['tasks'][task]
        assert base['reference_sha256']==current['reference_sha256']==old['reference_sha256']
        assert set(base['scaled'])==set(current['scaled'])==set(old['scaled']) and len(base['scaled'])==6
        mean=lambda r:sum(map(lambda x:F(float(x)),r['scaled'].values()),F())/6
        ok=mean(current)>mean(base);passed &= ok
        rows.append(dict(task=task,pass_vs_v5b=ok,v5b=base['scaled'],billmA=old['scaled'],billmB=current['scaled'],
            overall_v5b=float(mean(base)),overall_billmA=float(mean(old)),overall_billmB=float(mean(current)),
            exact_delta_vs_v5b=str(mean(current)-mean(base))))
    job.save(job.ROOT/'artifacts/acceptance.json',dict(status='PASS' if passed else 'FAIL',rule='Both contexts strictly exceed frozen v5b; exact stored binary64 rational six-metric means',tasks=rows,
        macro=scores['equal_context_macro'],scores_sha256=job.sha(job.ROOT/'artifacts/local_scores.json')))
    print('ACCEPTANCE',passed,flush=True)
if __name__=='__main__':main()
