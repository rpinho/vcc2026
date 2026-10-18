"""Store displayed integer-derived moment bounds outward as exact binary64 rationals.
The generator computes integer-rational endpoints then float(Fraction) rounds to nearest.
One nextafter step outward encloses that conversion's rounding interval. No predictions change.
"""
import runtime,csv,io,math
from fractions import Fraction as F
rt=runtime.ROOT;records=[]
for arm in ['a','b']:
 path=rt/f'artifacts/prediction_v4{arm}_moments.csv';original=path.read_text();rows=list(csv.DictReader(io.StringIO(original)));assert len(rows)==900
 original_path=path.with_name(path.stem+'_nearest.csv')
 if 'mean_L1_upper_exact' in rows[0]:
  records.append(dict(arm=arm,rows=len(rows),original_sha256=runtime.sha(original_path),outward_sha256=runtime.sha(path)));continue
 original_path.write_text(original)
 for row in rows:
  for name in ['mean_L1_lower','mean_L1_upper','bulk_L1_lower','bulk_L1_upper']:
   x=float(row[name]);assert math.isfinite(x) and x>=0
   y=max(0.,math.nextafter(x,math.inf if name.endswith('upper') else -math.inf));row[name]=repr(y);row[name+'_exact']=str(F(y))
  assert F(row['mean_L1_lower_exact'])<=F(row['mean_L1_upper_exact']) and F(row['bulk_L1_lower_exact'])<=F(row['bulk_L1_upper_exact'])
  if arm=='b':assert F(row['mean_L1_upper_exact'])<F(1,50)
 with path.open('w') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
 records.append(dict(arm=arm,rows=len(rows),original_sha256=runtime.sha(original_path),outward_sha256=runtime.sha(path)))
runtime.save(rt/'artifacts/moment_table_encoding_verification.json',dict(status='VERIFIED exact rational storage of outward-expanded moment endpoints; no biological improvement claimed',method='Integer-derived rational endpoint was rounded nearest by float(Fraction); nextafter outward encloses that conversion. Exact columns store Fraction(outward binary64). Original endpoints retained.',files=records))
print('Outward rational moment-table representation VERIFIED',flush=True)
