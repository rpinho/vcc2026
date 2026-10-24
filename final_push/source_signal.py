"""Cache the noise-corrected per-(source, target) signal strength of every screen, so the tier recipe can read it
without re-reducing whole screens. Writes artifacts/source_signal.csv."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import banks as B
from strength import screen_strength

HERE = Path(__file__).resolve().parent

rows = []
for name, bank in B.discover(include_shards=False).items():
    s, c, e, sig, n = screen_strength(bank)
    med = float(np.median(n))
    for t, v in s.items():
        rows.append(dict(source=name, family=bank.family, target=t, signal=v['signal'], energy=v['energy'],
                         n_cells=v['n_cells'], depletion=v['n_cells'] / max(med, 1), floor_c=c))
    print(name, len(s), 'floor', round(c, 4), flush=True)
pd.DataFrame(rows).to_csv(HERE / 'artifacts/source_signal.csv', index=False)
print('wrote', len(rows))
