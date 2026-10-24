"""Build Bill response features for any moment-bank directory, source means/masks unchanged."""
import common as job
import io_guard;io_guard.install()
import argparse,shutil
from pathlib import Path
from finalize_sources import convert
from response_union import append,raw_priority

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--banks',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
    dest=Path(a.out).resolve();assert dest.is_relative_to(job.ROOT)
    dest.mkdir(parents=True,exist_ok=True)
    # Preserve the complete expanded source backbone (including full H1/GWPS).
    # Future panel banks add only novel target responses, with matched raw GWPS taking precedence.
    for p in sorted((job.ROOT/'data').glob('*.npz')):
        if not (dest/p.name).exists():shutil.copyfile(p,dest/p.name)
    fresh=dest/'fresh';fresh.mkdir(exist_ok=True)
    for p in sorted(Path(a.banks).glob('*.npz')):
        if p.name=='CD4_DE_statistics.npz' or p.name.startswith('K562_GWPS_bulk'):continue
        convert(p,fresh);append(dest,p.stem,fresh/p.name)
        (fresh/p.name).unlink()  # own canonical conversion intermediate; source moments and provenance remain
    raw_priority(dest)
    job.save(dest/'manifest.json',{p.stem:dict(sha256=job.sha(p),bytes=p.stat().st_size) for p in dest.glob('*.npz')})
if __name__=='__main__':main()
