"""Prepare a corrected baseline independently; upgrade only after old scoring exits."""
import runtime as rt
from pathlib import Path
import argparse
import official_scoring as s,bounded_de,bounded_score
from benchmark import setup,response_baseline,anchors,rescale_existing
p=argparse.ArgumentParser();p.add_argument('--task',required=True);p.add_argument('--upgrade',action='store_true');a=p.parse_args();m,ctrl,real=setup(a.task)
assert m['split']=='search','HOLDOUT scales are prepared only within the final audit'
del ctrl;bounded_de.install();bounded_score.install_anchor();folder=Path('artifacts/scores')/a.task;folder.mkdir(parents=True,exist_ok=True)
if a.upgrade:
 r=anchors(a.task,real,folder);rescale_existing(folder,r);print(__import__('verify_results').calibrate())
else:print(response_baseline(a.task,real,folder))
