# Feedback proxy V2

Local research harness for the official six metrics, with complete-destination-family
exclusion and five-split/depth-matched controls. The scientific acceptance criteria
are in ACCEPTANCE_V2.md. Local proxy scores do not establish official improvements.

Run from the repository root with the Python environment described in the top-level
README. First regenerate the synthetic scorer parity receipts in rsi/ (README
commands), then prepare public tasks and source banks with rsi/prepare_tasks.py and
proxy/tasks.py. These historical acquisition scripts need their public-source input
files and metadata; no task matrices or historical receipts are distributed.

`proxy/run_v2.py --help` lists the scoring options. The queue, calibration and
independent verification entry points are queue_v2.py, calibrate_v2.py and verify_v2.py.
They require newly generated tasks/anchors/registration rather than the private
working directory's cached state. Missing receipts deliberately fail closed.

The sorted-reference MWU optimization is excluded because its upstream license was
unspecified. execution_v2.py retains the official installed backend and exact-result
caching. Re-run check_full_adapter.py before trusting a newly generated calibration.

The code reads public data and writes local results. It does not submit to VCC.
Keep every generated matrix, raw leaderboard response, cache and log out of Git.
