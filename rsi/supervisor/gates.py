"""Disk and concurrency gates for local background experiments.
Machine-specific workload hooks are disabled in this export.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

RSI = Path(__file__).resolve().parent.parent
SUP = RSI / "supervisor"
GB = 1000 ** 3

FREE_DOWNLOAD_FLOOR = 30 * GB   # pause new downloads below this
FREE_HARD_FLOOR = 20 * GB       # pause every heavy step below this
JOB_CAP = 20 * GB               # live-loop retention waiver cap for rsi/
MAX_HEAVY_SLOTS = 3             # concurrent heavy steps (16 GB RAM machine)
NICE = 19                       # same tier as the kids' quest flow, which runs nice 19 + taskpolicy -b
WORKER_PREFIX = ["nice", "-n", str(NICE)] + (["taskpolicy", "-b"] if shutil.which("taskpolicy") else [])

# Anchored on a real `node` argv0 so prompt text / shell command lines that merely mention these never match.
_NODE = r"^(\S*/)?node\s+"
KID_PATTERNS = []


def ps_table() -> list[tuple[int, int, str]]:
    out = subprocess.run(["ps", "-axww", "-o", "pid=,ppid=,command="],
                         capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        parts = line.strip().split(None, 2)
        if len(parts) < 3:
            continue
        try:
            rows.append((int(parts[0]), int(parts[1]), parts[2]))
        except ValueError:
            continue
    return rows


KID_STALE_S = 4 * 3600  # a kids' job older than this (draft cap is 150 min) is presumed hung and stops gating


def _etime_s(pid: int) -> int:
    out = subprocess.run(["ps", "-o", "etime=", "-p", str(pid)], capture_output=True, text=True).stdout.strip()
    if not out:
        return 0
    days, _, hms = out.rpartition("-")
    parts = [int(x) for x in hms.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    return (int(days) if days else 0) * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


def kid_jobs(rows=None) -> list[str]:
    """Return command lines of running kids' jobs that must win the machine."""
    rows = rows if rows is not None else ps_table()
    hits = []
    for pid, _ppid, cmd in rows:
        if pid == os.getpid():
            continue
        if any(p.search(cmd) for p in KID_PATTERNS) and _etime_s(pid) < KID_STALE_S:
            hits.append(f"{pid} {cmd[:160]}")
    return hits


def free_bytes() -> int:
    return shutil.disk_usage(RSI).free


def job_bytes(max_age_s: float = 900) -> int:
    """Size of rsi/ (cached by the supervisor in supervisor/disk.json; recomputed if stale)."""
    cache = SUP / "disk.json"
    try:
        d = json.loads(cache.read_text())
        if time.time() - d["t"] < max_age_s:
            return int(d["job_bytes"])
    except Exception:
        pass
    return measure_job_bytes()


def measure_job_bytes() -> int:
    try:
        out = subprocess.run(["du", "-sk", str(RSI)], capture_output=True, text=True).stdout
        kb = int(out.split()[0])
    except Exception:
        return 0
    b = kb * 1024
    try:
        tmp = SUP / "disk.json.tmp"
        tmp.write_text(json.dumps({"t": time.time(), "job_bytes": b, "free_bytes": free_bytes()}))
        tmp.replace(SUP / "disk.json")
    except Exception:
        pass
    return b


def gate_reason(kind: str = "compute") -> str | None:
    """None if a new heavy step of this kind may start now, else why it must wait."""
    kids = kid_jobs()
    if kids:
        return "kids' job running: " + kids[0]
    free = free_bytes()
    if free < FREE_HARD_FLOOR:
        return f"free disk {free / GB:.1f} GB < {FREE_HARD_FLOOR / GB:.0f} GB hard floor"
    if kind == "download":
        if free < FREE_DOWNLOAD_FLOOR:
            return f"free disk {free / GB:.1f} GB < {FREE_DOWNLOAD_FLOOR / GB:.0f} GB download floor"
        jb = job_bytes()
        if jb > JOB_CAP:
            return f"job dir {jb / GB:.1f} GB > {JOB_CAP / GB:.0f} GB live cap; prune regenerable intermediates"
    return None
