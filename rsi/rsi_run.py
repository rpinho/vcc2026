#!/usr/bin/env python3
"""Launch ONE heavy step of the live VCC RSI loop, the only sanctioned way to run long work.

    ../../.venv/bin/python rsi_run.py --name acquire_h1_rich --kind download [--detach] -- ../../.venv/bin/python acquire_h1_rich.py

- One instance per --name (flock on locks/<name>.lock). A duplicate exits 75 at once: no duplicate-worker races.
- At most 3 heavy steps at a time (slot locks); waits for a free slot.
- Waits (never exits) while a kids' job runs (quest draft / hunts / games stage-promote-deploy),
  and, for --kind download, while free disk < 30 GB or rsi/ exceeds its ~20 GB live cap.
- Runs the command under `nice -n 19 taskpolicy -b`; log -> logs/steps/<name>.log; start/end -> ledger/steps.jsonl.
- --detach: returns immediately; the step survives round boundaries (the supervisor keeps rsi_run trees
  alive between rounds and kills any other leftover process in this directory).
- The supervisor may SIGSTOP the step while kids' jobs run and SIGCONT it afterwards. Expect wall-clock gaps.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

RSI = Path(__file__).resolve().parent
sys.path.insert(0, str(RSI / "supervisor"))
import gates  # noqa: E402

LOCKS = RSI / "locks"
STEPLOG = RSI / "logs" / "steps"
LEDGER = RSI / "ledger" / "steps.jsonl"


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def record(row: dict) -> None:
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    with LEDGER.open("a") as f:
        f.write(json.dumps(row) + "\n")


def try_lock(path: Path):
    fd = os.open(path, os.O_CREAT | os.O_RDWR, 0o644)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return fd
    except BlockingIOError:
        os.close(fd)
        return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--kind", choices=["compute", "download"], default="compute")
    ap.add_argument("--detach", action="store_true")
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd[:1] == ["--"] else a.cmd
    if not cmd:
        ap.error("no command after --")
    if not all(c.isalnum() or c in "-_." for c in a.name):
        ap.error("--name must be [A-Za-z0-9._-]")
    LOCKS.mkdir(exist_ok=True)
    STEPLOG.mkdir(parents=True, exist_ok=True)
    log = STEPLOG / f"{a.name}.log"

    if a.detach:
        # Probe the name lock first so a duplicate is reported synchronously.
        fd = try_lock(LOCKS / f"{a.name}.lock")
        if fd is None:
            print(f"ALREADY RUNNING: {a.name}", file=sys.stderr)
            return 75
        os.close(fd)
        args = [sys.executable, str(Path(__file__).resolve()), "--name", a.name, "--kind", a.kind, "--", *cmd]
        with open(log, "ab") as out:
            p = subprocess.Popen(args, cwd=os.getcwd(), stdin=subprocess.DEVNULL, stdout=out,
                                 stderr=subprocess.STDOUT, start_new_session=True)
        print(f"DETACHED {a.name} pid={p.pid} log={log}")
        return 0

    name_fd = try_lock(LOCKS / f"{a.name}.lock")
    if name_fd is None:
        print(f"ALREADY RUNNING: {a.name}", file=sys.stderr)
        return 75
    os.ftruncate(name_fd, 0)
    os.write(name_fd, f"{os.getpid()}\n".encode())
    record({"t": now(), "event": "queued", "name": a.name, "kind": a.kind, "pid": os.getpid(), "cmd": cmd})

    slot_fd = None
    waited_since = None
    while True:
        reason = gates.gate_reason(a.kind)
        if reason is None and slot_fd is None:
            for i in range(gates.MAX_HEAVY_SLOTS):
                slot_fd = try_lock(LOCKS / f"slot-{i}.lock")
                if slot_fd is not None:
                    break
            if slot_fd is None:
                reason = f"all {gates.MAX_HEAVY_SLOTS} heavy slots busy"
        if reason is None:
            break
        if waited_since is None:
            waited_since = time.time()
            print(f"[{now()}] {a.name} waiting: {reason}", flush=True)
        time.sleep(30)

    print(f"[{now()}] {a.name} start: {' '.join(cmd)}", flush=True)
    record({"t": now(), "event": "start", "name": a.name, "kind": a.kind, "pid": os.getpid(),
            "cmd": cmd, "cwd": os.getcwd(), "waited_s": round(time.time() - waited_since) if waited_since else 0})
    t0 = time.time()
    env = dict(os.environ, VCC_RSI_STEP=a.name, VCC_RSI_STEP_KIND=a.kind)
    child = subprocess.Popen(gates.WORKER_PREFIX + cmd, cwd=os.getcwd(), stdin=subprocess.DEVNULL, env=env)
    def _fwd(sig, _f):
        try:
            child.send_signal(sig)
        except ProcessLookupError:
            pass
    signal.signal(signal.SIGTERM, _fwd)
    signal.signal(signal.SIGINT, _fwd)
    rc = child.wait()
    record({"t": now(), "event": "end", "name": a.name, "rc": rc, "wall_s": round(time.time() - t0)})
    print(f"[{now()}] {a.name} end rc={rc}", flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
