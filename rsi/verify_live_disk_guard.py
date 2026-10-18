"""No real waits or downloads: reproduce the obsolete reserve failure and live pauses."""
from collections import namedtuple
from datetime import datetime, timezone
from unittest.mock import patch

import runtime as rt
from supervisor import gates

Usage = namedtuple('Usage', 'total used free')
GB = 10**9
checks = []

def check_live(kind, free_values, job_values, extra=0, sleeps=0):
    with patch.dict(rt.os.environ, VCC_RSI_STEP='fixture', VCC_RSI_STEP_KIND=kind), \
         patch.object(rt.shutil, 'disk_usage', side_effect=[Usage(100*GB, 0, n*GB) for n in free_values]), \
         patch.object(gates, 'job_bytes', side_effect=[n*GB for n in job_values]), \
         patch.object(rt.time, 'sleep') as sleep:
        rt.guard(extra)
        assert sleep.call_count == sleeps
        assert all(call.args == (30,) for call in sleep.call_args_list)

check_live('download', [35], [2])
checks.append('35 GB permits a live download below the retired 40 GiB reserve')
check_live('download', [29, 30], [2, 2], sleeps=1)
checks.append('downloads wait below 30 GB and resume at the floor')
check_live('compute', [19, 20], [2, 2], sleeps=1)
checks.append('compute waits below 20 GB and resumes at the floor')
check_live('download', [50, 50], [21, 20], sleeps=1)
checks.append('job cap pauses and resumes without aborting')
check_live('download', [31, 32], [2, 2], extra=2*GB, sleeps=1)
check_live('compute', [50, 50], [19, 18], extra=2*GB, sleeps=1)
checks.append('pending writes reserve both free space and job capacity')
with patch.dict(rt.os.environ, clear=True), \
     patch.object(rt.shutil, 'disk_usage', return_value=Usage(100*GB, 0, 35*GB)):
    try:
        rt.guard()
    except RuntimeError as e:
        assert str(e) == '40 GiB disk reserve'
    else:
        raise AssertionError('Non-loop guard unexpectedly changed')
checks.append('non-loop one-shot guard retains its original reserve')
rt.save('artifacts/live_disk_guard_fixture.json', dict(
    status='PASS', utc=datetime.now(timezone.utc).isoformat(), checks=checks,
    runtime_sha256=rt.sha('runtime.py'), runner_sha256=rt.sha('rsi_run.py'),
    gates_sha256=rt.sha('supervisor/gates.py'),
    limitation='Mocked disk transitions; supervisor SIGSTOP/SIGCONT integration is not exercised'))
print(f'PASS live disk guard ({len(checks)} checks)')
