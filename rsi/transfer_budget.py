"""Durable pre-I/O charges: failed requests and process restarts spend budget too."""
import fcntl
import json
import os
import threading
from pathlib import Path


class TransferBudget:
    def __init__(self, path, limit):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.file = self.path.open('a+', encoding='utf-8')
        fcntl.flock(self.file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        self.lock = threading.Lock()
        self.limit = limit
        self.total = 0
        self.bindings = {}
        self.file.seek(0)
        for line in self.file:
            row = json.loads(line)  # Torn/corrupt journals fail closed.
            if row['event'] == 'charge':
                assert row['bytes'] > 0
                self.total += row['bytes']
                assert row['total'] == self.total <= self.limit
            elif row['event'] == 'bind':
                assert row['key'] not in self.bindings
                self.bindings[row['key']] = row['value']
            else:
                raise ValueError('Unknown transfer journal event')
        self.bind('limit_bytes', limit)

    def _append(self, row):
        self.file.write(json.dumps(row, sort_keys=True) + '\n')
        self.file.flush()
        os.fsync(self.file.fileno())

    def bind(self, key, value):
        with self.lock:
            if key in self.bindings:
                assert self.bindings[key] == value, f'Transfer identity changed: {key}'
            else:
                self._append(dict(event='bind', key=key, value=value))
                self.bindings[key] = value

    def charge(self, amount, reason):
        assert isinstance(amount, int) and amount > 0
        with self.lock:
            total = self.total + amount
            assert total <= self.limit, 'Cumulative transfer budget exhausted'
            self._append(dict(event='charge', bytes=amount, total=total, reason=reason))
            self.total = total

    def close(self):
        self.file.close()
