"""Disk-backed payload arrays; publish progress only after array bytes are durable."""
import json
import os
from pathlib import Path

import numpy as np


class DownloadCheckpoint:
    def __init__(self, root, specs, identity):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.path = self.root / 'state.json'
        self.identity = dict(identity=identity, arrays={
            key: dict(shape=list(shape), dtype=np.dtype(dtype).str)
            for key, (shape, dtype) in specs.items()})
        saved = json.loads(self.path.read_text()) if self.path.exists() else None
        if saved is not None:
            assert saved['identity'] == self.identity, 'Payload checkpoint identity changed'
        self.arrays = {}
        for key, (shape, dtype) in specs.items():
            assert key.isidentifier()
            path = self.root / (key + '.npy')
            if saved is None:
                array = np.lib.format.open_memmap(path, mode='w+', dtype=dtype, shape=shape)
            else:
                array = np.lib.format.open_memmap(path, mode='r+')
                assert array.shape == tuple(shape) and array.dtype == np.dtype(dtype)
            self.arrays[key] = array
        self.next_plan = saved['next_plan'] if saved is not None else 0
        assert isinstance(self.next_plan, int) and self.next_plan >= 0
        if saved is None:
            self.commit(0)

    def commit(self, next_plan):
        assert next_plan >= self.next_plan
        for key, array in self.arrays.items():
            array.flush()
            with (self.root / (key + '.npy')).open('rb') as file:
                os.fsync(file.fileno())
        temporary = self.path.with_suffix('.tmp')
        with temporary.open('w') as file:
            json.dump(dict(identity=self.identity, next_plan=next_plan), file)
            file.flush()
            os.fsync(file.fileno())
        temporary.replace(self.path)
        directory = os.open(self.root, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
        self.next_plan = next_plan

    def remove(self):
        # Only these reader-created files, after successful output publication.
        self.path.unlink()
        for key, array in self.arrays.items():
            array._mmap.close()
            (self.root / (key + '.npy')).unlink()
        self.root.rmdir()
