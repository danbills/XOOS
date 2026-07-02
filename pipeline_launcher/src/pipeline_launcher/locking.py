"""Cooperative file-based locking for serializing Nextflow runs.

Uses fcntl advisory locks to prevent concurrent pipeline executions
on the same output directory. The lock is non-blocking with periodic
polling to avoid busy-waiting.
"""

from __future__ import annotations

import fcntl
import logging
import time
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Any, Generator


def _try_lock(file: IO[Any]) -> bool:
    """Try to acquire an exclusive non-blocking lock."""
    try:
        fcntl.flock(file, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return True
    except BlockingIOError:
        return False


@contextmanager
def file_lock(
    path: Path, poll_interval_secs: int = 300
) -> Generator[IO[Any], None, None]:
    """Acquire a cooperative file lock, polling until available.

    poll_interval_secs controls how often to retry when the lock is
    held by another process. Defaults to 5 minutes to avoid log spam
    on long-running pipelines.
    """
    with path.open("w") as file:
        while not _try_lock(file):
            logging.info(f"Waiting for lock on {path}")
            time.sleep(poll_interval_secs)
        try:
            yield file
        finally:
            fcntl.flock(file, fcntl.LOCK_UN)
