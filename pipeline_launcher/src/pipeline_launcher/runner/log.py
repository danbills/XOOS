"""Log file rotation.

Rotates log files by appending numeric suffixes (.1, .2, etc.) to
existing files before creating a new one at the base path. This keeps
a history of previous runs without overwriting.
"""

from __future__ import annotations

from pathlib import Path


def _rotated(p: Path, index: int) -> Path:
    return p if index == 0 else p.parent / f"{p.name}.{index}"


def rotate(p: Path) -> Path:
    """Rotate the given file path, shifting existing files up by one suffix.

    Returns the original path (now free for writing).
    """
    oldest = 0
    while _rotated(p, oldest).exists():
        oldest += 1

    for i in range(oldest, 0, -1):
        _rotated(p, i - 1).rename(_rotated(p, i))

    return p
