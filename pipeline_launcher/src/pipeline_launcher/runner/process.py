"""Subprocess execution with log teeing.

Runs a command with stdout and stderr merged into a single stream that is
written to both the console and a log file for later inspection.
"""

from __future__ import annotations

import codecs
import errno
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import IO


def _stream_pipe(process: subprocess.Popen, lf: IO[str]) -> None:
    """Tee process output from a pipe to the console and a log file."""
    if process.stdout is None:
        raise RuntimeError("Failed to capture process stdout")
    for line in process.stdout:
        sys.stdout.buffer.write(line)
        sys.stdout.flush()
        lf.write(line.decode())
        lf.flush()


def _stream_pty(master_fd: int, lf: IO[str]) -> None:
    """Tee process output from a PTY master fd to the console and a log file.

    Reads until the slave side of the PTY is closed (EIO), which happens
    when the subprocess exits and all slave file descriptors are released.

    Uses an incremental UTF-8 decoder so that multi-byte characters split
    across read boundaries are reassembled correctly.
    """
    decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")
    while True:
        try:
            chunk = os.read(master_fd, 4096)
        except OSError as exc:
            if exc.errno == errno.EIO:
                break
            raise
        sys.stdout.buffer.write(chunk)
        sys.stdout.flush()
        text = decoder.decode(chunk)
        if text:
            lf.write(text)
            lf.flush()
    # Flush any remaining bytes in the decoder buffer.
    trailing = decoder.decode(b"", final=True)
    if trailing:
        lf.write(trailing)
        lf.flush()


def run_with_log(
    cmd: list[str],
    cwd: Path,
    env: dict[str, str],
    log: Path,
    use_pty: bool = False,
) -> None:
    """Run a command, teeing output to both the console and a log file.

    When *use_pty* is True the subprocess is attached to a pseudoterminal so
    that programs which inspect the file descriptor (e.g. jansi inside
    Nextflow) see a real TTY and emit ANSI escape sequences.

    Raises subprocess.CalledProcessError on non-zero exit.
    Raises KeyboardInterrupt if the user cancels (via SIGINT or SIGTERM),
    after attempting a graceful shutdown of the child process.
    """
    log.parent.mkdir(parents=True, exist_ok=True)
    master_fd: int | None = None
    slave_fd: int | None = None
    with log.open("w") as lf:
        try:
            if use_pty:
                master_fd, slave_fd = os.openpty()
                try:
                    process = subprocess.Popen(
                        cmd,
                        cwd=cwd,
                        env=env,
                        stdout=slave_fd,
                        stderr=slave_fd,
                        stdin=subprocess.DEVNULL,
                        close_fds=True,
                    )
                finally:
                    # Close the slave side regardless of whether Popen
                    # succeeded — the child inherits its own copy.
                    os.close(slave_fd)
                    slave_fd = None
            else:
                process = subprocess.Popen(
                    cmd,
                    cwd=cwd,
                    env=env,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                )
            try:
                if use_pty:
                    _stream_pty(master_fd, lf)
                else:
                    _stream_pipe(process, lf)
                returncode = process.wait()
                if returncode != 0:
                    raise subprocess.CalledProcessError(returncode, cmd)
            except KeyboardInterrupt:
                logging.info("Cancelling...")
                process.terminate()
                try:
                    process.wait(timeout=120)
                except subprocess.TimeoutExpired:
                    logging.warning("Nextflow did not terminate in time, killing...")
                    process.kill()
                    process.wait(timeout=120)
                raise
        finally:
            if master_fd is not None:
                os.close(master_fd)
