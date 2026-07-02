"""Parallel rclone copy with file partitioning.

Partitions files by size across N workers using a greedy bin-packing
approach (largest-first), then runs rclone copy concurrently for each
partition. This saturates network bandwidth better than a single
rclone process for large transfers with many files.
"""

from __future__ import annotations

import asyncio
import atexit
import logging
import shutil
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


def human_data_size(size: float | int) -> str:
    """Format a byte count as a human-readable string (binary units)."""
    for unit in ["B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB"]:
        if size < 1024.0:
            break
        size /= 1024.0
    return f"{size:.2f} {unit}"


def human_data_rate(rate: float | int) -> str:
    return human_data_size(rate) + "/s"


@dataclass
class FilePartition:
    files: list[str] = field(default_factory=list)
    size: int = 0


def determine_file_partitions(
    files: list[tuple[int, str]], num_partitions: int
) -> list[FilePartition]:
    """Partition files across N buckets using largest-first greedy bin packing.

    This balances total size across partitions so no single rclone
    worker is bottlenecked by a few large files.
    """
    files.sort(key=lambda x: x[0], reverse=True)
    partitions = [FilePartition() for _ in range(num_partitions)]

    for size, filepath in files:
        smallest = min(partitions, key=lambda p: p.size)
        smallest.files.append(filepath)
        smallest.size += size

    return partitions


async def _run_cmd(log_prefix: Path, cmd: list[str]) -> tuple[Path, Path]:
    logging.debug(f"Running command: {' '.join(cmd)}")
    stdout_file = log_prefix.with_suffix(".stdout.txt")
    stderr_file = log_prefix.with_suffix(".stderr.txt")
    with stdout_file.open("w") as stdout_f, stderr_file.open("w") as stderr_f:
        process = await asyncio.create_subprocess_exec(
            *cmd, stdout=stdout_f, stderr=stderr_f
        )
        await process.communicate()
        if process.returncode != 0:
            logging.error(f"Failed with return code {process.returncode}")
            logging.error(stdout_file.read_text())
            logging.error(stderr_file.read_text())
            raise RuntimeError(f"rclone failed with return code {process.returncode}")
    return stdout_file, stderr_file


async def _run_rclone(
    input_path: str,
    partition_file: Path,
    destination: str,
    log_prefix: Path,
    additional_args: list[str],
) -> None:
    cmd = [
        "rclone",
        "copy",
        input_path,
        destination,
        "--files-from",
        str(partition_file),
    ]
    cmd.extend(additional_args)
    await _run_cmd(log_prefix, cmd)


def _write_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fp:
        for line in lines:
            fp.write(f"{line}\n")


def _parse_files_from_output(output: Path) -> list[tuple[int, str]]:
    files = []
    for line in output.read_text().splitlines():
        size_str, filepath = line.split(";", 1)
        files.append((int(size_str), filepath))
    return files


async def _copy_async(
    src: str,
    dst: str,
    num_partitions: int,
    work_dir: Path,
    includes: list[str] | None,
    extra_args: list[str],
) -> None:
    # List files with sizes using rclone lsf.
    cmd = ["rclone", "lsf", "--format", "sp", "--recursive", "--files-only", src]
    if includes:
        cmd.extend(f"--include={inc}" for inc in includes)
    cmd.extend(extra_args)

    stdout, _ = await _run_cmd(work_dir / "rclone-lsf", cmd)
    files = _parse_files_from_output(stdout)
    if not files:
        logging.warning("No files found.")
        return

    partitions = determine_file_partitions(files, num_partitions)
    total_size = sum(size for size, _ in files)
    logging.info(f"Total size of all files: {human_data_size(total_size)}")
    logging.info(f"Total number of files: {len(files)}")

    tasks = []
    for index, partition in enumerate(partitions):
        if not partition.files:
            continue
        logging.debug(f"Size of partition {index}: {partition.size}")
        partition_file = work_dir / f"partition_{index}" / "files.txt"
        _write_lines(partition_file, partition.files)
        log_file = work_dir / f"partition_{index}" / "rclone-copy"
        tasks.append(_run_rclone(src, partition_file, dst, log_file, extra_args))

    start_time = time.time()
    await asyncio.gather(*tasks)
    total_time = time.time() - start_time

    logging.info(f"Total time taken: {total_time:.2f} seconds")
    logging.info(f"Transfer rate: {human_data_rate(total_size / total_time)}")


def rclone_copy(
    src: str,
    dst: str,
    num_partitions: int = 10,
    work_dir: Path | None = None,
    includes: list[str] | None = None,
    extra_args: list[str] | None = None,
) -> None:
    """Run a parallel rclone copy from src to dst.

    Partitions the source files across num_partitions concurrent rclone
    processes for higher throughput.
    """
    if work_dir is None:
        work_dir = Path(tempfile.mkdtemp(prefix="rclone_parallel_"))
        atexit.register(shutil.rmtree, work_dir, ignore_errors=True)

    if not work_dir.exists():
        work_dir.mkdir(parents=True)

    asyncio.run(
        _copy_async(src, dst, num_partitions, work_dir, includes, extra_args or [])
    )
