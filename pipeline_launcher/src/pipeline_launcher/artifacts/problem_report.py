"""Problem report generation on pipeline failure.

Collects Nextflow logs, work directory command files, pipeline info,
and metrics into a timestamped zip for troubleshooting.
"""

from __future__ import annotations

import logging
import os
import time
import zipfile
from pathlib import Path
from typing import Generator

from pipeline_launcher.paths import nextflow_dir, work_dir_target


def _collect_all(directory: Path) -> Generator[tuple[Path, Path], None, None]:
    """Yield all files and symlinks in a directory (non-recursive)."""
    if not directory.exists():
        return
    for item in directory.glob("*"):
        if item.is_file() or item.is_symlink():
            yield (item, item.relative_to(directory.parent))


def _collect_work_dir(work_dir: Path) -> Generator[tuple[Path, Path], None, None]:
    """Recursively collect Nextflow command files and symlinks from the work dir."""
    if not work_dir.exists():
        return

    to_include = {
        ".command.begin",
        ".command.err",
        ".command.log",
        ".command.out",
        ".command.run",
        ".command.sh",
        ".command.trace",
        ".exitcode",
    }
    for item in work_dir.rglob("*"):
        if (item.is_file() and item.name in to_include) or item.is_symlink():
            yield (item, item.relative_to(work_dir.parent))


def _collect_metrics(directory: Path) -> Generator[tuple[Path, Path], None, None]:
    """Collect metrics report files.

    Uses a fixed ``metrics/`` prefix with preserved subdirectory structure,
    so zip entries are stable and unique regardless of whether ``out_dir``
    equals ``analysis_dir``.
    """
    if not directory.exists():
        return

    to_include = {
        "sample-status-report-latest.tsv",
        "metrics-report.tsv",
    }
    for item in directory.rglob("*"):
        if item.is_file() and item.name in to_include:
            yield (item, Path("metrics") / item.relative_to(directory))


def _create_zip(files: list[tuple[Path, Path]], zip_path: Path) -> Path:
    """Create a zip file, storing symlink targets as .symlink text files.

    All entries are placed under a top-level directory matching the zip
    stem so that extracting the archive keeps files contained.
    """
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    prefix = zip_path.stem
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for abs_path, rel_path in files:
            arc_name = f"{prefix}/{rel_path}"
            if abs_path.is_symlink():
                # Zip doesn't preserve symlinks portably; store the target
                # so intent can be reconstructed when inspecting the report.
                link_target = os.readlink(abs_path)
                zf.writestr(f"{arc_name}.symlink", link_target)
            else:
                zf.write(abs_path, arc_name)
    return zip_path


def create_problem_report(analysis_dir: Path, out_dir: Path) -> Path:
    """Create a problem report zip with logs and artifacts for troubleshooting.

    ``analysis_dir`` locates Nextflow runtime artifacts (logs, work dir).
    ``out_dir`` is the actual ``--outdir`` where pipeline outputs live
    (pipeline_info, sample-status-report, metrics-report).  These may
    differ when the user explicitly sets ``--outdir``.
    """
    files: list[tuple[Path, Path]] = []
    files.extend(_collect_all(nextflow_dir(analysis_dir)))
    files.extend(_collect_work_dir(work_dir_target(analysis_dir)))
    files.extend(_collect_all(out_dir / "pipeline_info"))
    files.extend(_collect_metrics(out_dir))

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    zip_path = analysis_dir / "problem_reports" / f"problem_report_{timestamp}.zip"

    logging.error(
        f"Failure in one or more samples, problem report available: {zip_path}"
    )
    return _create_zip(files, zip_path)
