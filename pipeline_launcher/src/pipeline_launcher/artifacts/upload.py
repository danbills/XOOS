"""Upload pipeline outputs and Nextflow artifacts using rclone.

Uploads the output directory in full, then selectively uploads
Nextflow logs and reports. Drops a CopyComplete.txt marker in the
destination to signal that the upload finished.
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

from pipeline_launcher.paths import nextflow_dir


def _rclone(
    action: str,
    src: Path | str,
    dst: Path | str,
    options: list[str],
    includes: list[str] | None = None,
) -> None:
    """Invoke rclone directly with the given action, options, and includes."""
    cmd: list[str] = ["rclone", action]
    cmd.extend(options)
    for inc in includes or []:
        cmd.extend(["--include", inc])
    cmd.extend([str(src), str(dst)])
    subprocess.check_call(cmd)


def upload(
    analysis_dir: Path, out_dir: Path, upload_dst: str, rclone_options: list[str]
) -> None:
    """Upload pipeline outputs and selected Nextflow artifacts.

    ``analysis_dir`` locates Nextflow runtime artifacts (logs).
    ``out_dir`` is the actual ``--outdir`` where pipeline outputs live.
    """
    upload_output = f"{upload_dst}/output"
    logging.info(f"Uploading {out_dir} to {upload_output}")
    _rclone("copy", out_dir, upload_output, rclone_options)

    upload_nf = f"{upload_dst}/nextflow"
    logging.info(f"Uploading {nextflow_dir(analysis_dir)} to {upload_nf}")
    _rclone(
        "copy",
        nextflow_dir(analysis_dir),
        upload_nf,
        rclone_options,
        includes=[
            "*.html",
            "nextflow.txt",
            ".nextflow.log",
            ".nextflow.log.*",
            "nextflow.log",
            "nextflow.log.*",
        ],
    )

    # Drop a completion marker then copy it to the destination.
    copy_complete = out_dir / "CopyComplete.txt"
    copy_complete.touch()
    _rclone("copy", copy_complete, upload_dst, rclone_options)
