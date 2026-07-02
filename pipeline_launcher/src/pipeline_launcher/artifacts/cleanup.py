"""Work directory cleanup after pipeline execution.

Supports three policies:
  - delete-if-succeeded: remove only when all samples passed
  - delete-always: remove regardless of outcome
  - delete-never: preserve for debugging or resume
"""

from __future__ import annotations

import logging
import shutil
from enum import Enum
from pathlib import Path

from pipeline_launcher.paths import work_dir_target


class WorkDirDelete(Enum):
    DELETE_IF_SUCCEEDED = "delete-if-succeeded"
    DELETE_ALWAYS = "delete-always"
    DELETE_NEVER = "delete-never"

    def __str__(self) -> str:
        return self.value


def remove_work_dir(
    policy: WorkDirDelete, out_dir: Path, pipeline_succeeded: bool
) -> None:
    """Remove the Nextflow work directory according to the selected policy."""
    if policy == WorkDirDelete.DELETE_NEVER:
        logging.info("Preserving Nextflow work directory as per configuration.")
        return

    if policy == WorkDirDelete.DELETE_IF_SUCCEEDED and not pipeline_succeeded:
        logging.info(
            "Pipeline did not succeed; preserving work directory for debugging or resume."
        )
        return

    link = work_dir_target(out_dir)
    if link.exists() and link.is_symlink():
        actual = link.resolve()
    else:
        actual = link

    try:
        shutil.rmtree(actual)
    except OSError as e:
        logging.warning("Failed to remove work dir %s: %s", actual, e)
        # Keep the symlink so the location is still discoverable.
        return

    # Only remove the symlink after the target was successfully deleted.
    if link.exists() and link.is_symlink():
        link.unlink()
