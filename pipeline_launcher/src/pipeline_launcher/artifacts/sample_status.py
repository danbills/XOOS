"""Parse the pipeline's sample-status-report-latest.tsv to determine run outcome.

The pipeline writes a TSV with per-sample status (COMPLETED, FAILED,
CACHED, PARTIALLY_COMPLETED).  A run is considered successful when at
least one sample completed or was cached and none failed.

All public functions accept ``out_dir`` — the actual Nextflow ``--outdir``
where the pipeline wrote its outputs.  Do **not** pass ``analysis_dir``
and rely on ``output_dir()`` to derive the path; when the user sets
``--outdir`` equal to ``--analysis-dir`` there is no ``output/``
subdirectory.
"""

from __future__ import annotations

import csv
import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

_FAILURE_STATUSES = {"FAILED", "PARTIALLY_COMPLETED"}


@dataclass(frozen=True)
class SampleStatusReport:
    """Parsed sample status report — read once, query many times."""

    status_counts: Mapping[str, int] = MappingProxyType({})
    failed_sample_names: tuple[str, ...] = ()

    @property
    def succeeded(self) -> bool:
        """True when at least one sample succeeded and none failed."""
        if (
            self.status_counts.get("FAILED", 0) > 0
            or self.status_counts.get("PARTIALLY_COMPLETED", 0) > 0
        ):
            return False
        return (
            self.status_counts.get("COMPLETED", 0) > 0
            or self.status_counts.get("CACHED", 0) > 0
        )


def read_sample_status(out_dir: Path) -> SampleStatusReport:
    """Read sample-status-report-latest.tsv and return a parsed report.

    Returns an empty report when the file is missing or malformed.
    Rows with empty sample_id or status values are skipped with a warning.
    """
    report_path = out_dir / "sample-status-report-latest.tsv"
    if not report_path.exists():
        return SampleStatusReport()

    with report_path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or not {"sample_id", "status"} <= set(
            reader.fieldnames
        ):
            logging.warning(
                "Malformed sample status report (missing sample_id/status columns): %s",
                report_path,
            )
            return SampleStatusReport()

        counts: dict[str, int] = defaultdict(int)
        failed: list[str] = []
        for line_num, row in enumerate(reader, start=2):
            sample_id = row.get("sample_id") or ""
            status = row.get("status") or ""
            if not sample_id or not status:
                logging.warning(
                    "Skipping row %d with missing sample_id/status in %s",
                    line_num,
                    report_path,
                )
                continue
            counts[status] += 1
            if status in _FAILURE_STATUSES:
                failed.append(sample_id)

        return SampleStatusReport(
            status_counts=MappingProxyType(dict(counts)),
            failed_sample_names=tuple(failed),
        )


# ---------------------------------------------------------------------------
# Convenience wrappers — kept for backward compatibility with callers that
# don't need the full report object.
# ---------------------------------------------------------------------------


def parse_sample_status_report(out_dir: Path) -> dict[str, int]:
    """Parse sample-status-report-latest.tsv into counts by status."""
    return dict(read_sample_status(out_dir).status_counts)


def failed_samples(out_dir: Path) -> list[str]:
    """Return sample names with FAILED or PARTIALLY_COMPLETED status."""
    return list(read_sample_status(out_dir).failed_sample_names)


def is_pipeline_succeeded(out_dir: Path) -> bool:
    """Return True when at least one sample succeeded and none failed."""
    return read_sample_status(out_dir).succeeded
