"""Run sheet generation for the xoos-nf-core pipeline.

Generates a single-row run sheet CSV from CLI arguments, written to
the output directory so it is preserved alongside pipeline results.

Cloud paths (s3://, gs://, az://) are accepted for ``run_dir`` and
``samplesheet`` — filesystem validation is skipped for those.
"""

from __future__ import annotations

import csv
import io
import logging
from dataclasses import dataclass
from pathlib import Path

from cloudpathlib import CloudPath

from pipeline_launcher.location import LocationPath, to_location

HEADER = ["run_name", "run_type", "file_type", "run_dir", "samplesheet"]


@dataclass
class RunSheetRow:
    """A single row in the nf-core xoos run sheet.

    ``run_dir`` and ``samplesheet`` accept either local ``Path`` objects
    or cloud URI strings (``s3://…``, ``gs://…``).
    """

    run_name: str
    run_type: str
    file_type: str
    run_dir: LocationPath
    samplesheet: LocationPath

    def validate(self) -> list[str]:
        """Return a list of validation error messages (empty if valid)."""
        errors: list[str] = []
        if not self.run_name or " " in self.run_name:
            errors.append(
                f"run_name must be a non-empty string without spaces, "
                f"got '{self.run_name}'"
            )
        if not self.run_type or not self.run_type.strip():
            errors.append("run_type must be a non-empty string")
        if not self.file_type or not self.file_type.strip():
            errors.append("file_type must be a non-empty string")

        # Filesystem checks only apply to local paths.
        run_dir_loc = (
            to_location(str(self.run_dir))
            if not isinstance(self.run_dir, (Path, CloudPath))
            else self.run_dir
        )
        samplesheet_loc = (
            to_location(str(self.samplesheet))
            if not isinstance(self.samplesheet, (Path, CloudPath))
            else self.samplesheet
        )

        if not isinstance(run_dir_loc, CloudPath):
            if not Path(str(run_dir_loc)).is_dir():
                errors.append(f"run_dir does not exist: {self.run_dir}")

        if not isinstance(samplesheet_loc, CloudPath):
            if not Path(str(samplesheet_loc)).is_file():
                errors.append(f"samplesheet does not exist: {self.samplesheet}")

        return errors


def _resolve_path(value: str | LocationPath) -> str:
    """Resolve a path to an absolute string, preserving cloud URIs."""
    loc = to_location(str(value)) if isinstance(value, str) else value
    return str(loc)


def generate_run_sheet(row: RunSheetRow, out_dir: Path) -> Path:
    """Write a single-row run sheet CSV and return its path.

    The file is written to ``<out_dir>/run_sheet.csv``.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "run_sheet.csv"

    errors = row.validate()
    if errors:
        for err in errors:
            logging.error("Run sheet validation: %s", err)
        raise SystemExit(1)

    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(HEADER)
        writer.writerow(
            [
                row.run_name,
                row.run_type,
                row.file_type,
                _resolve_path(row.run_dir),
                _resolve_path(row.samplesheet),
            ]
        )

    logging.info("Generated run sheet: %s", path)
    return path


def generate_run_sheet_content(row: RunSheetRow) -> str:
    """Return a run sheet CSV as a string.

    Used by cloud executors that upload the content directly to object
    storage without writing to disk.
    """
    errors = row.validate()
    if errors:
        for err in errors:
            logging.error("Run sheet validation: %s", err)
        raise SystemExit(1)

    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(HEADER)
    writer.writerow(
        [
            row.run_name,
            row.run_type,
            row.file_type,
            _resolve_path(row.run_dir),
            _resolve_path(row.samplesheet),
        ]
    )
    return buf.getvalue()
