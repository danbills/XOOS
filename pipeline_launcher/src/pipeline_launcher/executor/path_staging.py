"""Allowlist-driven path staging for cloud executors.

Certain passthrough flags (``-c``, ``--input``, ``--samplesheet``) may
carry local file paths that the driver container cannot access.  This
module classifies and rewrites those values before the Nextflow command
is assembled:

* **config flags** (``-c``): local files are uploaded to the cloud staging
  area.  The driver downloads them to ``/pipeline/conf/`` before Nextflow
  starts.  The value in the arg list is rewritten to the container-local
  path.  Cloud URIs are added to the fetch list so the driver downloads
  them directly.  Opaque values (e.g. container-baked paths) are passed
  through unchanged.

* **data flags** (``--input``, ``--samplesheet``): local files are uploaded
  to the cloud staging area and the value is rewritten to the cloud URI
  (``s3://`` or ``gs://``), which Nextflow can consume directly.  Cloud
  URIs are passed through unchanged.  Opaque values are passed through
  unchanged.

Cloud URIs from the wrong cloud provider are rejected with a user-facing
error.

Usage::

    from pipeline_launcher.executor.path_staging import stage_allowlisted_args

    result = stage_allowlisted_args(context.extra_args, stager)
    if isinstance(result, str):
        return SubmitResult(succeeded=False, message=result, exit_code=1)
    context.extra_args = result.rewritten_args
    # result.driver_fetch_files → fed to the Jinja template
    # result.staged_data_files  → audit log of data uploads
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

# Flags whose values are Nextflow config files.  Local paths are staged to
# /pipeline/conf/ and added to driver_fetch_files; cloud URIs are fetched
# there by the driver.  The value in extra_args is rewritten to the
# container path.
CONFIG_FLAGS: frozenset[str] = frozenset({"-c"})

# Flags whose values are data files.  Local paths are uploaded and the value
# is rewritten to the cloud URI.  Cloud URIs are passed through unchanged.
DATA_FLAGS: frozenset[str] = frozenset({"--input", "--samplesheet"})

# Combined set used for quick membership checks.
_ALL_STAGED_FLAGS: frozenset[str] = CONFIG_FLAGS | DATA_FLAGS


# ---------------------------------------------------------------------------
# Flag value types
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class S3Uri:
    """An ``s3://`` URI."""

    uri: str


@dataclass(frozen=True)
class GcsUri:
    """A ``gs://`` URI."""

    uri: str


@dataclass(frozen=True)
class LocalFile:
    """A path that exists as a file on the local filesystem.

    ``path`` is already resolved (absolute, symlinks followed).
    """

    path: Path


@dataclass(frozen=True)
class OpaqueValue:
    """Any other value: a container-baked path, a pipeline reference, etc.

    Passed through unchanged by all stagers without any upload.
    """

    raw: str


FlagValue = S3Uri | GcsUri | LocalFile | OpaqueValue


def classify_value(value: str) -> FlagValue:
    """Classify a raw flag value into the appropriate :data:`FlagValue` type."""
    if value.startswith("s3://"):
        return S3Uri(value)
    if value.startswith("gs://"):
        return GcsUri(value)
    p = Path(value)
    if p.is_file():
        return LocalFile(p.resolve())
    return OpaqueValue(value)


# ---------------------------------------------------------------------------
# Staging result types
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DriverFetchFile:
    """A file the driver container must download before Nextflow starts.

    Produced by config-flag staging (``-c``).  The driver downloads
    ``source_uri`` and saves it to ``container_path``.
    """

    source_uri: str
    container_path: str


@dataclass(frozen=True)
class StagedDataFile:
    """Record of a data file that was staged to cloud storage.

    Produced by data-flag staging (``--input``, ``--samplesheet``).
    ``cloud_uri`` is the value written into the rewritten arg list.
    When the original value was already a cloud URI,
    ``original_value == cloud_uri``.
    """

    flag: str
    original_value: str
    cloud_uri: str


@dataclass(frozen=True)
class StagingResult:
    """The outcome of a successful :func:`stage_allowlisted_args` call."""

    rewritten_args: list[str]
    driver_fetch_files: list[DriverFetchFile]
    staged_data_files: list[StagedDataFile]


# ---------------------------------------------------------------------------
# Stager protocol
# ---------------------------------------------------------------------------


class Stager(Protocol):
    """Protocol for executor-specific staging implementations.

    Each cloud executor provides a concrete stager that knows its native
    cloud, upload helpers, and destination paths.
    """

    def stage_config(
        self, flag: str, value: FlagValue
    ) -> tuple[str, DriverFetchFile | None]:
        """Stage a config-flag value.

        Returns ``(rewritten_arg_value, fetch_file_or_None)``.
        :class:`OpaqueValue` inputs must return ``(raw, None)`` without
        performing any upload.
        Raises :exc:`ValueError` with a user-facing message on invalid input.
        """
        ...

    def stage_data(
        self, flag: str, value: FlagValue
    ) -> tuple[str, StagedDataFile | None]:
        """Stage a data-flag value.

        Returns ``(rewritten_arg_value, staged_record_or_None)``.
        :class:`OpaqueValue` inputs must return ``(raw, None)`` without
        performing any upload.
        Raises :exc:`ValueError` with a user-facing message on invalid input.
        """
        ...


# ---------------------------------------------------------------------------
# Core walker
# ---------------------------------------------------------------------------


def stage_allowlisted_args(
    args: list[str],
    stager: Stager,
) -> StagingResult | str:
    """Walk *args*, staging values for allowlisted flags via *stager*.

    Handles both ``--flag value`` (two-token) and ``--flag=value``
    (single-token) forms.  All other tokens are copied through unchanged.

    Returns a :class:`StagingResult` on success or an error string on
    failure so the caller can surface a
    :class:`~pipeline_launcher.executor.base.SubmitResult`.
    """
    out: list[str] = []
    fetch_files: list[DriverFetchFile] = []
    data_files: list[StagedDataFile] = []
    i = 0
    while i < len(args):
        token = args[i]

        # Handle --flag=value inline form.
        if "=" in token and token.startswith("-"):
            flag, _, raw = token.partition("=")
            if flag in _ALL_STAGED_FLAGS:
                value = classify_value(raw)
                try:
                    if flag in CONFIG_FLAGS:
                        new_val, record = stager.stage_config(flag, value)
                    else:
                        new_val, record = stager.stage_data(flag, value)
                except ValueError as exc:
                    return str(exc)
                _accumulate(record, fetch_files, data_files)
                out.append(f"{flag}={new_val}")
            else:
                out.append(token)
            i += 1
            continue

        # Handle --flag value (two-token) form.
        if token in _ALL_STAGED_FLAGS and i + 1 < len(args):
            raw = args[i + 1]
            value = classify_value(raw)
            try:
                if token in CONFIG_FLAGS:
                    new_val, record = stager.stage_config(token, value)
                else:
                    new_val, record = stager.stage_data(token, value)
            except ValueError as exc:
                return str(exc)
            _accumulate(record, fetch_files, data_files)
            out.append(token)
            out.append(new_val)
            i += 2
            continue

        out.append(token)
        i += 1

    return StagingResult(
        rewritten_args=out,
        driver_fetch_files=fetch_files,
        staged_data_files=data_files,
    )


def _accumulate(
    record: DriverFetchFile | StagedDataFile | None,
    fetch_files: list[DriverFetchFile],
    data_files: list[StagedDataFile],
) -> None:
    if isinstance(record, DriverFetchFile):
        fetch_files.append(record)
    elif isinstance(record, StagedDataFile):
        data_files.append(record)
