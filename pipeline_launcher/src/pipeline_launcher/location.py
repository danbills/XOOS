"""Location type for values that can be either a local path or a cloud URI.

Cloud URIs (s3://, gs://, az://) are represented as ``cloudpathlib.CloudPath``
subclasses (S3Path, GSPath, AzureBlobPath).  Local filesystem paths are
``pathlib.Path``.  ``LocationPath`` is the union type used throughout the
launcher wherever either form is valid.

Use ``to_location()`` to parse a raw string into the appropriate type.
Use ``isinstance(value, CloudPath)`` to branch on cloud vs local.
"""

from __future__ import annotations

from pathlib import Path

from cloudpathlib import AnyPath, CloudPath

# Union type for any value that can be a local path or cloud URI.
LocationPath = Path | CloudPath


def to_location(value: str | Path | CloudPath) -> LocationPath:
    """Parse *value* into a ``CloudPath`` or a resolved local ``Path``.

    - Already-typed ``CloudPath`` and ``Path`` values are returned as-is.
    - Strings beginning with ``s3://``, ``gs://``, or ``az://`` are
      parsed by ``cloudpathlib.AnyPath`` into the matching ``CloudPath``.
    - All other strings are resolved to absolute ``Path`` objects.
    """
    if isinstance(value, CloudPath):
        return value
    if isinstance(value, Path):
        return value
    # Raw string — dispatch on scheme.
    if value.startswith(("s3://", "gs://", "az://")):
        return AnyPath(value)  # type: ignore[return-value]
    return Path(value).resolve()
