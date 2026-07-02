"""Analysis name generation and sequencing run ID detection.

Generates stable, collision-resistant names for Nextflow runs by combining
a human-readable prefix (auto-detected run ID, user override, or date-based
fallback) with a short hash suffix. The hash incorporates the current
timestamp and full argv to reduce accidental collisions when multiple
runs launch near-simultaneously.
"""

from __future__ import annotations

import hashlib
import re
import sys
from datetime import datetime

_NF_NAME_INVALID = re.compile(r"[^a-z0-9_-]+")


def _sanitize_for_nextflow_name(name: str) -> str:
    r"""Return *name* made safe for use as a Nextflow ``-name`` value.

    Nextflow enforces ``^[a-z](?:[a-z\d]|[-_](?=[a-z\d])){0,79}$``.
    This helper:

    * Lowercases the entire string.
    * Replaces every character outside ``[a-z0-9_-]`` with ``_``.
    * Collapses consecutive separators (``_`` or ``-``) into a single ``_``.
    * Strips leading and trailing separator characters.
    """
    lowered = name.lower()
    replaced = _NF_NAME_INVALID.sub("_", lowered)
    collapsed = re.sub(r"[_-]+", "_", replaced)
    return collapsed.strip("_-")


# Matches sequencing run IDs like:
#   250801_XYZ-HTP_02_pt-015_ABY22R03C09_cycle01
#   20251105_PT-353_E2E-0103_PT-353_72000371R_Q1_R1
_RUN_ID_PATTERN = re.compile(r"\d{6,8}_[^/]*?_(?:cycle\d+|R\d+)")


def auto_detect_run_id(candidates: list[str]) -> str | None:
    """Detect a sequencing run ID from candidate strings.

    Returns the match only when exactly one distinct ID is found,
    avoiding ambiguity when multiple runs appear in the arguments.
    """
    matches: set[str] = set()
    for token in candidates:
        for m in _RUN_ID_PATTERN.finditer(token):
            matches.add(m.group(0))

    if len(matches) == 1:
        return next(iter(matches))
    return None


def generate_analysis_name(
    name_override: str | None = None,
    analysis_name_override: str | None = None,
) -> str:
    """Return the deterministic, hash-free base name for the work directory path.

    Unlike :func:`generate_analysis_name_unique`, no random suffix is appended.
    The returned value is identical across re-invocations with the same
    inputs, so Nextflow can resume from the same work directory without
    pinning ``--analysis-name``.

    Priority mirrors :func:`generate_analysis_name_unique`:
      1. ``analysis_name_override`` — returned as-is when provided
      2. Auto-detected sequencing run ID from ``sys.argv``
      3. ``name_override``
      4. Date-based fallback (``YYYYMMDD_run``)
    """
    if analysis_name_override is not None:
        return analysis_name_override
    return (
        auto_detect_run_id(sys.argv) or name_override or f"{datetime.now():%Y%m%d}_run"
    )


def generate_analysis_name_unique(
    name_override: str | None = None,
    analysis_name_override: str | None = None,
) -> str:
    """Generate a stable, unique name for this analysis (with hash suffix).

    Priority:
      1. analysis_name_override — sanitized and used as the base when provided
      2. Auto-detected run ID from sys.argv
      3. name_override
      4. Date-based fallback

    Appends 4 hex characters from a hash of datetime + argv to reduce
    collisions.
    """
    base = _sanitize_for_nextflow_name(
        analysis_name_override
        or auto_detect_run_id(sys.argv)
        or name_override
        or f"{datetime.now():%Y%m%d}_run"
    )
    hash_input = f"{datetime.now().isoformat()}\n" + "\n".join(sys.argv)
    suffix = hashlib.sha256(hash_input.encode("utf-8")).hexdigest()[:4]
    max_nextflow_name_length = 80
    max_prefix_length = max_nextflow_name_length - len(suffix) - 1
    truncated = f"xoos_{base}"[:max_prefix_length].rstrip("_-")
    return f"{truncated}_{suffix}"
