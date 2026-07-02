"""Shared command formatting and config generation for executor driver scripts.

Provides consistent quoted, backslash-continuation formatting across all
executor backends (AWS Batch, GCP Batch, Slurm), plus helpers for
generating Nextflow config snippets.

Values are single-quoted by default.  Values containing ``${`` are
double-quoted so the shell expands variable references (e.g.
``${PIPELINE_DIR}/main.nf`` when the driver sets ``PIPELINE_DIR`` at
runtime via ``mktemp -d``).
"""

from __future__ import annotations

import re
from importlib import resources
from pathlib import Path
from typing import TYPE_CHECKING

from jinja2 import Environment, FileSystemLoader, select_autoescape

if TYPE_CHECKING:
    from pipeline_launcher.executor.base import RunContext


def _single_quote(val: str) -> str:
    """Single-quote a value, escaping embedded single quotes."""
    return "'" + val.replace("'", "'\\''") + "'"


def _double_quote(val: str) -> str:
    """Double-quote a value, preserving shell variable expansion (``$var``, ``${var}``).

    Backslashes and double-quote characters in *val* are escaped so they
    are treated as literals; ``$`` is intentionally left unescaped so the
    shell expands variable references.
    """
    return '"' + val.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _quote(val: str) -> str:
    """Return *val* wrapped in the appropriate shell quotes.

    Values containing ``${`` are double-quoted so the shell expands the
    variable reference at runtime.  All other values are single-quoted.
    """
    if "${" in val:
        return _double_quote(val)
    return _single_quote(val)


def format_shell_command(args: list[str]) -> str:
    """Format an arbitrary shell command with backslash continuations.

    Same quoting rules as :func:`format_nf_command` but without
    prepending ``nextflow``.  The first element is the executable.
    Consecutive non-flag positional arguments (e.g. a subcommand such
    as ``run``) are placed on the same line as the executable, unquoted.
    """
    parts = _pair_args(args)
    if not parts:
        return ""

    # Collect the executable and any leading positional subcommand(s)
    # onto the first line, unquoted.  Stop at flags or at any positional
    # that looks like a path (contains "/" or "://"), which should be quoted.
    first_line: list[str] = []
    start = 0
    while start < len(parts) and not parts[start].startswith("-"):
        raw = parts[start].strip("'")
        if "/" in raw or "://" in raw:
            break
        first_line.append(raw)
        start += 1
    lines = [" ".join(first_line)]
    for part in parts[start:]:
        lines.append(f"    {part}")
    return " \\\n".join(lines)


def resolve_remote_path(
    local_path: Path | str | None,
    source_dir: Path | None,
    remote_base: str,
) -> str | None:
    """Map a local path to its remote equivalent inside the driver container.

    If *local_path* is relative to *source_dir*, returns the corresponding
    path under *remote_base*.  Otherwise returns the original path as a
    string (it may be a container-resident path).

    Returns ``None`` when *local_path* is ``None``.
    """
    if local_path is None:
        return None
    if source_dir is None:
        return str(local_path)
    try:
        rel = Path(local_path).relative_to(source_dir)
        return f"{remote_base}/{rel.as_posix()}"
    except ValueError:
        return str(local_path)


def build_nf_run_args(
    script: str,
    analysis_name_unique: str,
    *,
    config: str | None = None,
    extra_configs: list[str] | None = None,
    work_dir: str | None = None,
    output_dir: str | None = None,
    extra_args: list[str] | None = None,
) -> list[str]:
    """Build the Nextflow ``run`` argument list.

    Returns a list starting with ``"run"`` (no ``nextflow`` prefix).
    Callers should prepend ``"nextflow"`` as needed — either for subprocess
    use or before passing to :func:`format_shell_command` for formatted output.

    Assembles args in a consistent order: script, ``-name``, ``-resume``,
    user config, extra configs, ``-work-dir``, ``--outdir``, then any extra
    passthrough args.  ``-resume`` is always included.  ``--outdir`` is
    placed before ``extra_args`` so a passthrough ``--outdir`` overrides it
    (last-value-wins).
    """
    nf_run: list[str] = ["run", script, "-name", analysis_name_unique, "-resume"]
    if config is not None:
        nf_run.extend(["-c", config])
    for cfg_path in extra_configs or []:
        nf_run.extend(["-c", cfg_path])
    if work_dir is not None:
        nf_run.extend(["-work-dir", work_dir])
    if output_dir is not None:
        nf_run.extend(["--outdir", output_dir])
    if extra_args:
        nf_run.extend(extra_args)
    return nf_run


def _sanitize_label(value: str) -> str:
    """Sanitize a value for use as a cloud resource label.

    GCP labels must be lowercase, max 63 chars, only ``[a-z0-9_-]``.
    AWS tags are more permissive but we normalize for consistency.
    """
    sanitized = re.sub(r"[^a-z0-9_-]", "-", value.lower())[:63]
    return sanitized.strip("-")


def _make_jinja_env() -> Environment:
    """Return a Jinja2 :class:`~jinja2.Environment` for the bundled templates.

    Auto-escaping is explicitly disabled via ``select_autoescape([])``
    because every template in this package renders a shell script, not
    HTML or XML.  Enabling HTML auto-escaping would corrupt shell syntax
    by escaping ``&``, ``<``, ``>``, and quotes.  The explicit call makes
    the intent clear and resolves the SonarQube
    "Disabling auto-escaping is security-sensitive" finding.
    """
    templates = Path(str(resources.files("pipeline_launcher") / "templates"))
    return Environment(
        loader=FileSystemLoader(str(templates)),
        autoescape=select_autoescape([]),
    )


def render_labels_config(user: str, analysis_name: str) -> str:
    """Render a Nextflow config snippet that sets ``process.resourceLabels``.
    Uses the deterministic analysis_name (no hash).
    """
    env = _make_jinja_env()
    template = env.get_template("resource_labels.config.jinja")
    return template.render(
        user=_sanitize_label(user),
        analysis_name=_sanitize_label(analysis_name),
    )


def _pair_args(args: list[str]) -> list[str]:
    """Group flags with their values and quote all values.

    Values are single-quoted unless they contain ``${``, in which case
    they are double-quoted to allow shell variable expansion.
    """
    parts: list[str] = []
    i = 0
    while i < len(args):
        arg = args[i]
        if (
            arg.startswith("-")
            and i + 1 < len(args)
            and not args[i + 1].startswith("-")
        ):
            # Flag with a value — keep on same line.
            parts.append(f"{arg} {_quote(args[i + 1])}")
            i += 2
        elif arg.startswith("-"):
            # Bare flag (e.g. -resume).
            parts.append(arg)
            i += 1
        else:
            # Positional value (e.g. script path).
            parts.append(_quote(arg))
            i += 1
    return parts


def extract_work_dir_from_args(args: list[str]) -> str | None:
    """Extract the value of ``-work-dir`` from a Nextflow argument list.

    Handles both the two-token form (``-work-dir VALUE``) and the
    equals form (``-work-dir=VALUE``).  Returns ``None`` when the flag
    is absent.
    """
    for i, arg in enumerate(args):
        if arg == "-work-dir" and i + 1 < len(args):
            return args[i + 1]
        if arg.startswith("-work-dir="):
            return arg.split("=", 1)[1]
    return None
