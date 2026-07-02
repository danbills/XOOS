"""Slurm executor — renders a driver script and submits via sbatch.

The driver script re-invokes the launcher inside the driver environment
so the actual Nextflow run happens on the allocated compute node with
the correct scratch paths and Singularity cache.
"""

from __future__ import annotations

import json
import logging
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from pipeline_launcher.config.model import SlurmAttributes, SlurmConfig
from pipeline_launcher.executor.base import RunContext, SubmitResult
from pipeline_launcher.executor.formatting import format_shell_command, _make_jinja_env
from pipeline_launcher.location import to_location
from pipeline_launcher.paths import stage_dir
from pipeline_launcher.runner.log import rotate
from cloudpathlib import CloudPath


def _build_sbatch_directives(
    job_name: str,
    comment: str,
    attrs: SlurmAttributes,
) -> list[tuple[str, str]]:
    """Return an ordered list of (directive-flag, shell-quoted value) pairs.

    The required ``--job-name`` directive is always first.  Optional
    directives from ``attrs.driver_options`` are appended in iteration
    order; entries whose values are ``None`` or empty are skipped.
    """
    directives: list[tuple[str, str]] = [
        ("--job-name", shlex.quote(job_name)),
    ]
    if comment:
        directives.append(("--comment", shlex.quote(comment)))
    for flag, value in attrs.driver_options.items():
        if value:
            directives.append((f"--{flag}", shlex.quote(value)))
    return directives


def render_slurm_driver_script(
    job_name: str,
    command: list[str],
    comment: str,
    config: SlurmConfig,
    launcher_cwd: str | None = None,
    analysis_dir: str | None = None,
    out_dir: str | None = None,
) -> str:
    """Render a Slurm batch script from the Jinja template."""
    attrs = config.attributes

    data = {
        "sbatchDirectives": _build_sbatch_directives(job_name, comment, attrs),
        "preamble": config.preamble,
        "launcherCwd": shlex.quote(launcher_cwd) if launcher_cwd else None,
        "command": format_shell_command(command),
        "outDir": (
            shlex.quote(out_dir or analysis_dir) if (out_dir or analysis_dir) else None
        ),
    }

    env = _make_jinja_env()
    template = env.get_template("slurm_driver.jinja")
    return template.render(data=data)


def _convert_args_to_driver_command(
    analysis_name: str,
    extra_args: list[str],
) -> list[str]:
    """Convert current sys.argv into a command for the Slurm driver script.

    Preserves the original --env value and adds --driver-mode so the
    re-invoked launcher uses the local executor.  The --name flag is
    replaced with --analysis-name to lock the run name.

    ``extra_args`` supplies the passthrough arguments (everything after
    ``--``) instead of extracting them from ``sys.argv``.  This is
    necessary because the initial invocation may inject ``--input``
    into the passthrough list after generating a run sheet — an
    argument that never appeared in the original ``sys.argv``.

    Relative paths in launcher-managed arguments are canonicalized to
    absolute paths since the driver re-invokes from a different cwd.
    """
    # Run sheet flags are consumed by the initial invocation which
    # generates the run sheet and injects --input into passthrough.
    # They must not be forwarded to the driver or the sheet would be
    # generated a second time.
    run_sheet_flags: set[str] = {
        "--run-name",
        "--run-type",
        "--file-type",
    }
    launcher_args_needing_canonicalization: set[str] = {
        "--pipeline-script",
        "--resources-base",
    }
    # Run sheet flags whose values are paths that also need stripping.
    run_sheet_path_flags: set[str] = {
        "--run-dir",
        "--samplesheet",
    }
    launcher_args: list[str] = []
    skip_next = False

    # Skip the first two args (the launcher executable and the "run" subcommand).
    # Only process launcher args (before --); passthrough args come from extra_args.
    in_passthrough = False
    for i, arg in enumerate(sys.argv[2:], 2):
        if skip_next:
            skip_next = False
            continue

        if arg == "--":
            in_passthrough = True
            continue

        # Skip everything after -- in sys.argv; we use extra_args instead.
        if in_passthrough:
            continue

        if arg == "--name":
            skip_next = True
            continue

        # Strip run sheet flags — the initial invocation already
        # generated the sheet and added --input to passthrough.
        if arg in run_sheet_flags or arg in run_sheet_path_flags:
            skip_next = True
            continue

        if arg in launcher_args_needing_canonicalization:
            skip_next = True
            launcher_args.append(arg)
            value = sys.argv[i + 1]
            # Preserve cloud URIs as-is; canonicalize local paths.
            loc = to_location(value)
            if isinstance(loc, CloudPath):
                launcher_args.append(value)
            else:
                launcher_args.append(str(loc))
            continue

        if arg == "--env":
            # Only canonicalize the env value if it's a file path
            value = sys.argv[i + 1]
            if Path(value).is_file():
                skip_next = True
                launcher_args.append(arg)
                launcher_args.append(str(Path(value).resolve()))
                continue

        launcher_args.append(arg)

    # Use the absolute path to the xoos binary so the driver script
    # works even when the venv is not activated inside the SLURM job.
    argv0 = sys.argv[0]
    p = Path(argv0)
    if p.is_absolute():
        xoos_bin = str(p)
    else:
        resolved = shutil.which(argv0)
        xoos_bin = resolved if resolved else str(p.resolve())
    command = [xoos_bin, "run"]
    command.extend(launcher_args)
    command.extend(["--driver-mode"])
    command.extend(["--analysis-name", analysis_name])
    command.append("--")
    command.extend(extra_args)

    return command


def _write_executable_script(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(content)
    path.chmod(0o755)


def _submit_sbatch(script: Path) -> str:
    job_id = subprocess.check_output(
        ["sbatch", "--parsable", str(script)],
        cwd=script.parent,
    )
    return job_id.decode().strip()


class SlurmExecutor:
    """Generate a Slurm driver script and submit it via sbatch."""

    def submit(self, context: RunContext) -> SubmitResult:
        config = context.config
        assert isinstance(config, SlurmConfig)
        assert isinstance(context.out_dir, Path), (
            "Slurm executor requires a local --outdir, got a cloud URI: "
            f"{context.out_dir}"
        )

        if not context.project:
            logging.error("--project must be provided when submitting an HPC job.")
            return SubmitResult(
                succeeded=False,
                message="--project is required for Slurm submission.",
                exit_code=1,
            )

        comment = json.dumps({"project": context.project})

        driver_command = _convert_args_to_driver_command(
            context.analysis_name,
            context.extra_args,
        )

        script_content = render_slurm_driver_script(
            context.analysis_name,
            driver_command,
            comment,
            config,
            launcher_cwd=str(context.launcher_cwd),
            analysis_dir=str(context.analysis_dir),
            out_dir=str(context.out_dir),
        )

        script_path = rotate(stage_dir(context.analysis_dir) / "driver.sh")
        _write_executable_script(script_path, script_content)

        job_id = _submit_sbatch(script_path)
        logging.info(f"Submitted batch job {job_id} with name {context.analysis_name}")

        return SubmitResult(
            succeeded=True,
            message=f"Submitted batch job {job_id}",
            job_id=job_id,
        )
