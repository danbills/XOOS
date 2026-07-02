"""Local executor — runs Nextflow directly on the current machine.

Handles optional file locking, log rotation, problem report generation
on failure, rclone upload, and work directory cleanup.
"""

from __future__ import annotations

import logging
import shlex
import subprocess
import sys
from contextlib import nullcontext
from pathlib import Path

from pipeline_launcher.artifacts.cleanup import WorkDirDelete, remove_work_dir
from pipeline_launcher.artifacts.problem_report import create_problem_report
from pipeline_launcher.artifacts.sample_status import (
    SampleStatusReport,
    read_sample_status,
)
from pipeline_launcher.artifacts.upload import upload

from pipeline_launcher.config.model import LocalConfig, SlurmConfig
from pipeline_launcher.executor.base import RunContext, SubmitResult
from pipeline_launcher.executor.formatting import extract_work_dir_from_args
from pipeline_launcher.locking import file_lock

from pipeline_launcher.paths import nextflow_dir
from pipeline_launcher.runner.command import (
    build_nextflow_command,
    build_nextflow_env,
    create_work_dir,
    create_work_dir_symlink,
)
from pipeline_launcher.runner.log import rotate
from pipeline_launcher.runner.process import run_with_log


def _invoke_callback(
    callback: str | None,
    subcommand: str,
    analysis_name: str,
    analysis_dir: Path | None = None,
) -> None:
    """Best-effort callback invocation around pipeline lifecycle events."""
    if callback is None:
        return
    try:
        cmd = shlex.split(callback)
        full_cmd = cmd + [subcommand, "--analysis-name", analysis_name]
        if analysis_dir is not None and subcommand in ("complete", "fail"):
            full_cmd.extend(["--analysis-dir", str(analysis_dir)])
        subprocess.run(full_cmd, check=False)
    except (OSError, ValueError) as e:
        logging.warning("Callback failed for subcommand %s: %s", subcommand, e)


def _resolve_config_paths(
    config: LocalConfig | SlurmConfig, launcher_cwd: Path
) -> list[Path]:
    """Resolve config file paths, making relative paths absolute against launcher_cwd.

    Paths prefixed with ``@nextflow_config/`` are resolved from the bundled
    ``nextflow_config/`` directory inside the pipeline_launcher package.
    """
    from pipeline_launcher.config.loader import resolve_config_path

    return [resolve_config_path(p, launcher_cwd) for p in config.config_files]


def _upload_if_configured(context: RunContext) -> None:
    """Upload pipeline outputs if an upload destination is configured."""
    if context.upload_dst is not None:
        rclone_opts = shlex.split(context.rclone_options)
        upload(context.analysis_dir, context.out_dir, context.upload_dst, rclone_opts)


def _cleanup_work_dir(context: RunContext, succeeded: bool) -> None:
    """Remove the Nextflow work directory according to the configured policy."""
    wdd = WorkDirDelete(context.work_dir_delete)
    remove_work_dir(wdd, context.analysis_dir, succeeded)


def _log_sample_failures(report: SampleStatusReport) -> None:
    """Log which samples failed or partially completed."""
    names = report.failed_sample_names
    if names:
        listing = ", ".join(names)
        logging.error(
            "Nextflow exited successfully but %d sample(s) failed: %s",
            len(names),
            listing,
        )
    else:
        logging.error("Nextflow exited successfully but some samples failed.")


class LocalExecutor:
    """Execute Nextflow directly on the local machine."""

    def submit(self, context: RunContext) -> SubmitResult:
        config = context.config
        assert isinstance(config, (LocalConfig, SlurmConfig))
        assert isinstance(context.out_dir, Path), (
            "local/Slurm executor requires a local --outdir, got a cloud URI: "
            f"{context.out_dir}"
        )

        # For SlurmConfig in driver mode, we have scratch/singularity settings.
        driver = config.driver if isinstance(config, SlurmConfig) else None

        # When scratch is configured, create a work dir there and symlink it
        # from the analysis dir for discoverability. Without scratch, Nextflow
        # defaults to a subdirectory of the launch dir. A user-supplied
        # -work-dir in passthrough args takes precedence over both.
        if extract_work_dir_from_args(context.extra_args) is not None:
            work = None
        elif driver:
            work = create_work_dir(driver, context.launcher_cwd)
        else:
            work = None
        if work is not None:
            create_work_dir_symlink(context.analysis_dir, work)

        # Nextflow is launched from this directory so its logs and trace
        # files land under the analysis dir.
        nf_launch_dir = nextflow_dir(context.analysis_dir)
        nf_launch_dir.mkdir(parents=True, exist_ok=True)

        config_files = _resolve_config_paths(config, context.launcher_cwd)

        cmd = build_nextflow_command(
            definition=context.definition,
            analysis_name=context.analysis_name_unique,
            config_files=config_files,
            extra_args=context.extra_args,
            work_dir=work,
        )
        cmd_env = build_nextflow_env(driver, context.singularity_cache)
        cmd_env.setdefault("NXF_FILE_ROOT", str(context.launcher_cwd))

        lock_ctx = file_lock(context.lock_path) if context.lock_path else nullcontext()

        with lock_ctx:
            nextflow_log = rotate(nf_launch_dir / "nextflow.log")
            logging.info("Output directory: %s", context.out_dir)
            try:
                _invoke_callback(context.callback, "start", context.analysis_name)
                run_with_log(
                    cmd,
                    nf_launch_dir,
                    cmd_env,
                    nextflow_log,
                    use_pty=sys.stdout.isatty(),
                )

                report = read_sample_status(context.out_dir)
                cb_event = "complete" if report.succeeded else "fail"
                _invoke_callback(
                    context.callback,
                    cb_event,
                    context.analysis_name,
                    analysis_dir=context.out_dir,
                )

                if not report.succeeded:
                    _log_sample_failures(report)
                    create_problem_report(context.analysis_dir, context.out_dir)

            except subprocess.CalledProcessError as e:
                _invoke_callback(
                    context.callback,
                    "fail",
                    context.analysis_name,
                    analysis_dir=context.out_dir,
                )
                logging.error(
                    "Pipeline failed with exit code %d. See log: %s",
                    e.returncode,
                    nextflow_log,
                )
                create_problem_report(context.analysis_dir, context.out_dir)
                _upload_if_configured(context)
                _cleanup_work_dir(context, False)

                return SubmitResult(
                    succeeded=False,
                    message=f"Pipeline failed with exit code {e.returncode}",
                    exit_code=e.returncode,
                )

        _upload_if_configured(context)
        _cleanup_work_dir(context, report.succeeded)

        if not report.succeeded:
            return SubmitResult(
                succeeded=False,
                message="Nextflow exited successfully but some samples failed.",
                exit_code=1,
            )

        return SubmitResult(succeeded=True, message="Pipeline completed successfully.")
