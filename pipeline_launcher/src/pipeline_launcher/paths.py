"""Analysis directory layout helpers.

The launcher organizes all pipeline artifacts beneath --analysis-dir:
    analysis_dir/output/          Pipeline results produced by Nextflow (set as --outdir if not specified)
    analysis_dir/nextflow/        Nextflow runtime artifacts (logs, trace, etc.)
    analysis_dir/nextflow/work    Symlink to the actual Nextflow work directory
    analysis_dir/stage/           Batch driver scripts and staging artifacts
    analysis_dir/problem_reports/ Zipped logs collected on failure
    analysis_dir/.nextflow/       Nextflow state directory (for resume)
"""

from __future__ import annotations

from pathlib import Path


def output_dir(analysis_dir: Path) -> Path:
    return (analysis_dir / "output").resolve()


def nextflow_dir(analysis_dir: Path) -> Path:
    return (analysis_dir / "nextflow").resolve()


def stage_dir(analysis_dir: Path) -> Path:
    return (analysis_dir / "stage").resolve()


def work_dir_target(analysis_dir: Path) -> Path:
    """Canonical symlink location for the work directory under analysis_dir."""
    return nextflow_dir(analysis_dir) / "work"


# Nextflow state directory (for resume)
def state_dir(analysis_dir: Path) -> Path:
    return (analysis_dir / ".nextflow").resolve()
