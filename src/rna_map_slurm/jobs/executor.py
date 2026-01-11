"""Job execution using submitit for SLURM array submissions."""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import pandas as pd

from rna_map_slurm.models.config import SlurmOptions
from rna_map_slurm.utils.logging import get_logger

log = get_logger("jobs.executor")


@dataclass
class SubmitResult:
    """Result of batch job submission."""

    submitted: int = 0
    failed: int = 0
    job_ids: list[str] = field(default_factory=list)
    failed_jobs: list[str] = field(default_factory=list)
    elapsed_seconds: float = 0.0


@dataclass
class ArrayJobSpec:
    """Specification for a SLURM job array."""

    job_type: str
    job_scripts: list[Path]
    slurm_options: SlurmOptions
    array_size: int = 0

    def __post_init__(self) -> None:
        self.array_size = len(self.job_scripts)


def parse_slurm_time_to_minutes(time_str: str) -> int:
    """Convert SLURM time format to minutes.

    Args:
        time_str: Time in format HH:MM:SS, D-HH:MM:SS, or minutes.

    Returns:
        Time in minutes.
    """
    if time_str.isdigit():
        return int(time_str)

    if "-" in time_str:
        days, rest = time_str.split("-")
        parts = rest.split(":")
        hours = int(parts[0])
        minutes = int(parts[1]) if len(parts) > 1 else 0
        return int(days) * 24 * 60 + hours * 60 + minutes

    parts = time_str.split(":")
    hours = int(parts[0])
    minutes = int(parts[1]) if len(parts) > 1 else 0
    return hours * 60 + minutes


def parse_memory_to_mb(mem_str: str) -> int:
    """Convert memory string to MB.

    Args:
        mem_str: Memory in format like 2GB, 4000MB, or 4000.

    Returns:
        Memory in MB.
    """
    mem_str = mem_str.upper().strip()

    if mem_str.endswith("GB") or mem_str.endswith("G"):
        value = float(mem_str.rstrip("GB").rstrip("G"))
        return int(value * 1024)
    elif mem_str.endswith("MB") or mem_str.endswith("M"):
        return int(float(mem_str.rstrip("MB").rstrip("M")))
    else:
        return int(mem_str)


class JobExecutor:
    """Execute jobs using submitit with SLURM array support."""

    def __init__(
        self,
        log_folder: str = "logs/submitit",
        max_array_size: int = 1000,
        max_concurrent: int = 1000,
    ) -> None:
        """Initialize the job executor.

        Args:
            log_folder: Folder for submitit logs.
            max_array_size: Maximum jobs per array (SLURM limit).
            max_concurrent: Maximum concurrent array tasks.
        """
        self.log_folder = Path(log_folder)
        self.max_array_size = max_array_size
        self.max_concurrent = max_concurrent

    def _create_executor(
        self, slurm_options: SlurmOptions, array_parallelism: int | None = None
    ) -> Any:
        """Create a submitit executor with the given options.

        Args:
            slurm_options: SLURM configuration options.
            array_parallelism: Max concurrent array tasks.

        Returns:
            Configured submitit AutoExecutor.
        """
        import submitit

        executor = submitit.AutoExecutor(folder=str(self.log_folder))

        timeout_min = parse_slurm_time_to_minutes(slurm_options.time)
        mem_mb = parse_memory_to_mb(slurm_options.mem_per_cpu)

        params = {
            "timeout_min": timeout_min,
            "mem_gb": mem_mb / 1024,
            "cpus_per_task": slurm_options.cpus_per_task,
            "slurm_job_name": slurm_options.name,
        }

        if array_parallelism:
            params["slurm_array_parallelism"] = array_parallelism

        if slurm_options.extra_header_cmds:
            params["slurm_additional_parameters"] = {
                "setup": slurm_options.extra_header_cmds
            }

        executor.update_parameters(**params)
        return executor

    def submit_array(
        self,
        job_type: str,
        tasks: list[Callable[[], Any]],
        slurm_options: SlurmOptions,
    ) -> list[Any]:
        """Submit tasks as a SLURM job array.

        Args:
            job_type: Type of job for logging.
            tasks: List of callable tasks to execute.
            slurm_options: SLURM configuration.

        Returns:
            List of submitit Job objects.
        """
        if not tasks:
            log.warning(f"No tasks to submit for {job_type}")
            return []

        log.info(f"Submitting {len(tasks)} tasks for {job_type} as array")

        executor = self._create_executor(
            slurm_options, array_parallelism=self.max_concurrent
        )

        jobs = executor.map_array(lambda fn: fn(), tasks)
        log.info(f"Submitted array job for {job_type}: {len(jobs)} tasks")

        return jobs

    def submit_scripts_as_array(
        self,
        job_type: str,
        scripts: list[Path],
        slurm_options: SlurmOptions,
        job_dir: Path,
        dependency_job_ids: list[str] | None = None,
    ) -> SubmitResult:
        """Submit shell scripts as a SLURM job array.

        Creates a task file and array wrapper script, then submits.

        Args:
            job_type: Type of job.
            scripts: List of shell script paths.
            slurm_options: SLURM configuration.
            job_dir: Directory for job files.
            dependency_job_ids: Optional list of job IDs that must complete first.

        Returns:
            SubmitResult with submission statistics.
        """
        import time

        if not scripts:
            return SubmitResult()

        start = time.time()
        job_ids: list[str] = []
        failed_jobs: list[str] = []

        # Split into chunks if exceeding max array size
        chunks = [
            scripts[i : i + self.max_array_size]
            for i in range(0, len(scripts), self.max_array_size)
        ]

        for chunk_idx, chunk in enumerate(chunks):
            try:
                # Create task file listing all scripts
                task_file = job_dir / f"{job_type}-tasks-{chunk_idx}.txt"
                with open(task_file, "w") as f:
                    for script in chunk:
                        f.write(f"bash {script}\n")

                # Create array wrapper script
                array_script = job_dir / f"{job_type}-array-{chunk_idx}.sh"
                array_indices = f"0-{len(chunk) - 1}"
                if self.max_concurrent:
                    array_indices += f"%{self.max_concurrent}"

                script_content = self._generate_array_script(
                    slurm_options=slurm_options,
                    array_indices=array_indices,
                    task_file=task_file,
                    job_dir=job_dir,
                    job_type=job_type,
                    chunk_idx=chunk_idx,
                )

                with open(array_script, "w") as f:
                    f.write(script_content)

                # Submit the array job with dependencies
                success, job_id = self._submit_array_script(
                    array_script, dependency_job_ids
                )

                if success:
                    job_ids.append(job_id or "unknown")
                    log.info(
                        f"Submitted array {job_type}-{chunk_idx}: "
                        f"{len(chunk)} tasks, job_id={job_id}"
                    )
                else:
                    failed_jobs.extend([str(s) for s in chunk])
                    log.error(f"Failed to submit array {job_type}-{chunk_idx}")

            except Exception as e:
                log.error(f"Error submitting array {job_type}-{chunk_idx}: {e}")
                failed_jobs.extend([str(s) for s in chunk])

        elapsed = time.time() - start

        return SubmitResult(
            submitted=len(scripts) - len(failed_jobs),
            failed=len(failed_jobs),
            job_ids=job_ids,
            failed_jobs=failed_jobs,
            elapsed_seconds=elapsed,
        )

    def _generate_array_script(
        self,
        slurm_options: SlurmOptions,
        array_indices: str,
        task_file: Path,
        job_dir: Path,
        job_type: str,
        chunk_idx: int,
    ) -> str:
        """Generate a SLURM array job script.

        Args:
            slurm_options: SLURM configuration.
            array_indices: Array index specification (e.g., "0-99%50").
            task_file: Path to file listing tasks.
            job_dir: Directory for output files.
            job_type: Job type name.
            chunk_idx: Chunk index for naming.

        Returns:
            Shell script content.
        """
        output_path = job_dir / f"{job_type}-array-{chunk_idx}-%a.out"
        error_path = job_dir / f"{job_type}-array-{chunk_idx}-%a.err"

        script = f"""#!/bin/bash
#SBATCH --job-name={job_type}-array
#SBATCH --array={array_indices}
#SBATCH --time={slurm_options.time}
#SBATCH --mem-per-cpu={slurm_options.mem_per_cpu}
#SBATCH --cpus-per-task={slurm_options.cpus_per_task}
#SBATCH --output={output_path}
#SBATCH --error={error_path}

{slurm_options.extra_header_cmds}

# Get the task command from the task file
# SLURM_ARRAY_TASK_ID is 0-indexed
TASK_LINE=$((SLURM_ARRAY_TASK_ID + 1))
TASK=$(sed -n "${{TASK_LINE}}p" {task_file})

echo "Running task $SLURM_ARRAY_TASK_ID: $TASK"
eval $TASK
EXIT_CODE=$?

echo "Task $SLURM_ARRAY_TASK_ID completed with exit code $EXIT_CODE"
exit $EXIT_CODE
"""
        return script

    def _submit_array_script(
        self,
        script_path: Path,
        dependency_job_ids: list[str] | None = None,
    ) -> tuple[bool, str | None]:
        """Submit an array job script using sbatch.

        Args:
            script_path: Path to the array job script.
            dependency_job_ids: Optional list of job IDs that must complete first.

        Returns:
            Tuple of (success, job_id).
        """
        import re

        try:
            cmd = ["sbatch"]

            # Add dependency if specified
            if dependency_job_ids:
                dep_str = ":".join(dependency_job_ids)
                cmd.extend(["--dependency", f"afterok:{dep_str}"])

            cmd.append(str(script_path))

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )

            if result.returncode == 0:
                match = re.search(r"Submitted batch job (\d+)", result.stdout)
                job_id = match.group(1) if match else None
                return (True, job_id)
            else:
                log.error(f"sbatch failed: {result.stderr.strip()}")
                return (False, None)

        except FileNotFoundError:
            log.error("sbatch not found - not on a SLURM cluster")
            return (False, None)
        except Exception as e:
            log.error(f"Error submitting {script_path}: {e}")
            return (False, None)


def submit_jobs_as_arrays(
    df: pd.DataFrame,
    slurm_options_map: dict[str, SlurmOptions],
    job_dir: Path = Path("jobs"),
    max_concurrent: int = 1000,
) -> SubmitResult:
    """Submit all jobs from a DataFrame as SLURM arrays.

    Args:
        df: DataFrame with job_type and job_path columns.
        slurm_options_map: Mapping of job type to SLURM options.
        job_dir: Base directory for job files.
        max_concurrent: Maximum concurrent array tasks.

    Returns:
        Combined SubmitResult for all job types.
    """
    executor = JobExecutor(max_concurrent=max_concurrent)

    total_result = SubmitResult()

    for job_type, group in df.groupby("job_type"):
        scripts = [Path(p) for p in group["job_path"]]
        options = slurm_options_map.get(
            str(job_type), SlurmOptions(name=str(job_type))
        )

        type_job_dir = job_dir / str(job_type)
        type_job_dir.mkdir(parents=True, exist_ok=True)

        result = executor.submit_scripts_as_array(
            job_type=str(job_type),
            scripts=scripts,
            slurm_options=options,
            job_dir=type_job_dir,
        )

        total_result.submitted += result.submitted
        total_result.failed += result.failed
        total_result.job_ids.extend(result.job_ids)
        total_result.failed_jobs.extend(result.failed_jobs)
        total_result.elapsed_seconds += result.elapsed_seconds

    return total_result
