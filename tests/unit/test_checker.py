"""Tests for job checker functionality."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from rna_map_slurm.cli.check_cmd import check_jobs
from rna_map_slurm.jobs.checker import CheckReport, JobChecker, JobStatus


class TestJobStatus:
    """Tests for JobStatus dataclass."""

    def test_success_status(self) -> None:
        """Test creating a success status."""
        status = JobStatus(job_name="test-job", status="success")

        assert status.job_name == "test-job"
        assert status.status == "success"
        assert status.reason is None

    def test_failed_status_with_reason(self) -> None:
        """Test creating a failed status with reason."""
        status = JobStatus(
            job_name="test-job",
            status="failed",
            reason="Job exceeded wall time limit",
            error_line="CANCELLED AT ... DUE TO TIME LIMIT",
        )

        assert status.status == "failed"
        assert status.reason == "Job exceeded wall time limit"


class TestCheckReport:
    """Tests for CheckReport dataclass."""

    def test_empty_report(self) -> None:
        """Test empty report properties."""
        report = CheckReport(total=0)

        assert report.total == 0
        assert report.success_count == 0
        assert report.failure_count == 0
        assert report.missing_count == 0

    def test_report_with_jobs(self) -> None:
        """Test report with various job statuses."""
        succeeded = [JobStatus("job1", "success"), JobStatus("job2", "success")]
        failed = [JobStatus("job3", "failed", reason="time_limit")]
        missing = [JobStatus("job4", "missing")]

        report = CheckReport(total=4, succeeded=succeeded, failed=failed, missing=missing)

        assert report.success_count == 2
        assert report.failure_count == 1
        assert report.missing_count == 1

    def test_summary(self) -> None:
        """Test summary generation."""
        report = CheckReport(
            total=3,
            succeeded=[JobStatus("job1", "success")],
            failed=[JobStatus("job2", "failed")],
            missing=[JobStatus("job3", "missing")],
        )

        summary = report.summary()

        assert "Total jobs: 3" in summary
        assert "Succeeded:  1" in summary
        assert "Failed:     1" in summary
        assert "Missing:    1" in summary


class TestJobChecker:
    """Tests for JobChecker class."""

    def test_finds_output_files(self, temp_dir: Path) -> None:
        """Test finding output files in directory."""
        # Create some output files
        (temp_dir / "job-001.out").write_text("completed successfully")
        (temp_dir / "job-002.out").write_text("completed successfully")

        checker = JobChecker(job_dir=temp_dir)
        files = checker._find_output_files()

        assert len(files) == 2
        assert "job-001" in files
        assert "job-002" in files

    def test_detects_time_limit_error(self, temp_dir: Path) -> None:
        """Test detection of time limit exceeded error."""
        output_file = temp_dir / "job-001.out"
        output_file.write_text(
            "Running job...\n"
            "Processing data...\n"
            "slurmstepd: error: *** JOB 13356506 ON c2007 CANCELLED AT 2026-01-04T15:37:27 DUE TO TIME LIMIT ***\n"
        )

        checker = JobChecker(job_dir=temp_dir)
        success, error_type, error_line = checker._parse_output_file(output_file)

        assert success is False
        assert error_type == "time_limit"
        assert "TIME LIMIT" in error_line

    def test_detects_oom_error(self, temp_dir: Path) -> None:
        """Test detection of out of memory error."""
        output_file = temp_dir / "job-001.out"
        output_file.write_text(
            "Running job...\n"
            "slurmstepd: error: Detected 1 oom-kill event(s) in StepId=123.0\n"
        )

        checker = JobChecker(job_dir=temp_dir)
        success, error_type, error_line = checker._parse_output_file(output_file)

        assert success is False
        assert error_type == "out_of_memory"

    def test_detects_success(self, temp_dir: Path) -> None:
        """Test detection of successful job."""
        output_file = temp_dir / "job-001.out"
        output_file.write_text(
            "Running job...\n"
            "Processing complete.\n"
            "Job finished successfully.\n"
        )

        checker = JobChecker(job_dir=temp_dir)
        success, error_type, error_line = checker._parse_output_file(output_file)

        assert success is True
        assert error_type is None

    def test_check_all_with_mixed_results(self, temp_dir: Path) -> None:
        """Test check_all with mixed job results."""
        # Create output files with different statuses
        (temp_dir / "job-001.out").write_text("Job completed successfully")
        (temp_dir / "job-002.out").write_text(
            "slurmstepd: error: *** JOB CANCELLED DUE TO TIME LIMIT ***"
        )

        checker = JobChecker(job_dir=temp_dir)
        report = checker.check_all()

        assert report.total == 2
        assert report.success_count == 1
        assert report.failure_count == 1

    def test_check_with_jobs_csv(self, temp_dir: Path) -> None:
        """Test checking against jobs.csv."""
        # Create jobs.csv
        jobs_csv = temp_dir / "jobs.csv"
        df = pd.DataFrame({
            "job_path": ["jobs/job-001.sh", "jobs/job-002.sh", "jobs/job-003.sh"],
            "job_type": ["rna-map", "rna-map", "rna-map"],
            "job_requirement": ["", "", ""],
        })
        df.to_csv(jobs_csv, index=False)

        # Create output for only some jobs
        job_dir = temp_dir / "jobs"
        job_dir.mkdir()
        (job_dir / "job-001.out").write_text("completed")
        (job_dir / "job-002.out").write_text("completed")
        # job-003 has no output file

        checker = JobChecker(job_dir=job_dir, jobs_csv=jobs_csv)
        report = checker.check_all()

        assert report.total == 3
        assert report.success_count == 2
        assert report.missing_count == 1

    def test_suggest_fixes_time_limit(self, temp_dir: Path) -> None:
        """Test fix suggestions for time limit errors."""
        (temp_dir / "rna-map-001.out").write_text(
            "slurmstepd: error: *** JOB CANCELLED DUE TO TIME LIMIT ***"
        )

        checker = JobChecker(job_dir=temp_dir)
        report = checker.check_all()
        suggestions = checker.suggest_fixes(report)

        assert len(suggestions) >= 1
        assert any("time" in s.lower() for s in suggestions)

    def test_suggest_fixes_oom(self, temp_dir: Path) -> None:
        """Test fix suggestions for OOM errors."""
        (temp_dir / "rna-map-001.out").write_text(
            "slurmstepd: error: Detected 1 oom-kill event(s)"
        )

        checker = JobChecker(job_dir=temp_dir)
        report = checker.check_all()
        suggestions = checker.suggest_fixes(report)

        assert len(suggestions) >= 1
        assert any("memory" in s.lower() for s in suggestions)


class TestCheckJobsCommand:
    """Tests for the check-jobs CLI command."""

    def test_no_job_dir(self, temp_dir: Path) -> None:
        """Test error when job directory doesn't exist."""
        runner = CliRunner()

        result = runner.invoke(check_jobs, ["--job-dir", str(temp_dir / "nonexistent")])

        assert result.exit_code == 1
        assert "not found" in result.output

    def test_empty_job_dir(self, temp_dir: Path) -> None:
        """Test with empty job directory."""
        runner = CliRunner()

        result = runner.invoke(check_jobs, ["--job-dir", str(temp_dir)])

        assert result.exit_code == 0
        assert "No jobs found" in result.output

    def test_successful_jobs(self, temp_dir: Path) -> None:
        """Test with all successful jobs."""
        (temp_dir / "job-001.out").write_text("completed")
        (temp_dir / "job-002.out").write_text("completed")

        runner = CliRunner()
        result = runner.invoke(check_jobs, ["--job-dir", str(temp_dir)])

        assert result.exit_code == 0
        assert "Succeeded" in result.output

    def test_failed_jobs_exit_code(self, temp_dir: Path) -> None:
        """Test that failed jobs result in non-zero exit code."""
        (temp_dir / "job-001.out").write_text(
            "slurmstepd: error: *** JOB CANCELLED DUE TO TIME LIMIT ***"
        )

        runner = CliRunner()
        result = runner.invoke(check_jobs, ["--job-dir", str(temp_dir)])

        assert result.exit_code == 1
        assert "Failed" in result.output

    def test_fix_flag_shows_suggestions(self, temp_dir: Path) -> None:
        """Test that --fix flag shows suggestions."""
        (temp_dir / "rna-map-001.out").write_text(
            "slurmstepd: error: *** JOB CANCELLED DUE TO TIME LIMIT ***"
        )

        runner = CliRunner()
        result = runner.invoke(check_jobs, ["--job-dir", str(temp_dir), "--fix"])

        assert "Suggested Fixes" in result.output

    def test_verbose_flag(self, temp_dir: Path) -> None:
        """Test that --verbose flag shows more details."""
        (temp_dir / "job-001.out").write_text(
            "slurmstepd: error: *** JOB CANCELLED DUE TO TIME LIMIT ***"
        )

        runner = CliRunner()
        result = runner.invoke(check_jobs, ["--job-dir", str(temp_dir), "--verbose"])

        # Should show the error line in verbose mode
        assert "Error Line" in result.output or "TIME LIMIT" in result.output
