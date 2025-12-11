"""Tests for SLURM-related functions using mocks."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

from rna_map_slurm.jobs.slurm import get_current_user, get_user_jobs


class TestGetCurrentUser:
    """Tests for get_current_user function."""

    def test_from_env_var(self) -> None:
        """Test getting user from USER env variable."""
        with patch.dict("os.environ", {"USER": "testuser"}):
            user = get_current_user()
            assert user == "testuser"

    def test_fallback_to_getpass(self) -> None:
        """Test fallback to getpass when env var not set."""
        with (
            patch.dict("os.environ", {}, clear=True),
            patch("getpass.getuser", return_value="fallbackuser"),
        ):
            user = get_current_user()
            assert user == "fallbackuser"


class TestGetUserJobs:
    """Tests for get_user_jobs function with mocked SLURM."""

    def test_get_jobs_success(self) -> None:
        """Test successful job query."""
        mock_output = "12345 batch job1 user1 R 0:10 1:00:00 1 node1"

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                stdout=mock_output,
                returncode=0,
            )

            jobs = get_user_jobs("testuser")

            assert len(jobs) == 1
            assert jobs[0]["JobID"] == "12345"
            assert jobs[0]["Name"] == "job1"

    def test_get_jobs_empty(self) -> None:
        """Test empty job queue."""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout="", returncode=0)

            jobs = get_user_jobs("testuser")

            assert jobs == []

    def test_get_jobs_multiple(self) -> None:
        """Test multiple jobs."""
        mock_output = """12345 batch job1 user1 R 0:10 1:00:00 1 node1
12346 batch job2 user1 PD 0:00 2:00:00 1 (Resources)"""

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=mock_output, returncode=0)

            jobs = get_user_jobs("testuser")

            assert len(jobs) == 2
            assert jobs[0]["Name"] == "job1"
            assert jobs[1]["Name"] == "job2"

    def test_get_jobs_slurm_error(self) -> None:
        """Test handling of SLURM errors."""
        with patch("subprocess.run") as mock_run:
            from subprocess import CalledProcessError

            mock_run.side_effect = CalledProcessError(1, "squeue", stderr="error")

            jobs = get_user_jobs("testuser")

            assert jobs == []


class TestSubmitJobs:
    """Tests for job submission with mocked sbatch."""

    def test_submit_job_called(self) -> None:
        """Test that sbatch is called correctly."""
        import pandas as pd

        from rna_map_slurm.cli.utils import submit_jobs

        df = pd.DataFrame({
            "job_path": ["/path/to/job1.sh", "/path/to/job2.sh"],
        })

        with patch("os.system") as mock_system:
            submit_jobs(df)

            assert mock_system.call_count == 2
            mock_system.assert_any_call("sbatch /path/to/job1.sh")
            mock_system.assert_any_call("sbatch /path/to/job2.sh")
