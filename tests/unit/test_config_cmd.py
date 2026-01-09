"""Tests for generate-example-config command."""

from __future__ import annotations

from pathlib import Path

import yaml
from click.testing import CliRunner

from rna_map_slurm.cli.config_cmd import _build_commented_config, generate_example_config


class TestBuildCommentedConfig:
    """Tests for _build_commented_config function."""

    def test_returns_commented_map(self) -> None:
        """Test that function returns a CommentedMap."""
        config = _build_commented_config()

        # Should be able to access keys
        assert "fastq_chunks" in config
        assert "paths" in config
        assert "slurm_options" in config

    def test_contains_all_main_sections(self) -> None:
        """Test that all main config sections are present."""
        config = _build_commented_config()

        expected_sections = [
            "fastq_chunks",
            "paths",
            "construct_options",
            "tasks_per_job",
            "slurm_options",
        ]
        for section in expected_sections:
            assert section in config, f"Missing section: {section}"

    def test_paths_has_all_keys(self) -> None:
        """Test that paths section has all required keys."""
        config = _build_commented_config()
        paths = config["paths"]

        expected_keys = ["log", "jobs", "submits", "inputs", "tmp", "seq_path"]
        for key in expected_keys:
            assert key in paths, f"Missing path key: {key}"

    def test_slurm_options_has_job_types(self) -> None:
        """Test that slurm_options has all job type configurations."""
        config = _build_commented_config()
        slurm_opts = config["slurm_options"]

        expected_job_types = [
            "default",
            "split-fastq",
            "trim-galore",
            "demultiplex",
            "int-demultiplex",
            "join-fastq-files",
            "rna-map",
            "rna-map-combine",
            "int-demultiplex-rna-map",
            "int-demultiplex-rna-map-combine",
        ]
        for job_type in expected_job_types:
            assert job_type in slurm_opts, f"Missing job type: {job_type}"


class TestGenerateExampleConfigCommand:
    """Tests for the generate-example-config CLI command."""

    def test_generates_file(self, temp_dir: Path) -> None:
        """Test that command generates a YAML file."""
        runner = CliRunner()
        output_path = temp_dir / "test_config.yml"

        result = runner.invoke(generate_example_config, ["--output", str(output_path)])

        assert result.exit_code == 0
        assert output_path.exists()

    def test_output_is_valid_yaml(self, temp_dir: Path) -> None:
        """Test that generated file is valid YAML."""
        runner = CliRunner()
        output_path = temp_dir / "test_config.yml"

        runner.invoke(generate_example_config, ["--output", str(output_path)])

        # Should parse without error
        with open(output_path) as f:
            config = yaml.safe_load(f)

        assert isinstance(config, dict)
        assert "fastq_chunks" in config

    def test_output_contains_comments(self, temp_dir: Path) -> None:
        """Test that generated file contains comments."""
        runner = CliRunner()
        output_path = temp_dir / "test_config.yml"

        runner.invoke(generate_example_config, ["--output", str(output_path)])

        content = output_path.read_text()
        assert "#" in content  # Should have comment markers

    def test_default_output_path(self) -> None:
        """Test that default output path is used."""
        runner = CliRunner()

        with runner.isolated_filesystem():
            result = runner.invoke(generate_example_config)

            assert result.exit_code == 0
            assert Path("example_config.yml").exists()

    def test_displays_success_message(self, temp_dir: Path) -> None:
        """Test that success message is displayed."""
        runner = CliRunner()
        output_path = temp_dir / "test_config.yml"

        result = runner.invoke(generate_example_config, ["--output", str(output_path)])

        assert "Generated example configuration file" in result.output
