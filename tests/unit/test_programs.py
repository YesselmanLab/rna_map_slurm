"""Tests for checking external program availability.

These tests verify that all required external programs and Python packages
are available for the rna_map_slurm workflow to function correctly.

Tests are organized into:
- Required CLI programs: Must be available for core functionality
- Optional CLI programs: Enhance functionality but not strictly required
- Required Python packages: Core dependencies that may need manual installation
"""

from __future__ import annotations

import shutil
import subprocess

import pytest

# External command-line programs required by the workflow
REQUIRED_CLI_PROGRAMS = [
    "sabre",      # Demultiplexing (io/demultiplex.py)
    "bowtie2",    # Alignment (used by rna_map)
    "gzip",       # Compression (io/demultiplex.py)
    "cat",        # File concatenation (tasks/basic.py, runner/fastq_ops.py)
]

# Optional CLI programs that enhance functionality
OPTIONAL_CLI_PROGRAMS = [
    "fastqc",       # Quality control
    "cutadapt",     # Adapter trimming
    "trim_galore",  # Trimming wrapper
    "seqkit",       # Sequence manipulation
]

# Core Python packages (always required)
CORE_PYTHON_PACKAGES = [
    ("fastqsplitter", "fastqsplitter"),
]

# Full dependency Python packages (required for full functionality)
FULL_DEPENDENCY_PACKAGES = [
    ("rna_map", "rna_map"),
    ("barcode_demultiplex", "barcode_demultiplex"),
]

# Check if full dependencies are available
try:
    import rna_map  # noqa: F401
    import barcode_demultiplex  # noqa: F401
    HAS_FULL_DEPENDENCIES = True
except ImportError:
    HAS_FULL_DEPENDENCIES = False

requires_full_deps = pytest.mark.skipif(
    not HAS_FULL_DEPENDENCIES,
    reason="Requires full dependencies (rna_map, barcode_demultiplex)",
)


class TestRequiredCLIPrograms:
    """Tests for required command-line programs."""

    @pytest.mark.parametrize("program", REQUIRED_CLI_PROGRAMS)
    def test_required_program_available(self, program: str) -> None:
        """Test that required CLI programs are available in PATH."""
        path = shutil.which(program)
        assert path is not None, (
            f"Required program '{program}' not found in PATH. "
            f"Please install it (see environment.yml for installation instructions)."
        )

    @pytest.mark.parametrize("program", REQUIRED_CLI_PROGRAMS)
    def test_required_program_executable(self, program: str) -> None:
        """Test that required CLI programs can be executed."""
        path = shutil.which(program)
        if path is None:
            pytest.skip(f"Program '{program}' not found in PATH")

        # Try to get version or help output to verify executability
        # Different programs have different flags, so we try common ones
        version_flags = ["--version", "-v", "--help", "-h"]
        executed = False

        for flag in version_flags:
            try:
                result = subprocess.run(
                    [program, flag],
                    capture_output=True,
                    timeout=10,
                )
                # Some programs return non-zero for --help, that's OK
                executed = True
                break
            except (subprocess.TimeoutExpired, FileNotFoundError):
                continue

        assert executed, f"Could not execute '{program}' with any common flag"


class TestOptionalCLIPrograms:
    """Tests for optional command-line programs."""

    @pytest.mark.parametrize("program", OPTIONAL_CLI_PROGRAMS)
    def test_optional_program_available(self, program: str) -> None:
        """Test that optional CLI programs are available (warns if missing)."""
        path = shutil.which(program)
        if path is None:
            pytest.skip(
                f"Optional program '{program}' not found in PATH. "
                f"Some features may be unavailable."
            )
        assert path is not None


class TestCorePythonPackages:
    """Tests for core Python packages that are always required."""

    @pytest.mark.parametrize("package_name,import_name", CORE_PYTHON_PACKAGES)
    def test_core_python_package_importable(
        self, package_name: str, import_name: str
    ) -> None:
        """Test that core Python packages can be imported."""
        try:
            __import__(import_name)
        except ImportError:
            pytest.fail(
                f"Core Python package '{package_name}' could not be imported. "
                f"Please install it: pip install {package_name}"
            )


class TestFullDependencyPythonPackages:
    """Tests for Python packages required for full functionality."""

    @requires_full_deps
    @pytest.mark.parametrize("package_name,import_name", FULL_DEPENDENCY_PACKAGES)
    def test_full_dependency_package_importable(
        self, package_name: str, import_name: str
    ) -> None:
        """Test that full dependency Python packages can be imported."""
        try:
            __import__(import_name)
        except ImportError:
            pytest.fail(
                f"Full dependency package '{package_name}' could not be imported. "
                f"Please install it: pip install {package_name}"
            )


class TestSabreProgram:
    """Specific tests for the sabre demultiplexer."""

    def test_sabre_available(self) -> None:
        """Test that sabre is available."""
        path = shutil.which("sabre")
        assert path is not None, (
            "sabre not found. Install via: "
            "brew install sabre (macOS) or conda install -c bioconda sabre (Linux)"
        )

    def test_sabre_pe_subcommand(self) -> None:
        """Test that sabre pe subcommand is available."""
        path = shutil.which("sabre")
        if path is None:
            pytest.skip("sabre not found in PATH")

        # sabre pe --help should work
        result = subprocess.run(
            ["sabre", "pe", "--help"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        # sabre may return non-zero for help, check if it ran at all
        assert result.returncode == 0 or "usage" in result.stderr.lower() or "usage" in result.stdout.lower(), (
            "sabre pe subcommand not working as expected"
        )


class TestBowtie2Program:
    """Specific tests for bowtie2 aligner."""

    def test_bowtie2_available(self) -> None:
        """Test that bowtie2 is available."""
        path = shutil.which("bowtie2")
        assert path is not None, (
            "bowtie2 not found. Install via: conda install -c bioconda bowtie2"
        )

    def test_bowtie2_build_available(self) -> None:
        """Test that bowtie2-build is available."""
        path = shutil.which("bowtie2-build")
        assert path is not None, (
            "bowtie2-build not found. This usually comes with bowtie2 installation."
        )

    def test_bowtie2_version(self) -> None:
        """Test that bowtie2 can report its version."""
        path = shutil.which("bowtie2")
        if path is None:
            pytest.skip("bowtie2 not found in PATH")

        result = subprocess.run(
            ["bowtie2", "--version"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        assert result.returncode == 0, "bowtie2 --version failed"
        assert "bowtie2" in result.stdout.lower(), "Unexpected bowtie2 version output"


@requires_full_deps
class TestRnaMapPackage:
    """Specific tests for the rna_map Python package."""

    def test_rna_map_importable(self) -> None:
        """Test that rna_map can be imported."""
        try:
            import rna_map  # noqa: F401
        except ImportError:
            pytest.fail(
                "rna_map package not found. Install via: "
                "pip install git+https://github.com/jyesselm/rna_map.git"
            )

    def test_rna_map_run_module(self) -> None:
        """Test that rna_map.run module is available."""
        try:
            import rna_map.run  # noqa: F401
        except ImportError:
            pytest.fail("rna_map.run module not available")

    def test_rna_map_mutation_histogram_module(self) -> None:
        """Test that rna_map.mutation_histogram module is available."""
        try:
            from rna_map.mutation_histogram import (  # noqa: F401
                get_mut_histos_from_pickle_file,
                merge_mut_histo_dicts,
                write_mut_histos_to_pickle_file,
            )
        except ImportError:
            pytest.fail("rna_map.mutation_histogram module or functions not available")

    def test_rna_map_parameters_module(self) -> None:
        """Test that rna_map.parameters module is available."""
        try:
            from rna_map.parameters import (  # noqa: F401
                get_preset_params,
                parse_parameters_from_file,
            )
        except ImportError:
            pytest.fail("rna_map.parameters module or functions not available")


class TestFastqsplitterPackage:
    """Specific tests for the fastqsplitter Python package."""

    def test_fastqsplitter_importable(self) -> None:
        """Test that fastqsplitter can be imported."""
        try:
            import fastqsplitter  # noqa: F401
        except ImportError:
            pytest.fail(
                "fastqsplitter package not found. Install via: "
                "pip install fastqsplitter"
            )

    def test_fastqsplitter_split_function(self) -> None:
        """Test that fastqsplitter.split_fastqs function is available."""
        try:
            from fastqsplitter import split_fastqs  # noqa: F401
        except ImportError:
            pytest.fail("fastqsplitter.split_fastqs function not available")


@requires_full_deps
class TestBarcodeDemultiplexPackage:
    """Specific tests for the barcode_demultiplex Python package."""

    def test_barcode_demultiplex_importable(self) -> None:
        """Test that barcode_demultiplex can be imported."""
        try:
            import barcode_demultiplex  # noqa: F401
        except ImportError:
            pytest.fail(
                "barcode_demultiplex package not found. Install via: "
                "pip install git+https://github.com/jyesselm/barcode_demultiplex.git"
            )

    def test_find_helix_barcodes_function(self) -> None:
        """Test that find_helix_barcodes function is available."""
        try:
            from barcode_demultiplex.demultiplex import find_helix_barcodes  # noqa: F401
        except ImportError:
            pytest.fail(
                "barcode_demultiplex.demultiplex.find_helix_barcodes not available"
            )


def check_all_programs() -> dict[str, bool]:
    """Check availability of all programs and return status dict.

    Returns:
        Dictionary mapping program names to availability status.
    """
    status: dict[str, bool] = {}

    # Check CLI programs
    for program in REQUIRED_CLI_PROGRAMS + OPTIONAL_CLI_PROGRAMS:
        status[program] = shutil.which(program) is not None

    # Check Python packages
    all_packages = CORE_PYTHON_PACKAGES + FULL_DEPENDENCY_PACKAGES
    for package_name, import_name in all_packages:
        try:
            __import__(import_name)
            status[package_name] = True
        except ImportError:
            status[package_name] = False

    return status


def get_missing_programs() -> tuple[list[str], list[str]]:
    """Get lists of missing required and optional programs.

    Returns:
        Tuple of (missing_required, missing_optional) program lists.
    """
    missing_required: list[str] = []
    missing_optional: list[str] = []

    for program in REQUIRED_CLI_PROGRAMS:
        if shutil.which(program) is None:
            missing_required.append(program)

    for program in OPTIONAL_CLI_PROGRAMS:
        if shutil.which(program) is None:
            missing_optional.append(program)

    # Core packages are required
    for package_name, import_name in CORE_PYTHON_PACKAGES:
        try:
            __import__(import_name)
        except ImportError:
            missing_required.append(package_name)

    # Full dependency packages are optional for basic functionality
    for package_name, import_name in FULL_DEPENDENCY_PACKAGES:
        try:
            __import__(import_name)
        except ImportError:
            missing_optional.append(package_name)

    return missing_required, missing_optional
