"""Tests for utility functions."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from rna_map_slurm.utils.dataframe import get_data_row, split_dataframe_into_n
from rna_map_slurm.utils.files import get_file_size, random_string


class TestRandomString:
    """Tests for random_string function."""

    def test_correct_length(self) -> None:
        """Test that output has correct length."""
        result = random_string(10)
        assert len(result) == 10

    def test_all_letters(self) -> None:
        """Test that output contains only letters."""
        result = random_string(100)
        assert result.isalpha()

    def test_zero_length(self) -> None:
        """Test zero length string."""
        result = random_string(0)
        assert result == ""


class TestGetFileSize:
    """Tests for get_file_size function."""

    def test_file_size(self, temp_dir: Path) -> None:
        """Test getting file size."""
        test_file = temp_dir / "test.txt"
        test_file.write_text("hello world")

        size = get_file_size(str(test_file))

        assert size == 11  # "hello world" is 11 bytes


class TestGetDataRow:
    """Tests for get_data_row function."""

    def test_single_match(self) -> None:
        """Test finding a single matching row."""
        df = pd.DataFrame({
            "barcode_seq": ["ACGT", "TGCA"],
            "construct": ["RNA1", "RNA2"],
            "value": [1, 2],
        })

        row = get_data_row(df, "ACGT", "RNA1")

        assert row is not None
        assert row["value"] == 1

    def test_no_match(self) -> None:
        """Test when no match is found."""
        df = pd.DataFrame({
            "barcode_seq": ["ACGT"],
            "construct": ["RNA1"],
        })

        row = get_data_row(df, "XXXX", "RNA1")

        assert row is None

    def test_multiple_matches(self) -> None:
        """Test when multiple matches are found."""
        df = pd.DataFrame({
            "barcode_seq": ["ACGT", "ACGT"],
            "construct": ["RNA1", "RNA1"],
        })

        row = get_data_row(df, "ACGT", "RNA1")

        assert row is None


class TestSplitDataframeIntoN:
    """Tests for split_dataframe_into_n function."""

    def test_even_split(self) -> None:
        """Test splitting evenly divisible DataFrame."""
        df = pd.DataFrame({"a": range(10)})

        result = split_dataframe_into_n(df, 5)

        assert len(result) == 5
        assert all(len(chunk) == 2 for chunk in result)

    def test_uneven_split(self) -> None:
        """Test splitting non-evenly divisible DataFrame."""
        df = pd.DataFrame({"a": range(10)})

        result = split_dataframe_into_n(df, 3)

        assert len(result) == 3
        # 10 / 3 = 3 with remainder 1, so sizes should be 4, 3, 3
        sizes = [len(chunk) for chunk in result]
        assert sum(sizes) == 10

    def test_more_chunks_than_rows(self) -> None:
        """Test when n > number of rows."""
        df = pd.DataFrame({"a": range(3)})

        result = split_dataframe_into_n(df, 5)

        assert len(result) == 5
        # Some chunks will be empty
        assert sum(len(chunk) for chunk in result) == 3

    def test_zero_chunks(self) -> None:
        """Test with zero chunks."""
        df = pd.DataFrame({"a": range(5)})

        result = split_dataframe_into_n(df, 0)

        assert len(result) == 1
        assert len(result[0]) == 5
