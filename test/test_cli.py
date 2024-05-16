import pandas as pd
import pytest
import os
import logging
from unittest.mock import patch, mock_open


from rna_map_slurm.cli import format_sequencing_run_info

TEST_DIR = os.path.dirname(os.path.realpath(__file__))

# Sample data for testing
df = pd.read_csv(f"{TEST_DIR}/resources/test_fastqs/data.csv")


def test_format_sequencing_run_info(caplog):
    # Create a sample DataFrame for testing
    df = pd.DataFrame(
        {
            "code": ["A", "B", "C"],
            "exp_name": ["Experiment 1", "Experiment 2", "Experiment 3"],
            "construct": ["Construct 1", "Construct 2", "Construct 3"],
            "exp_type": ["Type 1", "Type 2", "Type 3"],
        }
    )

    # Create a sample sequence sheet DataFrame
    df_seq = pd.DataFrame(
        {"code": ["A", "B"], "demultiplex": ["Demultiplex 1", "Demultiplex 2"]}
    )

    # Call the function under test
    format_sequencing_run_info(df)


def test_format_sequencing_run_info_empty_dataframe(caplog):
    df_empty = pd.DataFrame(columns=["code", "exp_name", "construct", "exp_type"])
    with caplog.at_level(logging.ERROR):
        with pytest.raises(SystemExit):
            format_sequencing_run_info(df_empty)
    assert "No sequencing run information found" in caplog.text


# Test case: No corresponding sequence information
def test_format_sequencing_run_info_no_sequence_info(caplog):
    df_no_sequence_info = df.copy()
    df_no_sequence_info["code"] = (
        "XYZ789"  # Assigning codes that do not exist in `sequence_sheet.csv`
    )
    with caplog.at_level(logging.WARNING):
        format_sequencing_run_info(df_no_sequence_info)
    assert "No sequence information found for XYZ789" in caplog.text


# Test case: Mixed sequence information
def test_format_sequencing_run_info_mixed_sequence_info(caplog):
    df_mixed_sequence_info = df.copy()
    df_mixed_sequence_info.loc[0, "code"] = "XYZ789"
    df_mixed_sequence_info.loc[1, "code"] = (
        "ABC123"  # Assigning a code that exists in `sequence_sheet.csv`
    )
    sequence_data = pd.DataFrame(
        {
            "code": ["ABC123", "DEF456"],
            "demultiplex": ["demultiplex_cmd_1", "demultiplex_cmd_2"],
        }
    )
    with patch("gsheets.sheet.get_sequence_sheet", return_value=sequence_data):
        with caplog.at_level(logging.INFO):
            format_sequencing_run_info(df_mixed_sequence_info)
    assert "No sequence information found for XYZ789" in caplog.text
    assert "Found a demultiplexing command for" in caplog.text
