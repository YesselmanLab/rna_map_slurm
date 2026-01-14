"""C++ accelerated internal demultiplexing."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from seq_tools.sequence import get_reverse_complement

from rna_map_slurm.utils.logging import get_logger

log = get_logger("tasks.int_demultiplex_cpp")


def int_demultiplex_batch_cpp(
    r1_path: str | Path,
    r2_path: str | Path,
    output_dir: str | Path,
    barcode_json_path: str | Path,
    max_open_files: int = 200,
) -> dict[str, int]:
    """Perform batch internal demultiplexing using C++ implementation.

    This is ~50x faster than the seqkit-based implementation by:
    1. Processing all barcodes in a single pass through the input files
    2. Avoiding intermediate files and subprocess calls

    Args:
        r1_path: Path to R1 FASTQ file (contains 3' barcode, reverse complemented).
        r2_path: Path to R2 FASTQ file (contains 5' barcode, forward).
        output_dir: Directory for output files.
        barcode_json_path: Path to barcode JSON file with barcode specifications.
        max_open_files: Maximum number of file handles to open at once.

    Returns:
        Dictionary mapping barcode names to match counts.

    Raises:
        ImportError: If C++ module is not available.
        FileNotFoundError: If input files or barcode JSON not found.
    """
    from rna_map_slurm.cpp import process_fastq_files_batch

    # Load barcode JSON
    with open(barcode_json_path) as f:
        barcodes = json.load(f)

    if not barcodes:
        log.warning("No barcodes found in JSON file")
        return {}

    # Prepare barcode vectors for C++
    barcodes_5p: list[str] = []
    barcodes_3p_rc: list[str] = []
    starts_5p: list[int] = []
    ends_5p: list[int] = []
    starts_3p: list[int] = []
    ends_3p: list[int] = []
    output_names: list[str] = []

    for bc in barcodes:
        barcode_info = _extract_barcode_info(bc)
        if barcode_info is None:
            continue

        barcodes_5p.append(barcode_info["barcode_5p"])
        barcodes_3p_rc.append(barcode_info["barcode_3p_rc"])
        starts_5p.append(barcode_info["start_5p"])
        ends_5p.append(barcode_info["end_5p"])
        starts_3p.append(barcode_info["start_3p"])
        ends_3p.append(barcode_info["end_3p"])
        output_names.append(barcode_info["output_name"])

    if not barcodes_5p:
        log.warning("No valid barcodes extracted from JSON")
        return {}

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    log.info(f"Processing {len(barcodes_5p)} barcodes with C++ batch mode")

    # Call C++ batch function
    match_counts = process_fastq_files_batch(
        str(r1_path),
        str(r2_path),
        str(output_dir),
        barcodes_5p,
        barcodes_3p_rc,
        starts_5p,
        ends_5p,
        starts_3p,
        ends_3p,
        output_names,
        max_open_files,
    )

    # Build result dictionary
    result = dict(zip(output_names, match_counts))

    total_matches = sum(match_counts)
    log.info(f"Batch demultiplexing complete: {total_matches} total matches")

    return result


def _extract_barcode_info(bc: dict[str, Any], position_buffer: int = 2) -> dict[str, Any] | None:
    """Extract barcode information from a barcode record.

    Args:
        bc: Barcode record from JSON.
        position_buffer: Buffer to add to position ranges (default 2, matches seqkit behavior).

    Returns:
        Dictionary with extracted barcode info, or None if invalid.
    """
    try:
        # Get barcode sequences
        b1_seq = bc["barcodes"][0][0].replace("U", "T")  # 5' barcode
        b2_seq = bc["barcodes"][0][1].replace("U", "T")  # 3' barcode
        b2_seq_rc = get_reverse_complement(b2_seq)

        # Get position bounds
        bb1 = bc["barcode_bounds"][0][0]  # [start, end] for 5' barcode
        bb2 = bc["barcode_bounds"][0][1]  # [start, end] for 3' barcode

        # Map 3' barcode position (from end of sequence to read position)
        seq_len = len(bc["sequence"])
        b2_min_pos = seq_len - bb2[1]
        b2_max_pos = seq_len - bb2[0]

        # Add position buffer to match seqkit behavior (±2 bp)
        return {
            "barcode_5p": b1_seq,
            "barcode_3p_rc": b2_seq_rc,
            "start_5p": bb1[0] - position_buffer,
            "end_5p": bb1[1] + position_buffer,
            "start_3p": b2_min_pos - position_buffer,
            "end_3p": b2_max_pos + position_buffer,
            "output_name": bc["full_barcode"],
        }
    except (KeyError, IndexError, TypeError) as e:
        log.warning(f"Failed to extract barcode info: {e}")
        return None


def int_demultiplex_single_cpp(
    r1_path: str | Path,
    r2_path: str | Path,
    output_dir: str | Path,
    barcode_5p: str,
    barcode_3p: str,
    start_5p: int,
    end_5p: int,
    start_3p: int,
    end_3p: int,
) -> None:
    """Perform single barcode internal demultiplexing using C++.

    This is a wrapper around the legacy single-barcode C++ function.

    Args:
        r1_path: Path to R1 FASTQ file.
        r2_path: Path to R2 FASTQ file.
        output_dir: Directory for output files.
        barcode_5p: 5' barcode sequence (DNA).
        barcode_3p: 3' barcode sequence (DNA, will be reverse complemented).
        start_5p: Start position for 5' barcode search.
        end_5p: End position for 5' barcode search.
        start_3p: Start position for 3' barcode search.
        end_3p: End position for 3' barcode search.
    """
    from rna_map_slurm.cpp import process_fastq_files

    b3_seq_rc = get_reverse_complement(barcode_3p)

    os.makedirs(output_dir, exist_ok=True)

    # The legacy function searches barcode1 in file1, barcode2 in file2
    # R2 has 5' barcode, R1 has 3' barcode (RC)
    process_fastq_files(
        str(r2_path),  # R2 has 5' barcode
        str(r1_path),  # R1 has 3' barcode (RC)
        str(output_dir),
        barcode_5p,
        start_5p,
        end_5p,
        b3_seq_rc,
        start_3p,
        end_3p,
    )
