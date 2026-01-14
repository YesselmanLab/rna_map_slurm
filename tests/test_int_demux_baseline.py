"""Baseline test for int_demultiplex to measure seqkit performance."""

import json
import os
import shutil
import tempfile
import time
from pathlib import Path

from seq_tools.sequence import get_reverse_complement


def setup_test_dir(test_resources: Path, tmpdir: Path) -> None:
    """Set up test directory structure expected by int_demultiplex."""
    # Create demultiplexed directory with test fastq files
    demux_dir = tmpdir / "demultiplexed" / "ACGTACGT"
    demux_dir.mkdir(parents=True)

    # Copy test FASTQ files
    shutil.copy(test_resources / "R1.sub.fastq.gz", demux_dir / "test_R1.fastq.gz")
    shutil.copy(test_resources / "R2.sub.fastq.gz", demux_dir / "test_R2.fastq.gz")

    # Create output directory
    (tmpdir / "int-demultiplexed").mkdir(exist_ok=True)


def run_seqkit_baseline(test_resources: Path, tmpdir: Path, barcode: dict) -> float:
    """Run seqkit-based int_demultiplex and return elapsed time."""
    from rna_map_slurm.tasks.int_demultiplex import int_demultiplex

    # Extract barcode info
    b1_seq = barcode["barcodes"][0][0]  # 5' barcode
    b2_seq = barcode["barcodes"][0][1]  # 3' barcode
    bb1 = barcode["barcode_bounds"][0][0]  # [47, 54]
    bb2 = barcode["barcode_bounds"][0][1]  # [139, 146]

    # Map 3' barcode position (from end of sequence)
    seq_len = len(barcode["sequence"])
    b2_min_pos = seq_len - bb2[1]  # 170 - 146 = 24
    b2_max_pos = seq_len - bb2[0]  # 170 - 139 = 31

    construct_barcode = "ACGTACGT"  # Library barcode from data.csv

    start = time.time()
    int_demultiplex(
        construct_barcode=construct_barcode,
        b1_seq=b1_seq,
        b2_seq=b2_seq,
        b1_min_pos=bb1[0],
        b1_max_pos=bb1[1],
        b2_min_pos=b2_min_pos,
        b2_max_pos=b2_max_pos,
        params={"paths": {"tmp": str(tmpdir / "scratch")}},
    )
    elapsed = time.time() - start

    return elapsed


def run_cpp_baseline(test_resources: Path, tmpdir: Path, barcode: dict) -> float:
    """Run C++ fastq_filter (single barcode) and return elapsed time."""
    try:
        from rna_map_slurm.cpp import process_fastq_files
    except ImportError:
        print("C++ module not available, skipping")
        return -1

    # Extract barcode info
    b1_seq = barcode["barcodes"][0][0]  # 5' barcode
    b2_seq = barcode["barcodes"][0][1]  # 3' barcode
    bb1 = barcode["barcode_bounds"][0][0]  # [47, 54]
    bb2 = barcode["barcode_bounds"][0][1]  # [139, 146]

    # Map 3' barcode position (from end of sequence)
    seq_len = len(barcode["sequence"])
    b2_min_pos = seq_len - bb2[1]
    b2_max_pos = seq_len - bb2[0]

    # Get reverse complement for 3' barcode
    b2_seq_rc = get_reverse_complement(b2_seq.replace("U", "T"))
    b1_seq_dna = b1_seq.replace("U", "T")

    r1_path = str(tmpdir / "demultiplexed" / "ACGTACGT" / "test_R1.fastq.gz")
    r2_path = str(tmpdir / "demultiplexed" / "ACGTACGT" / "test_R2.fastq.gz")
    output_dir = str(tmpdir / "int-demultiplexed" / "ACGTACGT")
    os.makedirs(output_dir, exist_ok=True)

    start = time.time()
    # Note: Current C++ processes R1 first, R2 second (opposite of seqkit order)
    # barcode1 searches in R1, barcode2 searches in R2
    process_fastq_files(
        r2_path,  # R2 has the 5' barcode
        r1_path,  # R1 has the 3' barcode (RC)
        output_dir,
        b1_seq_dna,  # 5' barcode
        bb1[0], bb1[1],
        b2_seq_rc,  # 3' barcode reverse complement
        b2_min_pos, b2_max_pos,
    )
    elapsed = time.time() - start

    return elapsed


def run_cpp_batch(test_resources: Path, tmpdir: Path, barcodes: list) -> float:
    """Run C++ batch mode with all barcodes and return elapsed time."""
    try:
        from rna_map_slurm.cpp import process_fastq_files_batch
    except ImportError:
        print("C++ batch module not available, skipping")
        return -1

    # Prepare barcode vectors
    barcodes_5p = []
    barcodes_3p_rc = []
    starts_5p = []
    ends_5p = []
    starts_3p = []
    ends_3p = []
    output_names = []

    for bc in barcodes:
        b1_seq = bc["barcodes"][0][0].replace("U", "T")  # 5' barcode
        b2_seq = bc["barcodes"][0][1].replace("U", "T")  # 3' barcode
        b2_seq_rc = get_reverse_complement(b2_seq)

        bb1 = bc["barcode_bounds"][0][0]  # [47, 54]
        bb2 = bc["barcode_bounds"][0][1]  # [139, 146]

        # Map 3' barcode position (from end of sequence)
        seq_len = len(bc["sequence"])
        b2_min_pos = seq_len - bb2[1]
        b2_max_pos = seq_len - bb2[0]

        barcodes_5p.append(b1_seq)
        barcodes_3p_rc.append(b2_seq_rc)
        starts_5p.append(bb1[0])
        ends_5p.append(bb1[1])
        starts_3p.append(b2_min_pos)
        ends_3p.append(b2_max_pos)
        output_names.append(bc["full_barcode"])

    r1_path = str(tmpdir / "demultiplexed" / "ACGTACGT" / "test_R1.fastq.gz")
    r2_path = str(tmpdir / "demultiplexed" / "ACGTACGT" / "test_R2.fastq.gz")
    output_dir = str(tmpdir / "int-demultiplexed" / "ACGTACGT")
    os.makedirs(output_dir, exist_ok=True)

    start = time.time()
    match_counts = process_fastq_files_batch(
        r1_path,
        r2_path,
        output_dir,
        barcodes_5p,
        barcodes_3p_rc,
        starts_5p,
        ends_5p,
        starts_3p,
        ends_3p,
        output_names,
        max_open_files=200,
    )
    elapsed = time.time() - start

    print(f"Match counts: {dict(zip(output_names, match_counts))}")
    return elapsed


def main():
    # Use absolute path so it works after chdir
    script_dir = Path(__file__).parent.parent.resolve()
    test_resources = script_dir / "tests/resources/test_cases/C0098"

    with open(test_resources / "C0098_barcodes.json") as f:
        barcodes = json.load(f)

    print(f"Loaded {len(barcodes)} barcodes")

    # Test with first barcode
    barcode = barcodes[0]
    print(f"\nTesting with barcode: {barcode['full_barcode']}")
    print(f"  5' barcode: {barcode['barcodes'][0][0]} at {barcode['barcode_bounds'][0][0]}")
    print(f"  3' barcode: {barcode['barcodes'][0][1]} at {barcode['barcode_bounds'][0][1]}")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        os.chdir(tmpdir)

        setup_test_dir(test_resources, tmpdir)

        # Run seqkit baseline
        print("\n--- Running seqkit baseline ---")
        try:
            seqkit_time = run_seqkit_baseline(test_resources, tmpdir, barcode)
            print(f"Seqkit time: {seqkit_time:.3f}s")

            # Check output
            output_dir = tmpdir / "int-demultiplexed" / "ACGTACGT"
            outputs = list(output_dir.glob("*.fastq.gz"))
            print(f"Output files: {[f.name for f in outputs]}")
            for f in outputs:
                print(f"  {f.name}: {f.stat().st_size} bytes")
        except Exception as e:
            print(f"Seqkit failed: {e}")
            seqkit_time = -1

        # Run C++ baseline
        print("\n--- Running C++ baseline ---")
        try:
            cpp_time = run_cpp_baseline(test_resources, tmpdir, barcode)
            if cpp_time > 0:
                print(f"C++ time: {cpp_time:.3f}s")

                output_dir = tmpdir / "int-demultiplexed" / "ACGTACGT"
                outputs = list(output_dir.glob("*.fastq.gz"))
                print(f"Output files: {[f.name for f in outputs]}")
                for f in outputs:
                    print(f"  {f.name}: {f.stat().st_size} bytes")
        except Exception as e:
            print(f"C++ failed: {e}")
            cpp_time = -1

        if seqkit_time > 0 and cpp_time > 0:
            print(f"\n--- Comparison (single barcode) ---")
            print(f"Seqkit: {seqkit_time:.3f}s")
            print(f"C++:    {cpp_time:.3f}s")
            print(f"Speedup: {seqkit_time / cpp_time:.1f}x")

        # Run C++ batch mode with all barcodes
        print("\n--- Running C++ batch mode (all 24 barcodes) ---")
        try:
            # Need fresh output directory for batch
            batch_output = tmpdir / "int-demultiplexed-batch" / "ACGTACGT"
            batch_output.mkdir(parents=True, exist_ok=True)

            batch_time = run_cpp_batch(test_resources, tmpdir, barcodes)
            if batch_time > 0:
                print(f"C++ batch time: {batch_time:.3f}s for {len(barcodes)} barcodes")

                output_dir = tmpdir / "int-demultiplexed" / "ACGTACGT"
                outputs = list(output_dir.glob("*.fastq.gz"))
                print(f"Output files: {len(outputs)} files")
                non_empty = [f for f in outputs if f.stat().st_size > 20]
                print(f"Non-empty output files: {len(non_empty)}")

                # Compare to seqkit running all barcodes
                seqkit_all_time = seqkit_time * len(barcodes)
                print(f"\n--- Comparison (all {len(barcodes)} barcodes) ---")
                print(f"Seqkit estimated: {seqkit_all_time:.3f}s ({len(barcodes)} × {seqkit_time:.3f}s)")
                print(f"C++ batch:        {batch_time:.3f}s")
                print(f"Speedup:          {seqkit_all_time / batch_time:.1f}x")
        except Exception as e:
            print(f"C++ batch failed: {e}")
            import traceback
            traceback.print_exc()


if __name__ == "__main__":
    main()
