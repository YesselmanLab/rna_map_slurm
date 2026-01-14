"""Compare seqkit and C++ int-demultiplex outputs to verify correctness."""

import gzip
import json
import os
import shutil
import tempfile
from pathlib import Path


def read_fastq_names(path: Path) -> set[str]:
    """Read FASTQ file and return set of read names."""
    names = set()
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:  # Header line
                names.add(line.strip().split()[0])
    return names


def main():
    # Use absolute paths
    script_dir = Path(__file__).parent.parent.resolve()
    test_resources = script_dir / "tests/resources/test_cases/C0098"

    with open(test_resources / "C0098_barcodes.json") as f:
        barcodes = json.load(f)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create demultiplexed directory with test fastq files
        demux_dir = tmpdir / "demultiplexed" / "ACGTACGT"
        demux_dir.mkdir(parents=True)
        shutil.copy(test_resources / "R1.sub.fastq.gz", demux_dir / "test_R1.fastq.gz")
        shutil.copy(test_resources / "R2.sub.fastq.gz", demux_dir / "test_R2.fastq.gz")

        os.chdir(tmpdir)

        # Test first barcode
        bc = barcodes[0]
        b1_seq = bc["barcodes"][0][0]
        b2_seq = bc["barcodes"][0][1]
        bb1 = bc["barcode_bounds"][0][0]
        bb2 = bc["barcode_bounds"][0][1]
        seq_len = len(bc["sequence"])
        b2_min_pos = seq_len - bb2[1]
        b2_max_pos = seq_len - bb2[0]

        print(f"Testing barcode: {bc['full_barcode']}")
        print(f"  5p: {b1_seq} at {bb1}")
        print(f"  3p: {b2_seq} at {bb2} -> mapped to {b2_min_pos}-{b2_max_pos}")

        # Run seqkit
        print("\n--- Running seqkit ---")
        from rna_map_slurm.tasks.int_demultiplex import int_demultiplex

        int_demultiplex(
            construct_barcode="ACGTACGT",
            b1_seq=b1_seq,
            b2_seq=b2_seq,
            b1_min_pos=bb1[0],
            b1_max_pos=bb1[1],
            b2_min_pos=b2_min_pos,
            b2_max_pos=b2_max_pos,
            params={"paths": {"tmp": str(tmpdir / "scratch")}},
        )

        # Run C++
        print("\n--- Running C++ batch ---")
        cpp_dir = tmpdir / "cpp_out" / "ACGTACGT"
        cpp_dir.mkdir(parents=True)

        from rna_map_slurm.tasks.int_demultiplex_cpp import int_demultiplex_batch_cpp

        result = int_demultiplex_batch_cpp(
            r1_path=demux_dir / "test_R1.fastq.gz",
            r2_path=demux_dir / "test_R2.fastq.gz",
            output_dir=cpp_dir,
            barcode_json_path=test_resources / "C0098_barcodes.json",
        )

        # Compare results
        print("\n=== Comparing outputs ===")

        seqkit_out = tmpdir / "int-demultiplexed" / "ACGTACGT"

        # Read seqkit output
        seqkit_mate1 = seqkit_out / f"{b1_seq}_{b2_seq}_mate1.fastq.gz"
        seqkit_mate2 = seqkit_out / f"{b1_seq}_{b2_seq}_mate2.fastq.gz"

        seqkit_names_m1 = read_fastq_names(seqkit_mate1)
        seqkit_names_m2 = read_fastq_names(seqkit_mate2)

        print(f"Seqkit mate1: {len(seqkit_names_m1)} reads")
        print(f"Seqkit mate2: {len(seqkit_names_m2)} reads")

        # Read C++ output
        cpp_mate1 = cpp_dir / f"{bc['full_barcode']}_mate1.fastq.gz"
        cpp_mate2 = cpp_dir / f"{bc['full_barcode']}_mate2.fastq.gz"

        cpp_names_m1 = read_fastq_names(cpp_mate1)
        cpp_names_m2 = read_fastq_names(cpp_mate2)

        print(f"C++ mate1: {len(cpp_names_m1)} reads")
        print(f"C++ mate2: {len(cpp_names_m2)} reads")

        # Compare
        print("\nComparison:")
        print(f"  Seqkit reads: {len(seqkit_names_m1)}")
        print(f"  C++ reads:    {len(cpp_names_m1)}")

        common = seqkit_names_m1 & cpp_names_m1
        only_seqkit = seqkit_names_m1 - cpp_names_m1
        only_cpp = cpp_names_m1 - seqkit_names_m1

        print(f"  Common:       {len(common)}")
        print(f"  Only seqkit:  {len(only_seqkit)}")
        print(f"  Only C++:     {len(only_cpp)}")

        if only_seqkit:
            print(f"  Example seqkit-only: {list(only_seqkit)[:3]}")
        if only_cpp:
            print(f"  Example C++-only: {list(only_cpp)[:3]}")

        # Verify mate1 == mate2 (should be same reads in both)
        print("\nMate consistency check:")
        print(f"  Seqkit mate1 == mate2: {seqkit_names_m1 == seqkit_names_m2}")
        print(f"  C++ mate1 == mate2:    {cpp_names_m1 == cpp_names_m2}")

        if len(common) == len(seqkit_names_m1) == len(cpp_names_m1):
            print("\n✓ PASS: Outputs match exactly!")
        else:
            print("\n✗ FAIL: Outputs differ")
            return 1

        # Test a few more barcodes
        print("\n=== Testing additional barcodes ===")
        all_pass = True

        for bc_idx in [5, 10, 15, 20]:  # Test a few more
            bc = barcodes[bc_idx]
            b1_seq = bc["barcodes"][0][0]
            b2_seq = bc["barcodes"][0][1]
            bb1 = bc["barcode_bounds"][0][0]
            bb2 = bc["barcode_bounds"][0][1]
            seq_len = len(bc["sequence"])
            b2_min_pos = seq_len - bb2[1]
            b2_max_pos = seq_len - bb2[0]

            # Run seqkit for this barcode
            int_demultiplex(
                construct_barcode="ACGTACGT",
                b1_seq=b1_seq,
                b2_seq=b2_seq,
                b1_min_pos=bb1[0],
                b1_max_pos=bb1[1],
                b2_min_pos=b2_min_pos,
                b2_max_pos=b2_max_pos,
                params={"paths": {"tmp": str(tmpdir / "scratch")}},
            )

            # Compare
            seqkit_mate1 = seqkit_out / f"{b1_seq}_{b2_seq}_mate1.fastq.gz"
            cpp_mate1 = cpp_dir / f"{bc['full_barcode']}_mate1.fastq.gz"

            seqkit_names = read_fastq_names(seqkit_mate1)
            cpp_names = read_fastq_names(cpp_mate1)

            match = seqkit_names == cpp_names
            status = "✓" if match else "✗"
            print(f"  {status} {bc['full_barcode']}: seqkit={len(seqkit_names)}, cpp={len(cpp_names)}, match={match}")

            if not match:
                all_pass = False
                only_seqkit = seqkit_names - cpp_names
                only_cpp = cpp_names - seqkit_names
                if only_seqkit:
                    print(f"      Only in seqkit: {list(only_seqkit)[:2]}")
                if only_cpp:
                    print(f"      Only in C++: {list(only_cpp)[:2]}")

        if all_pass:
            print("\n✓ ALL TESTS PASSED!")
        else:
            print("\n✗ SOME TESTS FAILED")
            return 1

    return 0


if __name__ == "__main__":
    exit(main())
