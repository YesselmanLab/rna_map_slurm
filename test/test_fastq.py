"""
testing structured data module
"""

import os
import pytest

from rna_map_slurm.fastq import FastqFile, get_paired_fastqs, pair_fastq_files

TEST_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture
def base_dir():
    return os.path.dirname(os.path.realpath(__file__)) + "/resources/"


@pytest.fixture
def test_dirs(base_dir):
    return {
        "valid": os.path.join(base_dir, "test_fastqs"),
        "compressed": os.path.join(base_dir, "test_fastqs_gziped"),
        "messed": os.path.join(base_dir, "test_fastqs_messed"),
        "non_existent": os.path.join(base_dir, "non_existent"),
    }


def test_fastq():
    path = TEST_DIR + "/resources/test_fastqs/C0098_S1_L001_R1_001.fastq"
    fq = None
    try:
        fq = FastqFile(path)
    except:
        pytest.fail("FastqFile failed to initialize not a valid file")

    assert fq.path == path
    assert fq.is_compressed() is False
    assert fq.is_r1() is True
    assert fq.is_r2() is False


def test_fastq_compressed():
    path = TEST_DIR + "/resources/test_fastqs_gziped/C0098_S1_L001_R1_001.fastq.gz"
    fq = None
    try:
        fq = FastqFile(path)
    except:
        pytest.fail("FastqFile failed to initialize not a valid file")

    assert fq.path == path
    assert fq.is_compressed() is True
    assert fq.is_r1() is True
    assert fq.is_r2() is False


def test_pair_fastq_files():
    f1_paths = [
        "test/resources/test_fastqs/C0098_S1_L001_R1_001.fastq",
        "test/resources/test_fastqs/C0098_S1_L001_R1_002.fastq",
    ]
    f2_paths = [
        "test/resources/test_fastqs/C0098_S1_L001_R2_001.fastq",
        "test/resources/test_fastqs/C0098_S1_L001_R2_002.fastq",
    ]
    paired_files = pair_fastq_files(f1_paths, f2_paths)
    assert len(paired_files) == 2, "Should find 2 paired files"
    for pair in paired_files:
        assert pair.read_1.is_r1(), "First file should be R1"
        assert pair.read_2.is_r2(), "Second file should be R2"
        assert pair.is_compressed() is False, "Should not be compressed"


def test_valid_pairs(test_dirs):
    paired_files = get_paired_fastqs(test_dirs["valid"])
    assert len(paired_files) > 0, "Should find paired fastq files"
    for pair in paired_files:
        assert pair.read_1.is_r1(), "First file should be R1"
        assert pair.read_2.is_r2(), "Second file should be R2"
        assert pair.is_compressed() is False, "Should not be compressed"


def test_compressed_pairs(test_dirs):
    paired_files = get_paired_fastqs(test_dirs["compressed"])
    assert len(paired_files) > 0, "Should find paired fastq files"
    for pair in paired_files:
        assert pair.read_1.is_r1(), "First file should be R1"
        assert pair.read_2.is_r2(), "Second file should be R2"
        assert pair.is_compressed(), "Files should be compressed"


"""
TEST_DIR = os.path.dirname(os.path.realpath(__file__))







def test_fastq_messes():
    path = TEST_DIR + "/resources/test_fastqs_messed/C0098.fastq"
    with pytest.raises(ValueError):
        FastqFile(path)


def test_get_paired_fastqs():

    path = TEST_DIR + "/resources/test_fastqs"
    pfqs = get_paired_fastqs(path)
"""
