import os
import shutil
import pandas as pd

from rna_map_slurm.fastq import get_paired_fastqs
from rna_map_slurm.demultiplex import SabreDemultiplexer

TEST_DIR = os.path.dirname(os.path.realpath(__file__))


def test_sabre_demultiplex():
    """
    test demultiplex
    """
    # TODO make sure to do more file checks
    path = f"{TEST_DIR}/test_run/"
    os.makedirs(path, exist_ok=True)
    df = pd.read_csv(f"{TEST_DIR}/resources/test_fastqs/data.csv")
    pfqs = get_paired_fastqs(f"{TEST_DIR}/resources/test_fastqs/")
    demultiplexer = SabreDemultiplexer()
    demultiplexer.run(df, pfqs[0], path)
    assert os.path.exists(f"{TEST_DIR}/test_run/NC/test_R1.fastq")
    assert os.path.exists(f"{TEST_DIR}/test_run/NC/test_R2.fastq")
    dirs = [
        "ACAAAATGGTGG",
    ]
    # shutil.rmtree(f"{TEST_DIR}/test_run")


def _test_sabre_demultiplex_gziped():
    """
    test demultiplex
    """
    # TODO make sure to do more file checks
    path = f"{TEST_DIR}/test_run/"
    os.makedirs(path, exist_ok=True)
    df = pd.read_csv(f"{TEST_DIR}/resources/test_fastqs_gziped/data.csv")
    pfqs = get_paired_fastqs(f"{TEST_DIR}/resources/test_fastqs_gziped/")
    demultiplexer = SabreDemultiplexer()
    demultiplexer.run(df, pfqs[0], path)
    shutil.rmtree(f"{TEST_DIR}/test_run")
