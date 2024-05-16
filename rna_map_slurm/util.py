import random
import string
import os
import shutil
import gzip

from rna_map_slurm.logger import get_logger

log = get_logger(__name__)


import random
import string


def random_string(length):
    """
    Generate a random string of specified length.

    Parameters:
        length (int): The length of the random string to generate.

    Returns:
        str: A random string of the specified length.
    """
    return "".join(random.choices(string.ascii_letters, k=length))


import os
import gzip
import shutil


def gzip_files(directory):
    """
    Compresses all files in the specified directory (excluding already compressed files) using gzip compression.

    Args:
        directory (str): The directory path where the files are located.

    Returns:
        None
    """
    for root, dirs, files in os.walk(directory):
        for file in files:
            if not file.endswith(".gz"):  # Ignore already compressed files
                file_path = os.path.join(root, file)
                compressed_file_path = f"{file_path}.gz"
                with open(file_path, "rb") as f_in:
                    with gzip.open(compressed_file_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)

                os.remove(file_path)  # Remove the original file
