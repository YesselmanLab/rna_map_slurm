import random
import string
import os
import shutil
import gzip
import zipfile

from rna_map_slurm.logger import get_logger

log = get_logger("UTIL")


import random
import string


def get_file_size(file_path):
    file_path = os.path.realpath(file_path)
    return os.path.getsize(file_path)


def random_string(length):
    """
    Generate a random string of specified length.

    Parameters:
        length (int): The length of the random string to generate.

    Returns:
        str: A random string of the specified length.
    """
    return "".join(random.choices(string.ascii_letters, k=length))


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


def flatten_and_zip_directory(input_directory, output_zip):
    with zipfile.ZipFile(output_zip, "w") as zip_ref:
        for root, _, files in os.walk(input_directory):
            for file in files:
                file_path = os.path.join(root, file)
                # Save the file in the zip with only its base name
                zip_ref.write(file_path, os.path.basename(file))


def get_data_row(df, barcode_seq, rna_name):
    df_sub = df.query("barcode_seq == @barcode_seq and construct == @rna_name")
    if len(df_sub) == 0:
        log.error(
            f"barcode_seq {barcode_seq} with rna {rna_name} not found in data.csv"
        )
        return None
    if len(df_sub) > 1:
        log.warning(
            f"barcode_seq {barcode_seq} with rna {rna_name} has multiple entries in data.csv"
        )
        return None
    return df_sub.iloc[0]
