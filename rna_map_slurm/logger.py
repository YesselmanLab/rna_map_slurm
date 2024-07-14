import logging
import sys

APP_LOGGER_NAME = "rna-map-slurm"


def setup_logging(is_debug=False, file_name=None):
    """
    setup logging for everything.
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)  # Set the root logger level

    # Create a stream handler for output to console
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)  # Set the desired level for console output
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    if file_name:
        fh = logging.FileHandler(file_name)
        fh.setFormatter(formatter)
        root_logger.addHandler(fh)


def get_logger(module_name):
    """
    Get the logger for the module
    :param module_name: name of the module
    """
    return logging.getLogger(APP_LOGGER_NAME).getChild(module_name)
