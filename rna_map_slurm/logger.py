import logging
import sys

APP_LOGGER_NAME = "RNA-MAP-SLURM"


def setup_logging(is_debug=False, file_name=None):
    """
    setup logging for everything.
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)  # Set the root logger level

    # Create a stream handler for output to console
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)  # Set the desired level for console output
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    if file_name:
        fh = logging.FileHandler(file_name)
        fh.setFormatter(formatter)
        root_logger.addHandler(fh)


def setup_applevel_logger(logger_name=APP_LOGGER_NAME, is_debug=True, file_name=None):
    """
    setup the logger for the application
    :param logger_name: name of the logger
    :param is_debug: is the logger in debug mode
    :param file_name: name of the file to log to
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG if is_debug else logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    if file_name:
        fh = logging.FileHandler(file_name)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def get_logger(module_name):
    """
    Get the logger for the module
    :param module_name: name of the module
    """
    return logging.getLogger(APP_LOGGER_NAME).getChild(module_name)
