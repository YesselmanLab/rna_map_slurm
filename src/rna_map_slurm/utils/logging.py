"""Logging configuration and utilities."""

from __future__ import annotations

import logging
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

APP_LOGGER_NAME = "rna-map-slurm"


def setup_logging(
    is_debug: bool = False,
    file_name: str | Path | None = None,
) -> None:
    """Configure logging for the application.

    Sets up console output and optionally file logging.

    Args:
        is_debug: If True, set log level to DEBUG; otherwise INFO.
        file_name: Optional path to log file for persistent logging.
    """
    root_logger = logging.getLogger()
    log_level = logging.DEBUG if is_debug else logging.INFO
    root_logger.setLevel(log_level)

    # Clear existing handlers to avoid duplicates
    root_logger.handlers.clear()

    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    if file_name is not None:
        file_handler = logging.FileHandler(str(file_name))
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)


def get_logger(module_name: str) -> logging.Logger:
    """Get a logger for a specific module.

    Args:
        module_name: Name of the module requesting the logger.

    Returns:
        Logger instance configured as a child of the app logger.
    """
    return logging.getLogger(APP_LOGGER_NAME).getChild(module_name)
