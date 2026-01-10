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
    use_rich: bool = True,
) -> None:
    """Configure logging for the application.

    Sets up console output with optional rich formatting and file logging.

    Args:
        is_debug: If True, set log level to DEBUG; otherwise INFO.
        file_name: Optional path to log file for persistent logging.
        use_rich: If True, use rich for colored console output.
    """
    root_logger = logging.getLogger()
    log_level = logging.DEBUG if is_debug else logging.INFO
    root_logger.setLevel(log_level)

    # Clear existing handlers to avoid duplicates
    root_logger.handlers.clear()

    if use_rich:
        try:
            from rich.logging import RichHandler

            console_handler = RichHandler(
                level=logging.INFO,
                show_time=False,
                show_path=False,
                markup=True,
                rich_tracebacks=True,
            )
            # Use simpler format since RichHandler adds its own formatting
            console_handler.setFormatter(logging.Formatter("%(message)s"))
        except ImportError:
            # Fall back to standard handler if rich not available
            console_handler = _create_standard_handler()
    else:
        console_handler = _create_standard_handler()

    root_logger.addHandler(console_handler)

    if file_name is not None:
        # File handler always uses plain text format
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        file_handler = logging.FileHandler(str(file_name))
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)


def _create_standard_handler() -> logging.StreamHandler:
    """Create a standard (non-rich) console handler."""
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    handler.setFormatter(formatter)
    return handler


def get_logger(module_name: str) -> logging.Logger:
    """Get a logger for a specific module.

    Args:
        module_name: Name of the module requesting the logger.

    Returns:
        Logger instance configured as a child of the app logger.
    """
    return logging.getLogger(APP_LOGGER_NAME).getChild(module_name)
