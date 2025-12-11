"""Timing utilities and decorators."""

from __future__ import annotations

import datetime
from collections.abc import Callable
from functools import wraps
from typing import Any, TypeVar

from rna_map_slurm.utils.logging import get_logger

log = get_logger("timing")

F = TypeVar("F", bound=Callable[..., Any])


def time_it(func: F) -> F:
    """Decorator to measure and log function execution time.

    Args:
        func: Function to wrap with timing.

    Returns:
        Wrapped function that logs execution time.

    Example:
        @time_it
        def slow_function():
            ...
    """

    @wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        start_time = datetime.datetime.now()
        result = func(*args, **kwargs)
        end_time = datetime.datetime.now()
        elapsed = end_time - start_time
        log.info(f"Function '{func.__name__}' executed in {elapsed.total_seconds():.4f} seconds")
        return result

    return wrapper  # type: ignore[return-value]
