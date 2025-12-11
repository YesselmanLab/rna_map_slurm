"""Parameter loading and management."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from rna_map_slurm.config.paths import get_lib_path
from rna_map_slurm.utils.logging import get_logger

log = get_logger("config.parameters")


def get_default_parameters() -> dict[str, Any]:
    """Get default parameters from the bundled YAML file.

    Returns:
        Dictionary of default parameters.
    """
    path = Path(get_lib_path()) / "rna_map_slurm" / "resources" / "default.yml"
    return _load_yaml(path)


def get_parameters_from_file(file_path: str | Path) -> dict[str, Any]:
    """Load parameters from a YAML file, filling in defaults.

    Args:
        file_path: Path to the YAML configuration file.

    Returns:
        Dictionary of parameters with defaults filled in.
    """
    config_data = _load_yaml(Path(file_path))
    return fill_in_missing_default_params(config_data)


def fill_in_missing_default_params(config_data: dict[str, Any]) -> dict[str, Any]:
    """Fill in missing parameters with defaults.

    Args:
        config_data: User-provided configuration.

    Returns:
        Configuration with defaults filled in.
    """
    default_data = get_default_parameters()
    return _fill_in_missing_values(default_data, config_data)


def _load_yaml(path: Path) -> dict[str, Any]:
    """Load a YAML file.

    Args:
        path: Path to YAML file.

    Returns:
        Parsed YAML contents.
    """
    with open(path, encoding="utf-8") as f:
        data: dict[str, Any] = yaml.safe_load(f)
    return data


def _fill_in_missing_values(
    default: dict[str, Any],
    current: dict[str, Any],
) -> dict[str, Any]:
    """Recursively fill in missing values from defaults.

    Args:
        default: Default values dictionary.
        current: Current values dictionary.

    Returns:
        Updated current dictionary.
    """
    for key, value in default.items():
        if isinstance(value, dict):
            node = current.setdefault(key, {})
            _fill_in_missing_values(value, node)
        elif key not in current:
            current[key] = value
    return current
