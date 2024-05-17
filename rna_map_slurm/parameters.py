import yaml

from rna_map_slurm.paths import get_lib_path
from rna_map_slurm.logger import get_logger

log = get_logger(__name__)


def fill_in_missing_default_dict_values(default, current):
    """
    Recursively fill in missing values in the current dictionary with values
    from the default dictionary.

    Parameters:
    - default (dict): The default dictionary containing the values to fill in.
    - current (dict): The current dictionary to be updated with missing values.

    Returns:
    - dict: The updated current dictionary with missing values filled in.
    """
    for key, value in default.items():
        if isinstance(value, dict):
            # If the value is a dictionary, recurse into it
            node = current.setdefault(key, {})
            fill_in_missing_default_dict_values(value, node)
        elif key not in current:
            # Set value if key is missing
            current[key] = value
    return current


def fill_in_missing_default_params(config_data):
    """
    Fills in missing default parameters in the given config_data dictionary.

    Parameters:
    - config_data (dict): The dictionary containing the configuration data.

    Returns:
    - dict: The updated config_data dictionary with missing default parameters filled in.
    """
    path = get_lib_path() + "/rna_map_slurm/resources/default.yml"
    with open(path, "r") as yaml_file:
        default_data = yaml.safe_load(yaml_file)
    fill_in_missing_default_dict_values(default_data, config_data)
    return config_data


def get_default_parameters():
    """
    Get the default parameters from the default YAML file.

    Returns:
    - dict: The dictionary containing the default parameters.
    """
    path = get_lib_path() + "/rna_map_slurm/resources/default.yml"
    with open(path, "r") as yaml_file:
        default_data = yaml.safe_load(yaml_file)
    return default_data


def get_parameters_from_file(file_path):
    """
    Get parameters from a YAML file and fills in missing default values

    Parameters:
    - file_path (str): The path to the YAML file containing the parameters.

    Returns:
    - dict: The dictionary containing the parameters from the YAML file.
    """
    config_data = {}
    with open(file_path, "r") as yaml_file:
        config_data = yaml.safe_load(yaml_file)
    return fill_in_missing_default_params(config_data)
