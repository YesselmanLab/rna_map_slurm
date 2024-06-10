import os

from rna_map_slurm.paths import get_lib_path

def write_celery_config_py_file(hostname, path="celery_config.py"):
    template_path = get_lib_path() + "/rna_map_slurm/resources/celery_config.txt"
    with open(template_path, "r") as f:
        template = f.read()
    template = template.replace("{hostname}", hostname)
    with open(path, "w") as f:
        f.write(template)

def run_rabbitmq_in_background():
    pass 

def generate_fastq_splitter_tasks():
    pass