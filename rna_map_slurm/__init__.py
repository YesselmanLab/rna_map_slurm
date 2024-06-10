__author__ = 'Joe Yesselman'
__email__ = 'jyesselm@unl.edu'
__version__ = '0.1.0'

# rna_map_slurm/__init__.py
from celery import Celery

app = Celery('rna_map_slurm')
app.config_from_object('celery_config')

# Ensure tasks are autodiscovered from the tasks module
app.autodiscover_tasks(['rna_map_slurm'])

# Make the app available as an attribute of the module
__all__ = ('app',)