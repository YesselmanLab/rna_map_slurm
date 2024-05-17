#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()


with open("README.md", "r", encoding="UTF-8") as f:
    readme = f.read()

with open("requirements.txt", "r", encoding="UTF-8") as f:
    requirements = f.read().splitlines()

setup(
    name="rna_map_slurm",
    version="0.1.0",
    description="a tool that takes care of all processes required to run rna_map on a slurm cluster",
    long_description=readme,
    long_description_content_type="test/markdown",
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    url="https://github.com/jyesselm/rna_map_slurm",
    packages=[
        "rna_map_slurm",
    ],
    package_dir={"rna_map_slurm": "rna_map_slurm"},
    py_modules=[
        "rna_map_slurm/cli",
        "rna_map_slurm/demultiplex",
        "rna_map_slurm/fastq",
        "rna_map_slurm/generate_job",
        "rna_map_slurm/jobs",
        "rna_map_slurm/logger",
        "rna_map_slurm/parameters",
        "rna_map_slurm/paths",
        "rna_map_slurm/run",
        "rna_map_slurm/util",
    ],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords="rna_map_slurm",
    classifiers=[
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: PyPy",
    ],
    entry_points={
        "console_scripts": [
            "rna-map-slurm = rna_map_slurm.cli:cli",
            "rna-map-slurm-runner = rna_map_slurm.run:cli",
        ]
    },
)
