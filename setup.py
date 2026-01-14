"""Setup script for building C++ extension with pybind11."""

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

ext_modules = [
    Pybind11Extension(
        "rna_map_slurm.cpp",
        ["src/fastq_filter.cpp"],
        libraries=["z"],  # Link zlib
        extra_compile_args=["-O3"],
    ),
]

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
