# setup.py
from setuptools import setup
# import sys
# import os
# import pybind11

from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "c8_tracer_py",
        ["python/pywrapper.cpp"],
        include_dirs=["include"],
        cxx_std=17,
    ),
]

setup(
    name="c8_tracer",
    version="0.1.0",
    author="Alan Coleman",
    description="Standalone implementation of C8 Ray Tracer",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
