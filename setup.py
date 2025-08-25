from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

# Add python/ to sys.path so we can import bindings_sources
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "python", "bindings"))
from bindings_sources import get_sources

ext_modules = [
    Pybind11Extension(
        "c8_tracer",
        get_sources(),
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
