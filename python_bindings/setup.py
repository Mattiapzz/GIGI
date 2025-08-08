#!/usr/bin/env python3

import pybind11
from pybind11.setup_helpers import Pybind11Extension, build_ext
from pybind11 import get_cmake_dir
import pybind11.setup_helpers as sh
from setuptools import setup, Extension
import glob
import os

# Define the extension module
ext_modules = [
    Pybind11Extension(
        "pygigi",
        [
            "pybind_gigi.cpp",
        ] + glob.glob("../src/*.cc"),  # Include all GIGI source files
        include_dirs=[
            "../src",
            "../third_party/cxxopts/include",
            "../third_party/json/single_include/nlohmann", 
            "../third_party/rapidcsv/src",
        ],
        language='c++',
        cxx_std=17,
    ),
]

setup(
    name="pygigi",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.6",
    description="Python bindings for GIGI optimization library",
    long_description="Python bindings for the GIGI (G-G diagram based) optimization library for vehicle trajectory optimization",
    author="Mattia Piazza",
    author_email="mattia.piazza@unitn.it",
)
