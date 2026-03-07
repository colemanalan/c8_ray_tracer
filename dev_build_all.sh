#!/bin/bash

# rm -r build
cmake -B build
cmake --build build -j4

#uv pip install -e .
