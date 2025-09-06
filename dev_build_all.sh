#!/bin/bash

# rm -r build
cmake -B build
cmake --build build

#uv pip install -e .
