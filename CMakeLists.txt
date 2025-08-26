cmake_minimum_required(VERSION 3.15)
project(c8_tracer)

find_package(Python COMPONENTS Interpreter Development REQUIRED)
# Include pybind11
find_package(pybind11 REQUIRED)

# Add include directory
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add Python bindings
add_subdirectory(python)
