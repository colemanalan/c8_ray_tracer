cmake_minimum_required(VERSION 3.14)
project(c8_tracer LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add include directory
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add Python wrapper
add_subdirectory(python)
