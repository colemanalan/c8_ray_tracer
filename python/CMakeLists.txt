cmake_minimum_required(VERSION 3.15)
project(c8_tracer LANGUAGES CXX)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

execute_process(
    COMMAND python3 python/bindings_sources.py
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE SOURCE_LIST_RAW
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT SOURCE_LIST_RAW STREQUAL "")
    string(REPLACE " " ";" SOURCE_LIST ${SOURCE_LIST_RAW})
    message(STATUS "Source files: ${SOURCE_LIST}")
else()
    message(FATAL_ERROR "bindings_sources.py returned an empty source list.")
endif()

pybind11_add_module(c8_tracer ${SOURCE_LIST})

# Enable compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(c8_tracer PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Include project headers
target_include_directories(c8_tracer
    PRIVATE
    ${PROJECT_SOURCE_DIR}/../include # Adjust path if needed
)

# Link against Python and other dependencies
target_link_libraries(c8_tracer
    PRIVATE
    Python3::Python
)

# Optional: Set C++ standard
set_target_properties(c8_tracer PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

# Optional: Output to a specific directory
set_target_properties(c8_tracer PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/../build/python
)
