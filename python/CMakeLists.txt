cmake_minimum_required(VERSION 3.15)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

execute_process(
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/python/bindings/bindings_sources.py
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

pybind11_add_module(c8_tracer_ext ${SOURCE_LIST})

# Enable compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(c8_tracer_ext PRIVATE -Wall -Wextra -Wpedantic)
endif()

target_include_directories(c8_tracer_ext PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Link against Python and other dependencies
target_link_libraries(c8_tracer_ext
    PRIVATE
    Python3::Python
)

set_target_properties(c8_tracer_ext PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/c8_tracer
)
install(TARGETS c8_tracer_ext
    LIBRARY DESTINATION c8_tracer)
