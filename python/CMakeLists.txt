# only needed for building by hand

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

pybind11_add_module(c8_tracer_py pywrapper.cpp)

target_include_directories(c8_tracer_py PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(c8_tracer_py PRIVATE Python3::Python)