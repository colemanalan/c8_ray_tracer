#include <pybind11/pybind11.h>
#include "c8_tracer/c8_tracer.hpp"

namespace py = pybind11;

PYBIND11_MODULE(c8_tracer_py, m) {
    m.def("add", &c8_tracer::add, "Add a value");
}