#include <pybind11/pybind11.h>
#include "c8_tracer/c8_tracer.hpp"

namespace py = pybind11;
using namespace c8_tracer;

void bind_tracer(py::module_ &m)
{
  m.def("add", &c8_tracer::add, "Add a value");
}
