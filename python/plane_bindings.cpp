#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "c8_tracer/plane.hpp"

namespace py = pybind11;
using namespace c8_tracer;

void bind_plane(py::module_ &m)
{

  // Plane
  py::class_<Plane>(m, "Plane")
      .def(py::init<Vec3, Vec3>())
      .def("getNormal", &Plane::getNormal)
      .def("getCenter", &Plane::getCenter);
}
