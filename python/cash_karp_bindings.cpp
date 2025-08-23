#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "c8_tracer/transcribed/CashKarpIntegrator.hpp"

namespace py = pybind11;
using namespace c8_tracer;

void bind_cash_karp(py::module_ &m)
{
  py::module_ integrators = m.def_submodule("integrators", "Numerical integators");

  py::class_<c8_tracer::CashKarpIntegrator>(integrators, "CashKarpIntegrator")
      .def(py::init<float, float, float>())
      .def("AdaptiveStep", &c8_tracer::CashKarpIntegrator::AdaptiveStep)
      .def("Step", &c8_tracer::CashKarpIntegrator::Step);
}
