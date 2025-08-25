#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "c8_tracer/transcribed/CashKarpIntegrator.hpp"

namespace py = pybind11;
using namespace c8_tracer;

void bind_cash_karp(py::module_ &m)
{
  py::module_ integrators = m.def_submodule("integrators", "Numerical integators");

  py::class_<c8_tracer::CashKarpIntegrator>(integrators, "CashKarpIntegrator")
      .def(py::init<double, double, double>())
      .def("AdaptiveStep", [](CashKarpIntegrator &self, const Vec3 &startPos, const Vec3 &startDir, const EnvironmentBase &env, double h0, bool updateStep)
           {
        Vec3 endPos, endDir;
        double h0_copy = h0;
        self.AdaptiveStep(startPos, startDir, endPos, endDir, h0_copy, env, updateStep);
        return py::make_tuple(endPos, endDir, h0_copy); }, py::arg("startPos"), py::arg("startDir"), py::arg("env"), py::arg("h0"), py::arg("updateStep") = true)

      .def("Step", [](CashKarpIntegrator &self, const Vec3 &startPos, const Vec3 &startDir, double h0, const EnvironmentBase &env)
           {
        Vec3 endPos, endDir, dirError;
        self.Step(startPos, startDir, endPos, endDir, dirError, h0, env);
        return py::make_tuple(endPos, endDir, dirError); });
}
