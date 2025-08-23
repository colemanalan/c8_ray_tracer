#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "c8_tracer/environment.hpp"
#include "c8_tracer/vec3.hpp"

namespace py = pybind11;
using namespace c8_tracer;

// Trampoline class to allow Python subclasses of environment_base
class PyEnvironmentBase : public EnvironmentBase
{
public:
  using c8_tracer::EnvironmentBase::EnvironmentBase;

  double get_n(const Vec3 &position) const override
  {
    PYBIND11_OVERRIDE_PURE(
        double,                      // Return type
        c8_tracer::EnvironmentBase, // Parent class
        get_n,                      // Name of function
        position                    // Argument
    );
  }

  Vec3 get_grad_n(const Vec3 &position) const override
  {
    PYBIND11_OVERRIDE_PURE(
        Vec3,
        c8_tracer::EnvironmentBase,
        get_grad_n,
        position);
  }
};

void bind_environment(py::module_ &m)
{
  py::module_ environment = m.def_submodule("environment", "Environmental descriptions");

  // add underscore to class name to avoid exposing the base
  py::class_<c8_tracer::EnvironmentBase, PyEnvironmentBase>(environment, "_EnvironmentBase")
      .def(py::init<>())
      .def("get_n", &c8_tracer::EnvironmentBase::get_n)
      .def("get_grad_n", &c8_tracer::EnvironmentBase::get_grad_n);

  py::class_<c8_tracer::IsotropicEnvironment, c8_tracer::EnvironmentBase>(environment, "IsotropicEnvironment")
      .def(py::init<double>())
      .def("get_n", &c8_tracer::IsotropicEnvironment::get_n)
      .def("get_grad_n", &c8_tracer::IsotropicEnvironment::get_grad_n);

  py::class_<c8_tracer::LinearRadialEnvironment, c8_tracer::EnvironmentBase>(environment, "LinearRadialEnvironment")
      .def(py::init<Vec3, double, double>())
      .def("get_n", &c8_tracer::LinearRadialEnvironment::get_n)
      .def("get_grad_n", &c8_tracer::LinearRadialEnvironment::get_grad_n);

  py::class_<c8_tracer::CartesianLinearEnvironment, c8_tracer::EnvironmentBase>(environment, "CartesianLinearEnvironment")
      .def(py::init<Vec3, double, double>())
      .def("get_n", &c8_tracer::CartesianLinearEnvironment::get_n)
      .def("get_grad_n", &c8_tracer::CartesianLinearEnvironment::get_grad_n);
}
