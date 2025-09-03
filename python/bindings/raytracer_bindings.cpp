#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "c8_tracer/transcribed/RayTracer.hpp"
#include "c8_tracer/signal_path.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"

namespace py = pybind11;
using namespace c8_tracer;

// Trampoline class to allow Python subclasses of RayTracerBase
class PyRayTracerBase : public RayTracerBase
{
public:
  using RayTracerBase::RayTracerBase;

  std::vector<SignalPath> PropagateToPoint(
      const Point &start,
      const Point &end,
      const EnvironmentBase &env) override
  {
    PYBIND11_OVERRIDE_PURE(
        std::vector<SignalPath>,
        RayTracerBase,
        PropagateToPoint,
        start, end, env);
  }

  bool FindEmitAndReceive(
      const Point &start,
      const Point &end,
      const EnvironmentBase &env,
      const DirectionVector &seed,
      DirectionVector &emit,
      DirectionVector &receive) override
  {
    PYBIND11_OVERRIDE_PURE(
        bool,
        RayTracerBase,
        FindEmitAndReceive,
        start, end, env, seed, emit, receive);
  }

  SignalPath GetSignalPath(
      const Point &start,
      const DirectionVector &startDir,
      const Point &target,
      const EnvironmentBase &env) override
  {
    PYBIND11_OVERRIDE_PURE(
        SignalPath,
        RayTracerBase,
        GetSignalPath,
        start, startDir, target, env);
  }
};

void bind_raytracer(py::module_ &m)
{

  py::module_ tracer = m.def_submodule("tracer", "Ray tracers and related materials");

  py::class_<RayTracerBase, PyRayTracerBase>(tracer, "RayTracerBase")
      .def(py::init<DirectionVector>())
      .def("PropagateToPoint", &c8_tracer::RayTracerBase::PropagateToPoint)
      .def("FindEmitAndReceive", &c8_tracer::RayTracerBase::FindEmitAndReceive)
      .def("GetSignalPath", &c8_tracer::RayTracerBase::GetSignalPath)
      .def("GetAxis", &RayTracer2D::GetAxis);

  // RayTracer2D
  py::class_<RayTracer2D, RayTracerBase>(tracer, "RayTracer2D")
      .def(py::init<DirectionVector, LengthType, LengthType, double>(),
           py::arg("direction"),
           py::arg("minStep") = 0.0001,
           py::arg("maxStep") = 10.0,
           py::arg("tolerance") = 1e-8)
      .def("AddReflectionLayer", &RayTracer2D::AddReflectionLayer)
      .def("ShootOneRayToMaximumR", &RayTracer2D::ShootOneRayToMaximumR)
      .def("ShootOneRayToMinimumZ", &RayTracer2D::ShootOneRayToMinimumZ)
      .def("FindEmitAndReceiveBrent", &RayTracer2D::FindEmitAndReceiveBrent)
      .def("FindIntersectionWithPlane", &RayTracer2D::FindIntersectionWithPlane)
      .def("ReflectOffPlane", &RayTracer2D::ReflectOffPlane)
      .def("FindRadius", &RayTracer2D::FindRadius)
      .def("Get2DProjection", &RayTracer2D::Get2DProjection)
      .def("Get2DRadialDistance", &RayTracer2D::Get2DRadialDistance)
      .def("PrintProfiling", &RayTracer2D::PrintProfiling)
      .def("ResetProfiling", &RayTracer2D::ResetProfiling);
  ;
}
