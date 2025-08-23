#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "c8_tracer/transcribed/RayTracer.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"

namespace py = pybind11;
using namespace c8_tracer;

void bind_raytracer(py::module_ &m)
{

  // RayTracer2D
  py::class_<RayTracer2D>(m, "RayTracer2D")
      .def(py::init<DirectionVector, LengthType, LengthType, double>())
      .def("AddReflectionLayer", &RayTracer2D::AddReflectionLayer)
      .def("PropagateToPoint", &RayTracer2D::PropagateToPoint)
      // .def("FindEmitAndReceive", &RayTracer2D::FindEmitAndReceive)
      // .def("ShootOneRayToMaximumR", &RayTracer2D::ShootOneRayToMaximumR)
      // .def("ShootOneRayToMinimumZ", &RayTracer2D::ShootOneRayToMinimumZ)
      // .def("FindIntersectionWithPlane", &RayTracer2D::FindIntersectionWithPlane)
      // .def("ReflectOffPlane", &RayTracer2D::ReflectOffPlane)
      //     .def("FindRadius", &RayTracer2D::FindRadius)
      //     .def("Get2DProjection", &RayTracer2D::Get2DProjection)
      //     .def("Get2DRadialDistance", &RayTracer2D::Get2DRadialDistance)
      //     .def("GetSignalPath", &RayTracer2D::GetSignalPath)
      //     .def("GetAxis", &RayTracer2D::GetAxis)
      //     .def("PrintProfiling", &RayTracer2D::PrintProfiling)
      //     .def("ResetProfiling", &RayTracer2D::ResetProfiling);
      ;
}
