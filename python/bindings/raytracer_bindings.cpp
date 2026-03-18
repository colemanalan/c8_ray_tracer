#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "c8_tracer/transcribed/RayTracer.hpp"
#include "c8_tracer/signal_path.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"

namespace py = pybind11;
using namespace c8_tracer;

void bind_raytracer(py::module_ &m)
{

    py::module_ tracer = m.def_submodule("tracer", "Ray tracers and related materials");

    py::enum_<SolutionMethod>(tracer, "SolutionMethod")
        .value("Brent", SolutionMethod::Brent)
        .value("NGD", SolutionMethod::NGD)
        .export_values();

    // RayTracer2D
    py::class_<RayTracer2D>(tracer, "RayTracer2D")
        .def(py::init<DirectionVector, LengthType, LengthType, double, int>(),
             py::arg("direction"),
             py::arg("minStep") = 0.0001,
             py::arg("maxStep") = 10.0,
             py::arg("tolerance") = 1e-8,
             py::arg("nRays") = 13)
        .def("AddReflectionLayer", &RayTracer2D::AddReflectionLayer)
        .def("GetSignalPathsBrent", &RayTracer2D::GetSignalPathsBrent)
        .def("GetSignalPathsNGD", &RayTracer2D::GetSignalPathsNGD)
        .def("FindEmitAndReceiveBrent", &RayTracer2D::FindEmitAndReceiveBrent)
        .def("FindEmitAndReceiveNGD", &RayTracer2D::FindEmitAndReceiveNGD)
        .def("Get2DProjection", &RayTracer2D::Get2DProjection)
        .def("Get2DRadialDistance", &RayTracer2D::Get2DRadialDistance)
        .def("PrintProfiling", &RayTracer2D::PrintProfiling)
        .def("ResetProfiling", &RayTracer2D::ResetProfiling)
        .def("ShootOneRayToMaximumR", &RayTracer2D::ShootOneRayToMaximumR)
        .def("ShootOneRayToMinimumZ", &RayTracer2D::ShootOneRayToMinimumZ);
    ;
}
