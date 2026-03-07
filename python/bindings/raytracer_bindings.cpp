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

    std::vector<SignalPath> GetSignalPaths(
        const Point &start,
        const Point &end,
        const EnvironmentBase &env) override
    {
        PYBIND11_OVERRIDE_PURE(
            std::vector<SignalPath>,
            RayTracerBase,
            GetSignalPaths,
            start, end, env);
    }
};

void bind_raytracer(py::module_ &m)
{

    py::module_ tracer = m.def_submodule("tracer", "Ray tracers and related materials");

    py::enum_<SolutionMethod>(tracer, "SolutionMethod")
        .value("Brent", SolutionMethod::Brent)
        .value("NGD", SolutionMethod::NGD)
        .export_values();

    py::class_<RayTracerBase, PyRayTracerBase>(tracer, "RayTracerBase")
        .def(py::init<DirectionVector>())
        .def("PropagateToPoint", &c8_tracer::RayTracerBase::PropagateToPoint)
        .def("FindEmitAndReceive", &c8_tracer::RayTracerBase::FindEmitAndReceive)
        .def("GetSignalPath", &c8_tracer::RayTracerBase::GetSignalPath)
        .def("GetSignalPaths", &c8_tracer::RayTracerBase::GetSignalPaths)
        .def("GetAxis", &RayTracer2D::GetAxis);

    // RayTracer2D
    py::class_<RayTracer2D, RayTracerBase>(tracer, "RayTracer2D")
        .def(py::init<DirectionVector, LengthType, LengthType, double, int>(),
             py::arg("direction"),
             py::arg("minStep") = 0.0001,
             py::arg("maxStep") = 10.0,
             py::arg("tolerance") = 1e-8,
             py::arg("nRays") = 21)
        .def("AddReflectionLayer", &RayTracer2D::AddReflectionLayer)
        .def("FindEmitAndReceiveBrent", &RayTracer2D::FindEmitAndReceiveBrent)
        .def("Get2DProjection", &RayTracer2D::Get2DProjection)
        .def("Get2DRadialDistance", &RayTracer2D::Get2DRadialDistance)
        .def("GetSignalPaths", &RayTracer2D::GetSignalPaths)
        .def("PrintProfiling", &RayTracer2D::PrintProfiling)
        .def("ResetProfiling", &RayTracer2D::ResetProfiling)
        .def("ShootOneRayToMaximumR", &RayTracer2D::ShootOneRayToMaximumR)
        .def("ShootOneRayToMinimumZ",
             static_cast<bool (RayTracer2D::*)(
                 Point const &, DirectionVector const &, Point &,
                 DirectionVector &, Point const &,
                 EnvironmentBase const &)>(&RayTracer2D::ShootOneRayToMinimumZ))
        .def("ShootOneRayToMinimumZ",
             static_cast<bool (RayTracer2D::*)(
                 Point const &, DirectionVector const &, Point &,
                 DirectionVector &, Plane const &,
                 EnvironmentBase const &)>(&RayTracer2D::ShootOneRayToMinimumZ));
}
