#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "c8_tracer/transcribed/Path.hpp"
#include "c8_tracer/signal_path.hpp"
#include "c8_tracer/vec3.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"

namespace py = pybind11;
using namespace c8_tracer;

void bind_path(py::module_ &m)
{

    py::module_ path = m.def_submodule("path", "Radio propagation paths");

    // Path
    py::class_<Path>(path, "Path")
        .def(py::init<Point>())
        .def(py::init<std::deque<Point>>())
        .def("addToEnd", &Path::addToEnd)
        .def("removeFromEnd", &Path::removeFromEnd)
        .def("getLength", &Path::getLength)
        .def("getStart", &Path::getStart)
        .def("getEnd", &Path::getEnd)
        .def("getPoint", &Path::getPoint)
        .def("getNSegments", &Path::getNSegments);

    // SignalPath
    py::class_<SignalPath>(path, "SignalPath")
        .def(py::init<double, double, double, double, Vec3, Vec3, double, std::deque<Point>, double, double>())
        .def_readwrite("propagation_time", &SignalPath::propagation_time_)
        .def_readwrite("average_refractive_index", &SignalPath::average_refractive_index_)
        .def_readwrite("refractive_index_source", &SignalPath::refractive_index_source_)
        .def_readwrite("refractive_index_destination", &SignalPath::refractive_index_destination_)
        .def_readwrite("emit", &SignalPath::emit_);
}
