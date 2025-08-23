#include <pybind11/pybind11.h>
#include "c8_tracer/logger.hpp"

namespace py = pybind11;
using namespace c8_tracer;

void bind_logger(py::module_ &m)
{

    // Create a 'logging' submodule
    py::module_ logging = m.def_submodule("logging", "Logging utilities");

    // Bind LogLevel enum
    py::enum_<c8_tracer::LogLevel>(logging, "LogLevel")
        .value("DEBUG", c8_tracer::LogLevel::DEBUG)
        .value("INFO", c8_tracer::LogLevel::INFO)
        .value("WARNING", c8_tracer::LogLevel::WARNING)
        .value("ERROR", c8_tracer::LogLevel::ERROR);

    // Bind only the global logger instance
    py::class_<c8_tracer::Logger>(logging, "_Logger")
        .def("set_level", &c8_tracer::Logger::set_level)
        .def("log", &c8_tracer::Logger::log)
        .def("debug", &c8_tracer::Logger::debug)
        .def("info", &c8_tracer::Logger::info)
        .def("warning", &c8_tracer::Logger::warning)
        .def("error", &c8_tracer::Logger::error);

    // Expose the singleton as a module-level attribute
    logging.attr("logger") = py::cast(&c8_tracer::logger, py::return_value_policy::reference);
}