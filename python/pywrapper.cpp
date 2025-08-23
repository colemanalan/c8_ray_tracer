#include <pybind11/pybind11.h>

namespace py = pybind11;

// Forward declarations
void bind_cash_karp(py::module_ &);
void bind_environment(py::module_ &);
void bind_logger(py::module_ &);
void bind_tracer(py::module_ &);
void bind_vec(py::module_ &);

PYBIND11_MODULE(c8_tracer_py, m)
{
    m.doc() = "Python bindings for c8_tracer";
    bind_cash_karp(m);
    bind_environment(m);
    bind_logger(m);
    bind_tracer(m);
    bind_vec(m);
}