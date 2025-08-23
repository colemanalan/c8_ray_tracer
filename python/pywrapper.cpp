#include <pybind11/pybind11.h>

namespace py = pybind11;

// Forward declarations
void bind_vec(py::module_&);
void bind_tracer(py::module_&);

PYBIND11_MODULE(c8_tracer_py, m) {
    m.doc() = "Python bindings for c8_tracer";
    bind_vec(m);
    bind_tracer(m);
}