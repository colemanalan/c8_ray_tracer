#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include "c8_tracer/vec2.hpp"
#include "c8_tracer/vec3.hpp"

namespace py = pybind11;
using namespace c8_tracer;

void bind_vec(py::module_ &m)
{
    // Vec2
    py::class_<Vec2>(m, "Vec2")
        .def(py::init<float, float>())
        .def_readwrite("x", &Vec2::x)
        .def_readwrite("y", &Vec2::y)
        .def("dot", &Vec2::dot)
        .def("length", &Vec2::length)
        .def("normalized", &Vec2::normalized)
        .def("__eq__", [](const Vec2 &a, const Vec2 &b)
             { return a.x == b.x && a.y == b.y; })
        .def("__neg__", [](const Vec2 &v)
             { return Vec2(-v.x, -v.y); })
        .def("__getitem__", [](const Vec2 &v, int i)
             {
            if (i == 0) return v.x;
            if (i == 1) return v.y;
            throw py::index_error("Vec2 index out of range"); })
        .def("__setitem__", [](Vec2 &v, int i, float value)
             {
            if (i == 0) v.x = value;
            else if (i == 1) v.y = value;
            else throw py::index_error("Vec2 index out of range"); })
        .def("__len__", [](const Vec2 &)
             { return 2; })
        .def("__iter__", [](const Vec2 &v)
             { return py::make_iterator(&v.x, &v.x + 2); }, py::keep_alive<0, 1>())
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * float())
        .def(py::self / float())
        .def("__repr__", [](const Vec2 &v)
             { return "<Vec2 x=" + std::to_string(v.x) + " y=" + std::to_string(v.y) + ">"; });

    // Vec3
    py::class_<Vec3>(m, "Vec3")
        .def(py::init<float, float, float>())
        .def_readwrite("x", &Vec3::x)
        .def_readwrite("y", &Vec3::y)
        .def_readwrite("z", &Vec3::z)
        .def("dot", &Vec3::dot)
        .def("cross", &Vec3::cross)
        .def("length", &Vec3::length)
        .def("normalized", &Vec3::normalized)
        .def("__eq__", [](const Vec3 &a, const Vec3 &b)
             { return a.x == b.x && a.y == b.y && a.z == b.z; })
        .def("__neg__", [](const Vec3 &v)
             { return Vec3(-v.x, -v.y, -v.z); })
        .def("__getitem__", [](const Vec3 &v, int i)
             {
            if (i == 0) return v.x;
            if (i == 1) return v.y;
            if (i == 2) return v.z;
            throw py::index_error("Vec3 index out of range"); })
        .def("__setitem__", [](Vec3 &v, int i, float value)
             {
            if (i == 0) v.x = value;
            else if (i == 1) v.y = value;
            else if (i == 2) v.z = value;
            else throw py::index_error("Vec3 index out of range"); })
        .def("__len__", [](const Vec3 &)
             { return 3; })
        .def("__iter__", [](const Vec3 &v)
             { return py::make_iterator(&v.x, &v.x + 3); }, py::keep_alive<0, 1>())
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * float())
        .def(py::self / float())
        .def("__repr__", [](const Vec3 &v)
             { return "<Vec3 x=" + std::to_string(v.x) +
                      " y=" + std::to_string(v.y) +
                      " z=" + std::to_string(v.z) + ">"; });
}