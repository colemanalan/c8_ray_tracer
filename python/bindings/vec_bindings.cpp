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
         .def(py::init<double, double>())
         .def_readwrite("x", &Vec2::x)
         .def_readwrite("y", &Vec2::y)
         .def("dot", &Vec2::dot)
         .def("length", &Vec2::length)
         .def("norm", &Vec2::norm)
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
         .def("__setitem__", [](Vec2 &v, int i, double value)
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
         .def(py::self * double())
         .def(py::self / double())
         .def("__repr__", [](const Vec2 &v)
              { return "Vec2(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ")"; });

     // Vec3
     py::class_<Vec3>(m, "Vec3", "3D vector")
         .def(py::init<double, double, double>())
         .def(py::init([](py::iterable iterable)
                       {
        std::vector<double> values;
        for (auto item : iterable)
            values.push_back(py::cast<double>(item));
        if (values.size() != 3)
            throw std::runtime_error("Vec3 requires an iterable of length 3");
        return Vec3(values[0], values[1], values[2]); }),
              "Initialize from any iterable of 3 numbers")
         .def_readwrite("x", &Vec3::x)
         .def_readwrite("y", &Vec3::y)
         .def_readwrite("z", &Vec3::z)
         .def("dot", &Vec3::dot)
         .def("cross", &Vec3::cross)
         .def("length", &Vec3::length)
         .def("norm", &Vec3::norm)
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
         .def("__setitem__", [](Vec3 &v, int i, double value)
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
         .def(py::self * double())
         .def(py::self / double())
         .def("__repr__", [](const Vec3 &v)
              { return "Vec3(" + std::to_string(v.x) +
                       ", " + std::to_string(v.y) +
                       ", " + std::to_string(v.z) + ")"; });
}