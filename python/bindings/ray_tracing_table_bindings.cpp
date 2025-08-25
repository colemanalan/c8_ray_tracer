#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "c8_tracer/signal_path.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"
#include "c8_tracer/transcribed/RayTracingTable.hpp"
#include "c8_tracer/transcribed/InterpolationTableGenerator2D.hpp"

namespace py = pybind11;
using namespace c8_tracer;

void bind_ray_tracing_table(py::module_ &m)
{

     py::module_ tables = m.def_submodule("tables", "Tables to hold pre-calculated ray-tracing tables");

     py::class_<RayTracingTable>(tables, "RayTracingTable")
         .def(py::init<LengthType, LengthType, uint, LengthType, LengthType, uint, const Plane &>(),
              py::arg("minR"), py::arg("maxR"), py::arg("nRBins"),
              py::arg("minZ"), py::arg("maxZ"), py::arg("nZBins"),
              py::arg("plane"))

         // Setters
         .def("SetLaunch", &RayTracingTable::SetLaunch)
         .def("SetReceive", &RayTracingTable::SetReceive)
         .def("SetLength", &RayTracingTable::SetLength)
         .def("SetDuration", &RayTracingTable::SetDuration)
         .def("SetFresnelS", &RayTracingTable::SetFresnelS)
         .def("SetFresnelP", &RayTracingTable::SetFresnelP)

         // Getters (by index)
         .def("GetLaunch", py::overload_cast<uint, uint>(&RayTracingTable::GetLaunch, py::const_))
         .def("GetReceive", py::overload_cast<uint, uint>(&RayTracingTable::GetReceive, py::const_))
         .def("GetLength", py::overload_cast<uint, uint>(&RayTracingTable::GetLength, py::const_))
         .def("GetDuration", py::overload_cast<uint, uint>(&RayTracingTable::GetDuration, py::const_))
         .def("GetFresnelS", py::overload_cast<uint, uint>(&RayTracingTable::GetFresnelS, py::const_))
         .def("GetFresnelP", py::overload_cast<uint, uint>(&RayTracingTable::GetFresnelP, py::const_))

         // Getters (by value)
         .def("GetLaunch", py::overload_cast<const LengthType &, const LengthType &>(&RayTracingTable::GetLaunch, py::const_))
         .def("GetReceive", py::overload_cast<const LengthType &, const LengthType &>(&RayTracingTable::GetReceive, py::const_))
         .def("GetLength", py::overload_cast<const LengthType &, const LengthType &>(&RayTracingTable::GetLength, py::const_))
         .def("GetDuration", py::overload_cast<const LengthType &, const LengthType &>(&RayTracingTable::GetDuration, py::const_))
         .def("GetFresnelS", py::overload_cast<const LengthType &, const LengthType &>(&RayTracingTable::GetFresnelS, py::const_))
         .def("GetFresnelP", py::overload_cast<const LengthType &, const LengthType &>(&RayTracingTable::GetFresnelP, py::const_))

         // Grid access
         .def("GetR", &RayTracingTable::GetR)
         .def("GetZ", &RayTracingTable::GetZ)
         .def("GetNRBins", &RayTracingTable::GetNRBins)
         .def("GetNZBins", &RayTracingTable::GetNZBins)
         .def("GetIrIz", &RayTracingTable::GetIrIz)
         .def("GetPlane", &RayTracingTable::GetPlane, py::return_value_policy::reference)

         // Utility
         .def("ContainsPoint", &RayTracingTable::ContainsPoint)
         .def("GetSignalPath", &RayTracingTable::GetSignalPath)
         .def("SetIndexOfRefraction", &RayTracingTable::SetIndexOfRefraction)
         .def("GetIndexOfRefraction", &RayTracingTable::GetIndexOfRefraction)
         .def("Print", &RayTracingTable::Print)
         .def("ResetTables", &RayTracingTable::ResetTables);
}

void bind_interpolation_table_generator_2d(py::module_ &m)
{
     py::module_ tables = m.def_submodule("tables", "Tables to hold pre-calculated ray-tracing tables");

     py::class_<InterpolationTableGenerator2D>(tables, "InterpolationTableGenerator2D")
         .def(py::init<RayTracerBase &, const EnvironmentBase &>(),
              py::arg("ray_tracer"), py::arg("environment"),
              py::keep_alive<1, 2>(), py::keep_alive<1, 3>()) // Ensure references stay alive

         .def("GenerateTables", &InterpolationTableGenerator2D::GenerateTables,
              py::arg("minR"), py::arg("maxR"), py::arg("nRBins"),
              py::arg("minZ"), py::arg("maxZ"), py::arg("nZBins"),
              py::arg("antennaPos"));
}
