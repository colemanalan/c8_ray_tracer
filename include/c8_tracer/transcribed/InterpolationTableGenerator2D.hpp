#pragma once

#include "c8_tracer/signal_path.hpp"
#include "c8_tracer/environment.hpp"
#include "c8_tracer/ray_tracer_base.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"
#include "c8_tracer/transcribed/RayTracingTable.hpp"

namespace c8_tracer
{
  class InterpolationTableGenerator2D
  {
  public:
    InterpolationTableGenerator2D(RayTracerBase &rayTracer, EnvironmentBase const &env);
    ~InterpolationTableGenerator2D() {}

    std::vector<RayTracingTable> GenerateTables(LengthType minR, LengthType maxR,
                                                uint const nRBins, LengthType minZ,
                                                LengthType maxZ, uint const nZBins,
                                                Point const &antennaPos) const;

  private:
    uint AnyArrivedAtCorrectLocation_(std::vector<SignalPath> const &paths,
                                      Point const &target) const;
    bool ArrivedAtCorrectLocation_(SignalPath const &path, Point const &target,
                                   LengthType tolerance = 0.1) const;

    bool AllSolutionsUnique_(std::vector<SignalPath> const &paths,
                             double tolerance = 1e-5) const;

    void FillTables_(RayTracingTable &table, const uint ir,
                     const uint iz, SignalPath const &path) const;

    double EstimateLaunchAngle_(uint ir, uint iz, RayTracingTable &table) const;

    RayTracerBase &rayTracer_;
    DirectionVector const axis_;
    EnvironmentBase const &env_;

    LengthType GetZ_(Point const &p) const { return rayTracer_.GetAxis().dot(p); }
  };
} // namespace c8_tracer

#include "c8_tracer/transcribed/InterpolationTableGenerator2D.inl"
