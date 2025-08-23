#pragma once

#include "c8_tracer/vec3.hpp"
#include "c8_tracer/transcribed/Path.hpp"

namespace c8_tracer
{

  class SignalPath : public Path
  {
  public:
    double propagation_time_ = 0.0;
    double average_refractive_index_ = 0.0;
    double refractive_index_source_ = 0.0;
    double refractive_index_destination_ = 0.0;
    Vec3 emit_;
    Vec3 receive_;
    double R_distance_;
    double fresnelS_;
    double fresnelP_;

    SignalPath() = default;
    SignalPath(double dt, double avg_n, double n_source, double n_dest,
               Vec3 const &emit, Vec3 const &receive, double distance,
               Path const &points, double fresnelS, double fresnelP)
        : Path(points), propagation_time_(dt), average_refractive_index_(avg_n),
          refractive_index_source_(n_source),
          refractive_index_destination_(n_dest), emit_(emit),
          R_distance_(distance), receive_(receive),
          fresnelS_(fresnelS), fresnelP_(fresnelP) {}
    ~SignalPath() = default;
  };

} // namespace c8_tracer
