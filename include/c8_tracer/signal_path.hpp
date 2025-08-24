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
          receive_(receive), R_distance_(distance),
          fresnelS_(fresnelS), fresnelP_(fresnelP) {}
    ~SignalPath() = default;

    inline std::string to_string() const
    {
      return "SignalPath[emit: " + std::to_string(emit_) +
             " receive: " + std::to_string(receive_) +
             " " + Path::to_string() + "]";
    }
    const_iterator begin() const { return Path::begin(); }
    const_iterator end() const { return Path::end(); }
  };

  inline SignalPath FlipSignalPath(SignalPath const &inPath)
  {
    Path flippedPath(inPath.getEnd());

    for (int ibin = inPath.getNSegments() - 2; ibin >= 0; ibin--)
    {
      flippedPath.addToEnd(inPath.getPoint(ibin));
    }

    return SignalPath(inPath.propagation_time_, inPath.average_refractive_index_,
                      inPath.refractive_index_destination_,
                      inPath.refractive_index_source_, inPath.receive_ * -1,
                      inPath.emit_ * -1, inPath.R_distance_, flippedPath,
                      inPath.fresnelS_, inPath.fresnelP_);
  }

  inline std::ostream &operator<<(std::ostream &os, const SignalPath &sig_path)
  {
    os << sig_path.to_string();
    return os;
  }

} // namespace c8_tracer
