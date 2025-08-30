#pragma once

#include <vector>

#include "c8_tracer/plane.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"

namespace c8_tracer
{

  class RayTracingTable
  {
  public:
    RayTracingTable(LengthType minR, LengthType maxR, uint nRBins, LengthType minZ,
                    LengthType maxZ, uint nZBins, Plane const &plane);
    ~RayTracingTable();

    // Setters
    void SetLaunch(double launch, uint ir, uint iz);
    void SetReceive(double receive, uint ir, uint iz);
    void SetLength(LengthType length, uint ir, uint iz);
    void SetDuration(TimeType duration, uint ir, uint iz);
    void SetFresnelS(double fresnelS, uint ir, uint iz);
    void SetFresnelP(double fresnelP, uint ir, uint iz);

    // Getters
    double GetLaunch(uint ir, uint iz) const;
    double GetReceive(uint ir, uint iz) const;
    LengthType GetLength(uint ir, uint iz) const;
    TimeType GetDuration(uint ir, uint iz) const;
    double GetFresnelS(uint ir, uint iz) const;
    double GetFresnelP(uint ir, uint iz) const;
    LengthType GetR(uint ir) const { return r_[ir]; }
    LengthType GetZ(uint iz) const { return z_[iz]; }
    uint GetNRBins() const { return nRBins_; }
    uint GetNZBins() const { return nZBins_; }
    std::tuple<uint, uint> GetIrIz(LengthType const &r, LengthType const &z) const;
    auto const &GetPlane() const { return plane_; }

    bool IsValid(uint ir, uint iz, bool warn) const;

    // Interpolators
    double GetLaunch(LengthType const &r, LengthType const &z) const;
    double GetReceive(LengthType const &r, LengthType const &z) const;
    LengthType GetLength(LengthType const &r, LengthType const &z) const;
    TimeType GetDuration(LengthType const &r, LengthType const &z) const;
    double GetFresnelS(LengthType const &r, LengthType const &z) const;
    double GetFresnelP(LengthType const &r, LengthType const &z) const;

    // Utility
    bool ContainsPoint(Point const &x) const;
    SignalPath GetSignalPath(Point const &x0, double n0, double nf) const;

    void SetIndexOfRefraction(double x) { indexOfRefraction_ = x; }
    double GetIndexOfRefraction() const { return indexOfRefraction_; }

    void Print() const;
    void ResetTables();

  private:
    // Grid dimensions
    const uint nRBins_;
    const uint nZBins_;
    const LengthType maxR_, minR_, maxZ_, minZ_;
    InverseLengthType inverseDR_ = 1.0;
    InverseLengthType inverseDZ_ = 1.0;

    // Coordinate grids
    std::vector<LengthType> r_;
    std::vector<LengthType> z_;

    // Flat table accessors
    inline size_t index(uint ir, uint iz) const { return iz * nRBins_ + ir; }

    // Tables
    std::vector<double> launch_;
    std::vector<double> receive_;
    std::vector<LengthType> length_;
    std::vector<TimeType> duration_;
    std::vector<double> fresnelS_;
    std::vector<double> fresnelP_;

    Plane const plane_;
    double indexOfRefraction_ = 1.0;

    mutable long int totalCalls_ = 0;
    mutable long int untrackedCalls_ = 0;

    // Interpolation methods
    template <typename T>
    T InterpolateFromIndex_(std::vector<T> const &table, LengthType const &r,
                            LengthType const &z, uint ir, uint iz) const;

    template <typename T>
    T Interpolate_(std::vector<T> const &table, LengthType const &r,
                   LengthType const &z) const;

    std::tuple<LengthType, LengthType> GetRAndZ_(Point const &x) const;

    template <typename T1, typename T2>
    void PrintTable_(T1 const &table, T2 const unit) const;

    uint SigFigs_(double val, uint decimalPoints) const;
  };

} // namespace c8_tracer

#include "c8_tracer/transcribed/RayTracingTable.inl"
