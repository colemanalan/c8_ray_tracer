#pragma once

#include <vector>

#include "c8_tracer/plane.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"

namespace c8_tracer
{

  class RayTracingTable
  {
  public:
    /**
     * Container to hold pre-calculated ray paths from many points in space to a single location.
     * The table covers an equally spaced grid of points in radius and z for a medium that varies
     * only in one dimension. This dimension is denoted z and r is the distance in the other two
     * dimensions. Tables are created for all the values that are needed to construct a `SignalPath`
     *
     * The table can be queried for any arbitrary point within the bounds. The result will be a
     * binlinear interpolation using the four neighbors (+/- z and +/- r). This table does not
     * extrapolate and will return a NaN if any of the neighbors include a point that is not
     * filled. So this table will also not interpolate within max(dR, dZ) of a caustic
     *
     * @param minR minimum radius of the table
     * @param maxR maximum radius of the table
     * @param nRBins number of bins in the r-direction
     * @param minZ minimum z-value of the table
     * @param maxZ maximum z-value of the table
     * @param nZBins number of bins in the z-direction
     * @param plane Plane object with a center point that defines where the (z=0, r=0) corresponds
     * to and normal vector that defines the z-hat direction
     * @param nCenter the index of refraction corresponding to the center point of the coordinate
     * system (i.e. at `plane.getCenter()`)
     */
    RayTracingTable(LengthType minR, LengthType maxR, uint nRBins, LengthType minZ,
                    LengthType maxZ, uint nZBins, Plane const &plane, double nCenter);
    RayTracingTable(RayTracingTable &&) noexcept = default;
    RayTracingTable &operator=(RayTracingTable &&) noexcept = default;

    ~RayTracingTable();

    // Setters of interpolation-node values

    void SetLaunch(double launch, uint ir, uint iz);
    void SetReceive(double receive, uint ir, uint iz);
    void SetLength(LengthType length, uint ir, uint iz);
    void SetDuration(TimeType duration, uint ir, uint iz);
    void SetFresnelS(double fresnelS, uint ir, uint iz);
    void SetFresnelP(double fresnelP, uint ir, uint iz);

    // Getters of interpolation-node values

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

    // check to see if `ir` and `iz` are valid indicies
    bool IsValid(uint ir, uint iz, bool warn) const;

    // Interpolators

    double GetLaunch(LengthType const &r, LengthType const &z) const;
    double GetReceive(LengthType const &r, LengthType const &z) const;
    LengthType GetLength(LengthType const &r, LengthType const &z) const;
    TimeType GetDuration(LengthType const &r, LengthType const &z) const;
    double GetFresnelS(LengthType const &r, LengthType const &z) const;
    double GetFresnelP(LengthType const &r, LengthType const &z) const;

    // Utility

    // checks if the point `x` is contained within the table
    bool ContainsPoint(Point const &x) const;
    /**
     * Gets the signal path for the ray starting at `x0` and ending at
     * the center point of the table (i.e. `plane.getCenter()`). The parameters
     * of the SignalPath will be interpolated from the table values.
     *
     * @param x0 starting location of the path
     * @param n0 index of refraction at `x0` (needed to fill signal path)
     *
     * @return interpolated signal path
     */
    SignalPath GetSignalPath(Point const &x0, double n0) const;

    double GetIndexOfRefraction() const { return indexOfRefraction_; }

    // Prints the table values
    void Print() const;
    // Reset the entries of the table to the default values (NaN)
    void ResetTables();

  private:
    // Grid dimensions

    uint const nRBins_;
    uint const nZBins_;
    LengthType const maxR_, minR_, maxZ_, minZ_;
    InverseLengthType inverseDR_ = 1.0;
    InverseLengthType inverseDZ_ = 1.0;

    // Coordinate grids

    std::vector<LengthType> r_;
    std::vector<LengthType> z_;

    // Flat table accessors
    inline size_t index(uint ir, uint iz) const noexcept { return iz * nRBins_ + ir; }

    // Tables

    std::unique_ptr<double[]> launch_;
    std::unique_ptr<double[]> receive_;
    std::unique_ptr<LengthType[]> length_;
    std::unique_ptr<TimeType[]> duration_;
    std::unique_ptr<double[]> fresnelS_;
    std::unique_ptr<double[]> fresnelP_;

    Plane const plane_;
    double const indexOfRefraction_ = 1.0;

    mutable long int totalCalls_ = 0;
    mutable long int untrackedCalls_ = 0;

    template <typename T>
    inline void FillWithNaN_(std::unique_ptr<T[]> &ptr, size_t N);

    // Interpolation methods

    template <typename T>
    T InterpolateFromIndex_(T const *table, LengthType const &r,
                            LengthType const &z, uint ir, uint iz) const noexcept;

    template <typename T>
    T Interpolate_(T const *table, LengthType const &r,
                   LengthType const &z) const noexcept;

    std::tuple<LengthType, LengthType> GetRAndZ_(Point const &x) const;

    template <typename T1, typename T2>
    void PrintTable_(T1 const &table, T2 const unit) const;

    uint SigFigs_(double val, uint decimalPoints) const;
  };

} // namespace c8_tracer

#include "c8_tracer/transcribed/RayTracingTable.inl"
