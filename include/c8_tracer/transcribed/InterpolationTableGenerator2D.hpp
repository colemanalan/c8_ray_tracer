#pragma once

#include "c8_tracer/signal_path.hpp"
#include "c8_tracer/environment.hpp"
#include "c8_tracer/transcribed/RayTracer.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"
#include "c8_tracer/transcribed/RayTracingTable.hpp"

namespace c8_tracer
{

  class InterpolationTableGenerator2D
  {
  public:
    /**
     * Generates `RayTracingTable` for specified geometries using the general ray-tracing infrastructure.
     * Like the ray-tracing tables that it creates, this class assumes that the media only varies in
     * one dimension (z) such that a 2D table is sufficient to descibe a volume
     *
     * @param rayTracer instance of a ray tracer to be used
     * @param env the environment object that can be queried for n and Grad[n]
     */
    InterpolationTableGenerator2D(RayTracer2D &rayTracer, EnvironmentBase const &env);
    ~InterpolationTableGenerator2D() {}

    /**
     * Generates a vector of ray-tracing tables for the specified volume and binning.
     * The solutions are found in two separate passes. In the first, the table is blindly filled
     * using the chosen solution-finding method. In the second pass, an attempt is made to
     * fill in missing values by using the solutions from neighboring points as a seed. The second
     * pass will always use numerical gradient descent, regardless of the chosen `method` type
     *
     * Currently assumes that exactly two solutions exist!
     *
     * @param minR minimum radius of the table
     * @param maxR maximum radius of the table
     * @param nRBins number of bins in the r-direction
     * @param minZ minimum z-value of the table
     * @param maxZ maximum z-value of the table
     * @param nZBins number of bins in the z-direction
     * @param centerPos center point that defines where the (z=0, r=0) corresponds to
     * @param method the method by which the first round of solutions is found
     *
     * @return vector of filled RayTracingTable objects, one for each ray solution
     */
    std::vector<RayTracingTable> GenerateTables(LengthType minR, LengthType maxR,
                                                uint const nRBins, LengthType minZ,
                                                LengthType maxZ, uint const nZBins,
                                                Point const &centerPos, SolutionMethod method) const;

  private:
    // Returns the number of paths that end at the `target`
    uint AnyArrivedAtCorrectLocation_(std::vector<SignalPath> const &paths,
                                      Point const &target) const;
    // does the signal path end at `target` within `tolerance`
    bool ArrivedAtCorrectLocation_(SignalPath const &path, Point const &target,
                                   LengthType tolerance = 0.1) const;

    // no solutions are the same within tolerance
    bool AllSolutionsUnique_(std::vector<SignalPath> const &paths,
                             double tolerance = 1e-5) const;

    // fill the table with the values in the signal path
    void FillTables_(RayTracingTable &table, const uint ir,
                     const uint iz, SignalPath const &path) const;

    // use the current entries in the table to find the solution
    double EstimateLaunchAngle_(uint ir, uint iz, RayTracingTable &table) const;

    RayTracer2D &rayTracer_;
    DirectionVector const axis_;
    EnvironmentBase const &env_;

    LengthType GetZ_(Point const &p) const { return rayTracer_.GetAxis().dot(p); }
  };
} // namespace c8_tracer

#include "c8_tracer/transcribed/InterpolationTableGenerator2D.inl"
