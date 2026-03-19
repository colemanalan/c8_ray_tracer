#pragma once

#include "c8_tracer/environment.hpp"
#include "c8_tracer/logger.hpp"
#include "c8_tracer/signal_path.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"
#include "c8_tracer/transcribed/RayTracer.hpp"
#include "c8_tracer/transcribed/RayTracingTable.hpp"

namespace c8_tracer
{

  inline InterpolationTableGenerator2D::InterpolationTableGenerator2D(RayTracer2D &rayTracer,
                                                                      EnvironmentBase const &env)
      : rayTracer_(rayTracer), axis_(rayTracer.GetAxis().normalized()), env_(env) {}

  inline std::vector<RayTracingTable>
  InterpolationTableGenerator2D::GenerateTables(
      LengthType minR, LengthType maxR, uint const nRBins, LengthType minZ,
      LengthType maxZ, uint const nZBins, Point const &antennaPos, SolutionMethod method) const
  {

    /*
     * Creates a set RayTracingTable objects for a geometry where the solutions
     * are translationally symmetric.
     * A cylinder of solutions is created which has extent minZ to maxZ in length
     * and minR to maxR in radius.
     * the physical space covered by this geometry is given with respect to `antennaPos`
     * and the symmetry axis from the `rayTracer` that is supplied to the constructor
     */

    if (maxR <= minR)
    {
      throw std::runtime_error("Maximum radius is <= minimum radius. Cannot make table");
    }

    if (maxZ <= minZ)
    {
      throw std::runtime_error("Maximum z is <= minimum z. Cannot make table");
    }

    if (nRBins < 2)
    {
      throw std::runtime_error("Must have at least 2 radial bins.");
    }

    if (nZBins < 2)
    {
      throw std::runtime_error("Must have at least 2 z bins.");
    }

    // auto const coordSys = antennaPos.getCoordinateSystem();
    auto const antennaPlane = Plane(antennaPos, axis_);
    auto const refIndex = env_.get_n(antennaPos);

    // TODO: hard coded number of solutions
    std::vector<RayTracingTable> tables;
    tables.reserve(2);
    tables.emplace_back(minR, maxR, nRBins, minZ, maxZ, nZBins, antennaPlane, refIndex);
    tables.emplace_back(minR, maxR, nRBins, minZ, maxZ, nZBins, antennaPlane, refIndex);

    // Get a unit vector that is perpendicular to the symmetry axis
    DirectionVector perpUnit = DirectionVector({0, 1, 0}).cross(axis_);
    if (perpUnit.getNorm() == 0)
    {
      perpUnit = DirectionVector({1, 0, 0}).cross(axis_);
    }
    perpUnit = perpUnit.normalized();

    auto const &start = antennaPos;

    uint constexpr nProgressBars = 50; // Number of steps in the progress bar
    double const totalCells = double(nRBins) * nZBins;
    uint progressSteps = 0;

    auto const startZ = GetZ_(start);
    LOG_INFO("Generating table for antenna at " + std::to_string(antennaPos));
    LOG_INFO("--> min R: " + std::to_string(minR) + ", max R: " + std::to_string(maxR) + ", N: " + std::to_string(nRBins));
    LOG_INFO("--> min Z: " + std::to_string(minZ + startZ) + ", max Z: " + std::to_string(maxZ + startZ) + ", N: " + std::to_string(nZBins));

    // Try once to find unique solutions for each point
    bool atLeastOneFound = false;
    uint nMissing = 0;
    for (uint iz = 0; iz < nZBins; iz++)
    {
      for (uint ir = 0; ir < nRBins; ir++)
      {

        // Loading bar
        uint indexFlat = iz * nRBins + ir;
        uint nextStep = static_cast<uint>((indexFlat / totalCells) * nProgressBars);
        if (nextStep > progressSteps)
        {
          progressSteps = nextStep;
          std::cout << '\r' << "Generated: [" << std::string(progressSteps, '-') << '>'
                    << std::string(nProgressBars - progressSteps, ' ') << "] "
                    << (indexFlat * 100 / int(totalCells)) << "% " << indexFlat << "/"
                    << int(totalCells) << std::flush;
        }

        Point const target =
            start + perpUnit * tables[0].GetR(ir) + axis_ * tables[0].GetZ(iz);

        LOG_DEBUG("Start " + std::to_string(start));
        LOG_DEBUG("Target " + std::to_string(target));
        LOG_DEBUG("N " + std::to_string(env_.get_n(start)));

        std::vector<SignalPath> signalPaths;
        if (method == SolutionMethod::Brent)
          signalPaths = rayTracer_.GetSignalPathsBrent(start, target, env_);

        if (method == SolutionMethod::NGD)
          signalPaths = rayTracer_.GetSignalPathsNGD(start, target, env_);

        if (signalPaths.empty())
        {
          nMissing++;
          continue;
        }

        if (tables.size() < signalPaths.size())
        {
          LOG_DEBUG("Table size (" + std::to_string(tables.size()) +
                    ") is less than the number of solutions (" + std::to_string(signalPaths.size()) + ")");
          throw std::runtime_error("Table/solution size mismatch!");
        }

        const size_t nValid = AnyArrivedAtCorrectLocation_(signalPaths, target);

        if (nValid == 0)
        {
          nMissing++;
          continue;
        }

        if (!AllSolutionsUnique_(signalPaths))
        {
          nMissing++;
          continue;
        }

        atLeastOneFound = true;

        for (size_t i = 0; i < signalPaths.size(); ++i)
        {
          FillTables_(tables[i], ir, iz, signalPaths[i]);
        }
      } // End for r bin
    } // End for z bin

    std::cout << '\r' << "Generated: [" << std::string(nProgressBars, '-') << "] 100%\n";

    if (!atLeastOneFound)
    {
      LOG_WARNING("No valid solutions were found to generate the table. Quitting...");
      return tables;
    }

    if (!nMissing)
    {
      LOG_INFO("All entries created. Table generated.");
      return tables;
    }

    LOG_INFO("After first pass, number of missing entries " + std::to_string(nMissing));

    // Loop 2, look again for solutions by using estimates from nearby nodes

    nMissing = 0;

    for (uint iz = 0; iz < tables[0].GetNZBins(); iz++)
    {

      bool shadowZone = false;

      for (uint ir = 0; ir < tables[0].GetNRBins(); ir++)
      {

        Point const target =
            start + perpUnit * tables[0].GetR(ir) + axis_ * tables[0].GetZ(iz);

        for (size_t isol = 0; isol < tables.size(); isol++)
        {

          // skip existing
          if (tables[isol].GetLength(ir, iz) > 0.0)
          {
            LOG_TRACE("Solution exists at " + std::to_string(ir) + ' ' + std::to_string(iz));
            continue;
          }

          if (shadowZone)
          {
            nMissing++;
            LOG_TRACE("In shadow zone " + std::to_string(ir) + ' ' + std::to_string(iz));
            continue;
          }

          double costheta = EstimateLaunchAngle_(ir, iz, tables[isol]);

          LOG_TRACE("Working on " + std::to_string(tables[0].GetZ(iz)) + ' ' +
                    std::to_string(tables[0].GetR(ir)) + ' ' + std::to_string(costheta));

          if (!costheta)
          {
            nMissing++;
            LOG_TRACE("costheta is zero");
            continue;
          }

          // Construct the seed vector from the launch angle
          costheta = std::min(std::max(costheta, -1.0), 1.0); // clamp on +/-1
          DirectionVector seed({sqrt(1 - costheta * costheta), 0.0, costheta});
          DirectionVector emit = seed; // init values
          DirectionVector receive = seed;
          bool const foundSolution =
              rayTracer_.FindEmitAndReceiveNGD(start, target, env_, seed, emit, receive);

          auto const path = (GetZ_(start) > GetZ_(target))
                                ? rayTracer_.GetSignalPath(start, emit, target, env_)
                                : FlipSignalPath(rayTracer_.GetSignalPath(
                                      target, receive * -1.0, start, env_));

          if (!foundSolution || (path.getEnd() - target).getNorm() > 0.01)
          {

            LOG_DEBUG("Bad path, Z " + std::to_string(tables[0].GetZ(iz)) + " R " +
                      std::to_string(tables[0].GetR(ir)) + ", Sol " + std::to_string(isol) +
                      ", success " + std::string(foundSolution ? "True" : "False") +
                      ", diff " + std::to_string((path.getEnd() - target).getNorm()));
            shadowZone = true; // All further radii on this row will be failures
            nMissing++;
          }
          else
          {
            LOG_TRACE("FILLING TABLE FOR SOLUTION " + std::to_string(isol));
            FillTables_(tables[isol], ir, iz, path);
          }

        } // End for solution
      } // End for r bin
    } // End for z bin

    LOG_INFO("After second pass, number of missing entries " + std::to_string(nMissing / 2));

    return tables;
  }

  inline uint
  InterpolationTableGenerator2D::AnyArrivedAtCorrectLocation_(
      std::vector<SignalPath> const &paths, Point const &target) const
  {
    uint count = 0;
    for (auto const &path : paths)
    {
      count += ArrivedAtCorrectLocation_(path, target);
    }
    return count;
  }

  inline bool
  InterpolationTableGenerator2D::ArrivedAtCorrectLocation_(
      SignalPath const &path, Point const &target, LengthType tolerance) const
  {
    return (path.getEnd() - target).getNorm() < tolerance;
  }

  inline bool InterpolationTableGenerator2D::AllSolutionsUnique_(
      std::vector<SignalPath> const &paths, double tolerance) const
  {
    /*
     * Checks if all the solutions are different from the others. This can happen
     * when the ray tracer finds the same minimum
     */
    for (size_t i = 0; i < paths.size(); i++)
    {
      for (size_t j = i + 1; j < paths.size(); j++)
      {
        if (acos(paths[i].emit_.dot(paths[j].emit_)) < tolerance)
        {
          return false;
        }
      }
    }

    return true;
  }

  inline void InterpolationTableGenerator2D::FillTables_(
      RayTracingTable &table, const uint ir, const uint iz,
      SignalPath const &path) const
  {
    /*
     * Fills the individual tables from the ray tracing solution
     */

    table.SetLaunch(-1 * GetZ_(path.receive_), ir,
                    iz); // swapped to go from table to ant
    table.SetReceive(GetZ_(path.emit_), ir,
                     iz); // swapped to go from table to ant, no neg to point inward
    table.SetLength(path.getLength(), ir, iz);
    table.SetDuration(path.propagation_time_, ir, iz);
    table.SetFresnelS(path.fresnelS_, ir, iz);
    table.SetFresnelP(path.fresnelP_, ir, iz);
  }

  // Helper function to cleanly interpolate over the nearby neighbors of a given bin
  template <typename GetCoord, typename GetValue, typename IsValid>
  static std::optional<double> Interpolate(LengthType v, int i0, int i1,
                                           GetCoord &&getCoord, GetValue &&getValue,
                                           IsValid &&isValid)
  {
    if (!isValid(i0) || !isValid(i1))
      return std::nullopt;
    auto const denom = getCoord(i1) - getCoord(i0);
    if (denom == 0.0)
      return std::nullopt;
    double f = (v - getCoord(i0)) / denom;
    return getValue(i0) * (1 - f) + getValue(i1) * f;
  }

  inline double
  InterpolationTableGenerator2D::EstimateLaunchAngle_(
      uint ir, uint iz, RayTracingTable &table) const
  {

    /*
     * Fairly simplistic interpolation scheme for the launch angle. If possible,
     * the solution will be found using the surrounding nodes. Otherwise, any
     * nearby values will be used instead. The average of all these bilinear
     * interpolations are used as the final guess
     */

    struct Estimate
    {
      double value = 0.0;
      double weight = 0.0;
      void add(std::optional<double> v, double w = 1.0)
      {
        if (v)
        {
          value += *v * w;
          weight += w;
        }
      }
    };

    Estimate est;

    auto getR = [&](int i)
    { return table.GetR(i); };
    auto getZ = [&](int i)
    { return table.GetZ(i); };
    auto getLaunchR = [&](int i)
    { return table.GetLaunch(i, iz); };
    auto getLaunchZ = [&](int i)
    { return table.GetLaunch(ir, i); };
    auto validR = [&](int i)
    { return table.IsValid(i, iz); };
    auto validZ = [&](int i)
    { return table.IsValid(ir, i); };

    // r-direction
    est.add(Interpolate(table.GetR(ir), ir - 1, ir + 1, getR, getLaunchR, validR));
    est.add(Interpolate(table.GetR(ir), ir - 2, ir - 1, getR, getLaunchR, validR));
    est.add(Interpolate(table.GetR(ir), ir + 1, ir + 2, getR, getLaunchR, validR));

    // z-direction
    est.add(Interpolate(table.GetZ(iz), iz - 1, iz + 1, getZ, getLaunchZ, validZ));
    est.add(Interpolate(table.GetZ(iz), iz - 2, iz - 1, getZ, getLaunchZ, validZ));
    est.add(Interpolate(table.GetZ(iz), iz + 1, iz + 2, getZ, getLaunchZ, validZ));

    return est.weight > 0.0 ? est.value / est.weight : 0.0;
  }

} // namespace c8_tracer
