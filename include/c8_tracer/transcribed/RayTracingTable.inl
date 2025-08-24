#pragma once

#include <iostream>
#include "c8_tracer/transcribed/c8_typedefs.hpp"
#include "c8_tracer/logger.hpp"

namespace c8_tracer
{

  inline RayTracingTable::RayTracingTable(LengthType minR, LengthType maxR, uint nRBins,
                                          LengthType minZ, LengthType maxZ, uint nZBins,
                                          Plane const &plane)
      : nRBins_(nRBins), nZBins_(nZBins), maxR_(maxR), minR_(minR), maxZ_(maxZ), minZ_(minZ), r_(nRBins), z_(nZBins), launch_(nRBins * nZBins), receive_(nRBins * nZBins), length_(nRBins * nZBins), duration_(nRBins * nZBins), fresnelS_(nRBins * nZBins), fresnelP_(nRBins * nZBins), plane_(plane)
  {

    if (nRBins_ < 2)
    {
      logger.error("Initializing with radius bins n = " + std::to_string(nRBins_));
      throw std::invalid_argument("Must have >1 bins");
    }
    if (nZBins_ < 2)
    {
      logger.error("Initializing with z bins n = " + std::to_string(nZBins_));
      throw std::invalid_argument("Must have >1 bins");
    }

    if (maxZ_ <= minZ_)
    {
      logger.error("Max z " + std::to_string(maxZ_) + ", Min z " + std::to_string(minZ_));
      throw std::invalid_argument("Required: max > min");
    }

    if (maxR_ <= minR_)
    {
      logger.error("Max r " + std::to_string(maxR_) + ", Min r " + std::to_string(minR_));
      throw std::invalid_argument("Required: max > min");
    }

    if (minR_ < 0.0)
    {
      logger.error("Min r " + std::to_string(minR_) + " must be > 0");
      throw std::invalid_argument("Min r must be > 0");
    }

    LengthType dr = (maxR_ - minR_) / (nRBins_ - 1);
    LengthType dz = (maxZ_ - minZ_) / (nZBins_ - 1);
    inverseDR_ = 1.0 / dr;
    inverseDZ_ = 1.0 / dz;

    for (uint i = 0; i < nRBins_; ++i)
    {
      r_[i] = minR_ + i * dr;
    }
    for (uint i = 0; i < nZBins_; ++i)
    {
      z_[i] = minZ_ + i * dz;
    }

    totalCalls_ = 0;
    untrackedCalls_ = 0;

    ResetTables();
  }

  inline RayTracingTable::~RayTracingTable()
  {
    if (untrackedCalls_)
    {
      double const percentUntracked =
          static_cast<float>(untrackedCalls_) / static_cast<float>(totalCalls_) * 100;
      logger.info("There were " + std::to_string(percentUntracked) + " % of tracks that were not contained within this RayTracingTable");
      logger.info("Total tracks: " + std::to_string(totalCalls_) + ", not contained tracks: " + std::to_string(untrackedCalls_));
      if (percentUntracked > 0.2)
      {
        logger.warning(
            "This percentage is large and may indicate that the table"
            "is not large enough to contain your cascade");
      }
    }
  }

  inline void RayTracingTable::ResetTables()
  {
    std::fill(launch_.begin(), launch_.end(), 0.0);
    std::fill(receive_.begin(), receive_.end(), 0.0);
    std::fill(length_.begin(), length_.end(), 0.0);
    std::fill(duration_.begin(), duration_.end(), 0.0);
    std::fill(fresnelS_.begin(), fresnelS_.end(), 1.0);
    std::fill(fresnelP_.begin(), fresnelP_.end(), 1.0);

    totalCalls_ = 0;
    untrackedCalls_ = 0;
  }

  ///////////////////////
  /// Setters
  ///////////////////////

  inline void RayTracingTable::SetLaunch(double launch, uint ir, uint iz)
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      logger.error("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                   " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    launch_[index(ir, iz)] = launch;
  }

  inline void RayTracingTable::SetReceive(double receive, uint ir, uint iz)
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      logger.error("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                   " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    receive_[index(ir, iz)] = receive;
  }

  inline void RayTracingTable::SetLength(LengthType length, uint ir, uint iz)
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      logger.error("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                   " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    length_[index(ir, iz)] = length;
  }

  inline void RayTracingTable::SetDuration(TimeType duration, uint ir, uint iz)
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      logger.error("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                   " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    duration_[index(ir, iz)] = duration;
  }

  inline void RayTracingTable::SetFresnelS(double fresnelS, uint ir, uint iz)
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      logger.error("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                   " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    fresnelS_[index(ir, iz)] = fresnelS;
  }

  inline void RayTracingTable::SetFresnelP(double fresnelP, uint ir, uint iz)
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      logger.error("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                   " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    fresnelP_[index(ir, iz)] = fresnelP;
  }

  ///////////////////////
  /// Getters
  ///////////////////////

  inline double RayTracingTable::GetLaunch(uint ir, uint iz) const
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      logger.error("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                   " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    return launch_[index(ir, iz)];
  }

  inline double RayTracingTable::GetReceive(uint ir, uint iz) const
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      logger.error("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                   " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    return receive_[index(ir, iz)];
  }

  inline LengthType RayTracingTable::GetLength(uint ir, uint iz) const
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      logger.error("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                   " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    return length_[index(ir, iz)];
  }

  inline TimeType RayTracingTable::GetDuration(uint ir, uint iz) const
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      logger.error("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                   " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    return duration_[index(ir, iz)];
  }

  inline double RayTracingTable::GetFresnelS(uint ir, uint iz) const
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      logger.error("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                   " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    return fresnelS_[index(ir, iz)];
  }

  inline double RayTracingTable::GetFresnelP(uint ir, uint iz) const
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      logger.error("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                   " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    return fresnelP_[index(ir, iz)];
  }

  inline bool RayTracingTable::ContainsPoint(Point const &x0) const
  {
    auto const [r, z] = GetRAndZ_(x0);
    if (r < r_[0] || r > r_[nRBins_ - 1] || z < z_[0] || z > z_[nZBins_ - 1])
    {
      return false;
    }
    return true;
  }

  inline SignalPath RayTracingTable::GetSignalPath(Point const &x0, double n0,
                                                   double nf) const
  {
    totalCalls_++;
    // Do not extrapolate
    if (!ContainsPoint(x0))
    {
      untrackedCalls_++;
      return SignalPath(std::numeric_limits<double>::infinity(), 1.0, n0, nf,
                        plane_.getNormal(), plane_.getNormal(),
                        std::numeric_limits<double>::infinity(), Path(x0), 1.0,
                        1.0);
    }

    auto const [r, z] = GetRAndZ_(x0);
    auto const [ir, iz] = GetIrIz(r, z);
    auto const length = InterpolateFromIndex_(length_, r, z, ir, iz);
    auto const duration = InterpolateFromIndex_(duration_, r, z, ir, iz);
    auto const fresnelS = InterpolateFromIndex_(fresnelS_, r, z, ir, iz);
    auto const fresnelP = InterpolateFromIndex_(fresnelP_, r, z, ir, iz);

    auto const launch = InterpolateFromIndex_(launch_, r, z, ir, iz);
    auto const receive = InterpolateFromIndex_(receive_, r, z, ir, iz);

    auto const dr = (plane_.getCenter() - x0);
    auto const xyComponents =
        (dr - plane_.getNormal() * dr.dot(plane_.getNormal())).normalized();

    auto const startDir =
        plane_.getNormal() * launch + xyComponents * sqrt(1 - launch * launch);
    auto const endDir =
        plane_.getNormal() * receive - xyComponents * sqrt(1 - receive * receive);

    return SignalPath(duration,
                      duration * constants::c / length, // average n
                      n0, nf, startDir, endDir, length, Path(x0), fresnelS, fresnelP);
  }

  inline void RayTracingTable::Print() const
  {
    std::cout << "\nTables for observer at " << plane_.getCenter();
    std::cout << "\n";
    std::cout << "Path lengths / m\n";
    PrintTable_(length_, 1.0);
    std::cout << "\n";
    std::cout << "Duration / ns\n";
    PrintTable_(duration_, 1.0);
    std::cout << "\n";
    std::cout << "Launch vectors (cos)\n";
    PrintTable_(launch_, 1.0);
    std::cout << "\n";
    std::cout << "Receive vectors (cos)\n";
    PrintTable_(receive_, 1.0);
    std::cout << "\n";
    std::cout << "Fresnel S\n";
    PrintTable_(fresnelS_, 1.0);
    std::cout << "\n";
    std::cout << "Fresnel P\n";
    PrintTable_(fresnelP_, 1.0);
  }

  inline std::tuple<uint, uint> RayTracingTable::GetIrIz(LengthType const &r,
                                                         LengthType const &z) const
  {
    uint ir = 0;
    if (r >= r_.front())
    {
      ir = static_cast<uint>((r - r_.front()) * inverseDR_);
      ir = std::min(ir, nRBins_ - 2);
    }

    uint iz = 0;
    if (z >= z_.front())
    {
      iz = static_cast<uint>((z - z_.front()) * inverseDZ_);
      iz = std::min(iz, nZBins_ - 2);
    }

    return {ir, iz};
  }

  template <typename T>
  inline T RayTracingTable::InterpolateFromIndex_(std::vector<T> const &table,
                                                  LengthType const &r,
                                                  LengthType const &z, uint ir,
                                                  uint iz) const
  {
    double const fracZ = static_cast<double>((z - z_[iz]) * inverseDZ_);
    double const fracR = static_cast<double>((r - r_[ir]) * inverseDR_);
    double const invFracZ = 1.0 - fracZ;
    double const invFracR = 1.0 - fracR;

    return table[index(ir, iz)] * (invFracZ * invFracR) +
           table[index(ir, iz + 1)] * (fracZ * invFracR) +
           table[index(ir + 1, iz)] * (invFracZ * fracR) +
           table[index(ir + 1, iz + 1)] * (fracZ * fracR);
  }

  template <typename T>
  inline T RayTracingTable::Interpolate_(std::vector<T> const &table, LengthType const &r,
                                         LengthType const &z) const
  {
    auto const [ir, iz] = GetIrIz(r, z);
    return InterpolateFromIndex_(table, r, z, ir, iz);
  }

  inline double RayTracingTable::GetLaunch(LengthType const &r,
                                           LengthType const &z) const
  {
    return Interpolate_(launch_, r, z);
  }

  inline double RayTracingTable::GetReceive(LengthType const &r,
                                            LengthType const &z) const
  {
    return Interpolate_(receive_, r, z);
  }

  inline LengthType RayTracingTable::GetLength(LengthType const &r,
                                               LengthType const &z) const
  {
    return Interpolate_(length_, r, z);
  }

  inline TimeType RayTracingTable::GetDuration(LengthType const &r,
                                               LengthType const &z) const
  {
    return Interpolate_(duration_, r, z);
  }

  inline double RayTracingTable::GetFresnelS(LengthType const &r,
                                             LengthType const &z) const
  {
    return Interpolate_(fresnelS_, r, z);
  }

  inline double RayTracingTable::GetFresnelP(LengthType const &r,
                                             LengthType const &z) const
  {
    return Interpolate_(fresnelP_, r, z);
  }

  inline uint RayTracingTable::SigFigs_(double val, uint decimalPoints = 3) const
  {
    if (val == 0)
      return 1;
    if (abs(val) < 1)
    {
      if (val < 0)
      {
        return 3 + decimalPoints;
      }
      else
      {
        return 2 + decimalPoints;
      }
    }
    if (val < 0)
      val *= -10; // Flip sign, add one for neg sign

    double lgVal = log10(val);
    return std::max(1, int(lgVal + 1));
  }

  template <typename T1, typename T2>
  inline void RayTracingTable::PrintTable_(T1 const &table, T2 const unit) const
  {
    uint maxLen = 4;

    for (uint ir = 0; ir < r_.size(); ir++)
    {
      maxLen = std::max(maxLen, SigFigs_(r_[ir]));
    }

    for (uint iz = 0; iz < z_.size(); iz++)
    {
      maxLen = std::max(maxLen, SigFigs_(z_[iz]));
    }

    for (uint iz = 0; iz < nZBins_; iz++)
    {
      for (uint ir = 0; ir < nRBins_; ir++)
      {
        maxLen =
            std::max(maxLen, SigFigs_(table[index(ir, iz)] / unit));
      }
    }

    // Print the r labels
    std::cout << " dz ";
    for (uint ispace = 0; ispace < maxLen - 4; ispace++)
    {
      std::cout << ' ';
    }
    std::cout << "|| ";
    for (uint ir = 0; ir < r_.size(); ir++)
    {
      // Add leading spaces
      if (ir)
        std::cout << " | ";

      for (uint ispace = 0; ispace < maxLen - SigFigs_(r_[ir]); ispace++)
      {
        std::cout << " ";
      }
      std::cout << int(r_[ir]);
    }
    std::cout << " |\n";

    for (uint ispace = 0; ispace < maxLen; ispace++)
    {
      std::cout << '=';
    }
    std::cout << "||=";
    for (uint ir = 0; ir < r_.size(); ir++)
    {
      // Add leading spaces
      if (ir)
        std::cout << "=|=";

      for (uint ispace = 0; ispace < maxLen; ispace++)
      {
        std::cout << "=";
      }
    }
    std::cout << "=|\n";

    for (uint iz = 0; iz < table.size(); iz++)
    {

      for (uint ispace = 0; ispace < maxLen - SigFigs_(z_[iz]); ispace++)
      {
        std::cout << ' ';
      }

      std::cout << int(z_[iz]) << "|| ";
      for (uint ir = 0; ir < nRBins_; ir++)
      {
        if (ir)
          std::cout << " | ";

        for (uint ispace = 0; ispace < maxLen - SigFigs_(table[index(ir, iz)] / unit);
             ispace++)
        {
          std::cout << " ";
        }

        if (abs(table[index(ir, iz)] / unit) > 1)
        {
          std::cout << int(table[index(ir, iz)] / unit);
        }
        else
        {
          std::cout << int(table[index(ir, iz)] / unit * pow(10, 3)) / pow(10, 3);
        }
      }
      std::cout << " |\n";
    }
  }

  inline std::tuple<LengthType, LengthType> RayTracingTable::GetRAndZ_(
      Point const &x) const
  {
    auto const dr = x - plane_.getCenter();
    auto const z = dr.dot(plane_.getNormal());
    // r = dr - dr_z;  z = dr_z
    return std::make_tuple((dr - plane_.getNormal() * z).getNorm(), z);
  }

} // namespace c8_tracer
