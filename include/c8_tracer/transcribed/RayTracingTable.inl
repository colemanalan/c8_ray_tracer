#pragma once

#include <iomanip>
#include <iostream>
#include "c8_tracer/transcribed/c8_typedefs.hpp"
#include "c8_tracer/logger.hpp"

namespace c8_tracer
{

  inline RayTracingTable::RayTracingTable(LengthType minR, LengthType maxR, uint nRBins,
                                          LengthType minZ, LengthType maxZ, uint nZBins,
                                          Plane const &plane, double nCenter)
      : nRBins_(nRBins), nZBins_(nZBins), maxR_(maxR), minR_(minR),
        maxZ_(maxZ), minZ_(minZ), r_(nRBins), z_(nZBins), plane_(plane), indexOfRefraction_(nCenter)
  {

    if (nRBins_ < 2)
    {
      LOG_ERROR("Initializing with radius bins n = " + std::to_string(nRBins_));
      throw std::invalid_argument("Must have >1 bins");
    }
    if (nZBins_ < 2)
    {
      LOG_ERROR("Initializing with z bins n = " + std::to_string(nZBins_));
      throw std::invalid_argument("Must have >1 bins");
    }

    if (maxZ_ <= minZ_)
    {
      LOG_ERROR("Max z " + std::to_string(maxZ_) + ", Min z " + std::to_string(minZ_));
      throw std::invalid_argument("Required: max > min");
    }

    if (maxR_ <= minR_)
    {
      LOG_ERROR("Max r " + std::to_string(maxR_) + ", Min r " + std::to_string(minR_));
      throw std::invalid_argument("Required: max > min");
    }

    if (minR_ < 0.0)
    {
      LOG_ERROR("Min r " + std::to_string(minR_) + " must be > 0");
      throw std::invalid_argument("Min r must be > 0");
    }

    // construct the tables here
    size_t N = nRBins_ * nZBins_;
    launch_ = std::make_unique<double[]>(N);
    receive_ = std::make_unique<double[]>(N);
    length_ = std::make_unique<LengthType[]>(N);
    duration_ = std::make_unique<TimeType[]>(N);
    fresnelS_ = std::make_unique<double[]>(N);
    fresnelP_ = std::make_unique<double[]>(N);

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
      LOG_INFO("There were " + std::to_string(percentUntracked) + " % of tracks that were not contained within this RayTracingTable");
      LOG_INFO("Total tracks: " + std::to_string(totalCalls_) + ", not contained tracks: " + std::to_string(untrackedCalls_));
      if (percentUntracked > 0.2)
      {
        LOG_WARNING(
            "This percentage is large and may indicate that the table"
            "is not large enough to contain your cascade");
      }
    }
  }

  inline void RayTracingTable::ResetTables()
  {
    size_t N = nRBins_ * nZBins_;

    FillWithNaN_(launch_, N);
    FillWithNaN_(receive_, N);
    FillWithNaN_(length_, N);
    FillWithNaN_(duration_, N);
    FillWithNaN_(fresnelS_, N);
    FillWithNaN_(fresnelP_, N);

    totalCalls_ = 0;
    untrackedCalls_ = 0;
  }

  ///////////////////////
  /// Setters
  ///////////////////////

  inline void RayTracingTable::SetLaunch(double launch, uint ir, uint iz)
  {
    if (!IsValid(ir, iz, true))
    {
      LOG_ERROR("Index out of range");
    }
    launch_[index(ir, iz)] = launch;
  }

  inline void RayTracingTable::SetReceive(double receive, uint ir, uint iz)
  {
    if (!IsValid(ir, iz, true))
    {
      LOG_ERROR("Index out of range");
    }
    receive_[index(ir, iz)] = receive;
  }

  inline void RayTracingTable::SetLength(LengthType length, uint ir, uint iz)
  {
    if (!IsValid(ir, iz, true))
    {
      LOG_ERROR("Index out of range");
    }
    length_[index(ir, iz)] = length;
  }

  inline void RayTracingTable::SetDuration(TimeType duration, uint ir, uint iz)
  {
    if (!IsValid(ir, iz, true))
    {
      LOG_ERROR("Index out of range");
    }
    duration_[index(ir, iz)] = duration;
  }

  inline void RayTracingTable::SetFresnelS(double fresnelS, uint ir, uint iz)
  {
    if (!IsValid(ir, iz, true))
    {
      LOG_ERROR("Index out of range");
    }
    fresnelS_[index(ir, iz)] = fresnelS;
  }

  inline void RayTracingTable::SetFresnelP(double fresnelP, uint ir, uint iz)
  {
    if (!IsValid(ir, iz, true))
    {
      LOG_ERROR("Index out of range");
    }
    fresnelP_[index(ir, iz)] = fresnelP;
  }

  ///////////////////////
  /// Getters
  ///////////////////////

  inline double RayTracingTable::GetLaunch(uint ir, uint iz) const
  {
    if (!IsValid(ir, iz, true))
    {
      LOG_ERROR("Index out of range");
    }
    return launch_[index(ir, iz)];
  }

  inline double RayTracingTable::GetReceive(uint ir, uint iz) const
  {
    if (!IsValid(ir, iz, true))
    {
      LOG_ERROR("Index out of range");
    }
    return receive_[index(ir, iz)];
  }

  inline LengthType RayTracingTable::GetLength(uint ir, uint iz) const
  {
    if (!IsValid(ir, iz, true))
    {
      LOG_ERROR("Index out of range");
    }
    return length_[index(ir, iz)];
  }

  inline TimeType RayTracingTable::GetDuration(uint ir, uint iz) const
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      LOG_ERROR("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    return duration_[index(ir, iz)];
  }

  inline double RayTracingTable::GetFresnelS(uint ir, uint iz) const
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      LOG_ERROR("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
    return fresnelS_[index(ir, iz)];
  }

  inline double RayTracingTable::GetFresnelP(uint ir, uint iz) const
  {
    if (ir >= nRBins_ || iz >= nZBins_)
      LOG_ERROR("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
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

  inline SignalPath RayTracingTable::GetSignalPath(Point const &x0, double n0) const
  {
    totalCalls_++;
    // Do not extrapolate
    if (!ContainsPoint(x0))
    {
      untrackedCalls_++;
      return SignalPath(std::numeric_limits<double>::infinity(), 1.0, n0, indexOfRefraction_,
                        plane_.getNormal(), plane_.getNormal(),
                        std::numeric_limits<double>::infinity(), Path(x0), 1.0,
                        1.0);
    }

    auto const [r, z] = GetRAndZ_(x0);
    auto const [ir, iz] = GetIrIz(r, z);
    auto const length = InterpolateFromIndex_(length_.get(), r, z, ir, iz);
    if (std::isnan(length))
    {
      return SignalPath(std::numeric_limits<double>::infinity(), 1.0, n0, indexOfRefraction_,
                        plane_.getNormal(), plane_.getNormal(),
                        std::numeric_limits<double>::infinity(), Path(x0), 1.0,
                        1.0);
    }
    auto const duration = InterpolateFromIndex_(duration_.get(), r, z, ir, iz);
    auto const fresnelS = InterpolateFromIndex_(fresnelS_.get(), r, z, ir, iz);
    auto const fresnelP = InterpolateFromIndex_(fresnelP_.get(), r, z, ir, iz);

    auto const launch = InterpolateFromIndex_(launch_.get(), r, z, ir, iz);
    auto const receive = InterpolateFromIndex_(receive_.get(), r, z, ir, iz);

    auto const dr = (plane_.getCenter() - x0);
    auto const xyComponents =
        (dr - plane_.getNormal() * dr.dot(plane_.getNormal())).normalized();

    auto const startDir =
        plane_.getNormal() * launch + xyComponents * sqrt(1 - launch * launch);
    auto const endDir =
        plane_.getNormal() * receive - xyComponents * sqrt(1 - receive * receive);

    return SignalPath(duration,
                      duration * constants::c / length, // average n
                      n0, indexOfRefraction_, startDir, endDir, length, Path(x0), fresnelS, fresnelP);
  }

  inline bool RayTracingTable::IsValid(uint ir, uint iz, bool warn = false) const
  {
    if (ir >= nRBins_ || iz >= nZBins_)
    {
      if (warn)
      {
        LOG_ERROR("Index out of range R " + std::to_string(ir) + "/" + std::to_string(nRBins_) +
                  " Z " + std::to_string(iz) + "/" + std::to_string(nZBins_));
      }
      return false;
    }
    return true;
  }

  inline void RayTracingTable::Print() const
  {
    std::cout << "\nTables for observer at " << plane_.getCenter();
    std::cout << "\n";
    std::cout << "Path lengths / m\n";
    PrintTable_(length_, 1.0);
    std::cout << "\n";
    std::cout << "Duration / ns\n";
    PrintTable_(duration_, 1e-9);
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
  inline T RayTracingTable::InterpolateFromIndex_(
      T const *table,
      LengthType const &r,
      LengthType const &z,
      uint ir,
      uint iz) const noexcept
  {
    double fracZ = (z - z_[iz]) * inverseDZ_;
    double fracR = (r - r_[ir]) * inverseDR_;

    uint base = iz * nRBins_ + ir;

    T const v00 = table[base];
    T const v01 = table[base + 1];
    T const v10 = table[base + nRBins_];
    T const v11 = table[base + nRBins_ + 1];

    if (std::isnan(v00) || std::isnan(v01) ||
        std::isnan(v10) || std::isnan(v11))
    {
      return std::numeric_limits<T>::quiet_NaN();
    }

    T const top = v00 + fracR * (v01 - v00);
    T const bottom = v10 + fracR * (v11 - v10);
    return top + fracZ * (bottom - top);
  }

  template <typename T>
  inline T RayTracingTable::Interpolate_(T const *table, LengthType const &r,
                                         LengthType const &z) const noexcept
  {
    auto const [ir, iz] = GetIrIz(r, z);
    return InterpolateFromIndex_(table, r, z, ir, iz);
  }

  inline double RayTracingTable::GetLaunch(LengthType const &r,
                                           LengthType const &z) const
  {
    return Interpolate_(launch_.get(), r, z);
  }

  inline double RayTracingTable::GetReceive(LengthType const &r,
                                            LengthType const &z) const
  {
    return Interpolate_(receive_.get(), r, z);
  }

  inline LengthType RayTracingTable::GetLength(LengthType const &r,
                                               LengthType const &z) const
  {
    return Interpolate_(length_.get(), r, z);
  }

  inline TimeType RayTracingTable::GetDuration(LengthType const &r,
                                               LengthType const &z) const
  {
    return Interpolate_(duration_.get(), r, z);
  }

  inline double RayTracingTable::GetFresnelS(LengthType const &r,
                                             LengthType const &z) const
  {
    return Interpolate_(fresnelS_.get(), r, z);
  }

  inline double RayTracingTable::GetFresnelP(LengthType const &r,
                                             LengthType const &z) const
  {
    return Interpolate_(fresnelP_.get(), r, z);
  }

  inline uint RayTracingTable::SigFigs_(double val, uint decimalPoints = 3) const
  {
    if (std::isnan(val))
    {
      return 3;
    }
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
    using std::abs;
    using std::setfill;
    using std::setw;

    //////////////////////////////
    /////// compute sig figs /////
    //////////////////////////////
    uint maxLen = 4;

    std::vector<uint> rWidth(r_.size());
    for (uint ir = 0; ir < r_.size(); ir++)
    {
      rWidth[ir] = SigFigs_(r_[ir]);
      maxLen = std::max(maxLen, rWidth[ir]);
    }

    std::vector<uint> zWidth(z_.size());
    for (uint iz = 0; iz < z_.size(); iz++)
    {
      zWidth[iz] = SigFigs_(z_[iz]);
      maxLen = std::max(maxLen, zWidth[iz]);
    }

    std::vector<uint> cellWidth(nRBins_ * nZBins_);
    for (uint iz = 0; iz < nZBins_; iz++)
    {
      for (uint ir = 0; ir < nRBins_; ir++)
      {
        double val = table[index(ir, iz)] / unit;
        uint w = SigFigs_(val);
        cellWidth[index(ir, iz)] = w;
        maxLen = std::max(maxLen, w);
      }
    }

    //////////////////////////////
    ///////// print header ///////
    //////////////////////////////
    std::cout << " dz " << std::setw(maxLen - 4) << "" << "|| ";

    for (uint ir = 0; ir < r_.size(); ir++)
    {
      if (ir)
        std::cout << " | ";
      std::cout << std::setw(maxLen) << int(r_[ir]);
    }
    std::cout << " |\n";

    //////////////////////////////
    ////// print separator ///////
    //////////////////////////////
    std::cout << std::setw(maxLen) << std::setfill('=') << "" << setfill(' ')
              << "||=";

    for (uint ir = 0; ir < r_.size(); ir++)
    {
      if (ir)
        std::cout << "=|=";
      std::cout << std::setw(maxLen) << std::setfill('=') << "" << setfill(' ');
    }
    std::cout << "=|\n";

    //////////////////////////////
    ////// print table rows //////
    //////////////////////////////
    for (uint iz = 0; iz < z_.size(); iz++)
    {

      std::cout << std::setw(maxLen) << int(z_[iz]) << "|| ";

      for (uint ir = 0; ir < nRBins_; ir++)
      {
        if (ir)
          std::cout << " | ";

        double val = table[index(ir, iz)] / unit;

        if (std::isnan(val))
        {
          std::cout << setw(maxLen) << "NaN";
          continue;
        }

        // Round small values to 3 decimals
        if (abs(val) > 1.0)
        {
          std::cout << std::setw(maxLen) << std::lround(val);
        }
        else
        {
          double rounded = std::round(val * 1000.0) / 1000.0;
          std::cout << std::setw(maxLen) << rounded;
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

  template <typename T>
  inline void RayTracingTable::FillWithNaN_(std::unique_ptr<T[]> &ptr, size_t N)
  {
    static_assert(std::is_floating_point_v<T>,
                  "FillWithNaN_ requires a floating‑point type");

    T nan = std::numeric_limits<T>::quiet_NaN();
    std::fill_n(ptr.get(), N, nan);
  }

} // namespace c8_tracer
