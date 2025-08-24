#pragma once

#include "c8_tracer/logger.hpp"
#include "c8_tracer/vec3.hpp"

namespace c8_tracer
{
  // Base class for all environment classes
  class EnvironmentBase
  {
  public:
    EnvironmentBase() = default;
    virtual ~EnvironmentBase() = default;

    // Returns a scalar field value at a given position
    virtual double get_n(const Vec3 &position) const = 0;

    // Returns the gradient of the scalar field at a given position
    virtual Vec3 get_grad_n(const Vec3 &position) const = 0;
  };

  /**
   * Construct a IsotropicEnvironment.
   * Implements an environment where the gradient is zero everywhere
   *
   *
   * @param nRefrac    index of refraction at all points in space
   */
  class IsotropicEnvironment : public EnvironmentBase
  {
  private:
    double nRefrac_;

  public:
    IsotropicEnvironment(double nRefrac) : nRefrac_(nRefrac) {}

    double get_n(const Vec3 &) const override
    {
      return nRefrac_;
    }

    Vec3 get_grad_n(const Vec3 &) const override
    {
      return Vec3(0.0, 0.0, 0.0);
    }
  };

  /**
   * Construct a LinearRadialEnvironment.
   * Implements an environment where the gradient is constant and
   * points along r_hat
   *
   *
   * @param center     center of the spherical field
   * @param nCenter    index of refraction at `center`
   * @param dNdR       the rate of change of the index of refraction
   */
  class LinearRadialEnvironment : public EnvironmentBase
  {
  private:
    Vec3 center_;
    double nCenter_;
    double dNdR_;

  public:
    LinearRadialEnvironment(Vec3 const &center, double nCenter, double dNdR)
        : center_(center), nCenter_(nCenter), dNdR_(dNdR) {}

    double get_n(const Vec3 &position) const override
    {
      return (position - center_).norm() * dNdR_ + nCenter_;
    }

    Vec3 get_grad_n(const Vec3 &position) const override
    {
      return (position - center_).normalized() * dNdR_;
    }
  };

  /**
   * Construct a CartesianLinearEnvironment.
   * Implements an environment where the gradient has a constant cartesian
   * vector. i.e. you can test intoductory kinematics with this
   *
   *
   * @param gradVec     gradient of n
   * @param refPoint    defines where `nRef` is known
   * @param nRef        index of refraction at `refPoint`
   */
  class CartesianLinearEnvironment : public EnvironmentBase
  {
  private:
    Vec3 gradVec_;
    Vec3 refPoint_;
    double nRef_;

  public:
    CartesianLinearEnvironment(Vec3 const &gradVec, Vec3 const &refPoint, double nRef)
        : gradVec_(gradVec), refPoint_(refPoint), nRef_(nRef) {}

    double get_n(const Vec3 &position) const override
    {
      return gradVec_.dot(position - refPoint_) + nRef_;
    }

    Vec3 get_grad_n(const Vec3 &) const override
    {
      return gradVec_;
    }
  };

  /**
   * Construct a CartesianSingleExponentialEnvironment.
   * Implements an environment where the gradient is exponential along
   * a single, cartesian direciton
   *
   *
   * @param nAsymptotic     asymtotic value of n
   * @param deltaN          amplitude of the exponential dependence
   * @param lambda          length-scale of the exponential
   * @param gradientAxis    axis along which the gradient changes points towards
   *                        decreasing values of n (typically towards the zenith)
   * @param refPoint        reference point for the exponential
   */
  class CartesianSingleExponentialEnvironment : public EnvironmentBase
  {
  private:
    double nAsymptotic_;
    double deltaN_;
    double inverseLambda_;
    Vec3 gradDir_;
    Vec3 refPoint_;

  public:
    CartesianSingleExponentialEnvironment(double nAsymptotic, double deltaN,
                                          double lambda, Vec3 const &gradientAxis,
                                          Vec3 const &refPoint)
        : nAsymptotic_(nAsymptotic),
          deltaN_(deltaN),
          inverseLambda_(1.0 / lambda),
          gradDir_(gradientAxis.normalized()),
          refPoint_(refPoint) {}

    double get_n(const Vec3 &position) const override
    {
      // n0 - dn exp((z - z0) / lambda)
      auto const z = (position - refPoint_).dot(gradDir_);
      return nAsymptotic_ - deltaN_ * exp(z * inverseLambda_);
    }

    Vec3 get_grad_n(const Vec3 &position) const override
    {
      // - dn / lambda exp((z - z0) / lambda)
      auto const z = (position - refPoint_).dot(gradDir_);
      return gradDir_ * -deltaN_ * inverseLambda_ * exp(z * inverseLambda_);
    }
  };

}
