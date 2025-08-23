#pragma once

#include "c8_tracer/vec3.hpp"

namespace c8_tracer
{
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

  class IsotropicEnvironment : public EnvironmentBase
  {
  private:
    double n_refrac_;

  public:
    IsotropicEnvironment(double n_refrac) : n_refrac_(n_refrac) {}

    double get_n(const Vec3 &position) const override
    {
      return n_refrac_;
    }

    Vec3 get_grad_n(const Vec3 &position) const override
    {
      return Vec3(0.0, 0.0, 0.0);
    }
  };

  class LinearRadialEnvironment : public EnvironmentBase
  {
  private:
    Vec3 center_;
    double nCenter_;
    double dNdR_;

  public:
    LinearRadialEnvironment(Vec3 const &center, double nCenter, double dNdR)
        : center_(center), nCenter_(nCenter) {}

    double get_n(const Vec3 &position) const override
    {
      return (position - center_).norm() * dNdR_ + nCenter_;
    }

    Vec3 get_grad_n(const Vec3 &position) const override
    {
      return (position - center_).normalized() * dNdR_;
    }
  };

  class CartesianLinearEnvironment : public EnvironmentBase
  {
  private:
    Vec3 gradVec_;  // gradient of n
    Vec3 refPoint_; // defines where `nRef_` is known
    double nRef_;   // index of refraction at `refPoint`

  public:
    CartesianLinearEnvironment(Vec3 const &gradVec, Vec3 const &refPoint, double nRef)
        : gradVec_(gradVec), refPoint_(refPoint), nRef_(nRef) {}

    double get_n(const Vec3 &position) const override
    {
      return gradVec_.dot(position - refPoint_) + nRef_;
    }

    Vec3 get_grad_n(const Vec3 &position) const override
    {
      return gradVec_;
    }
  };

}
