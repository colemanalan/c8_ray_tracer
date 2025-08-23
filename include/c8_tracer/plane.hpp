#pragma once

#include "c8_tracer/vec3.hpp"

namespace c8_tracer
{
  class Plane
  {
  private:
    Vec3 center_;
    Vec3 normal_;

  public:
    Plane(Vec3 const &center, Vec3 const &normal)
        : center_(center), normal_(normal) {}

    Vec3 getNormal() const { return normal_; }
    Vec3 getCenter() const { return center_; }
  };

}