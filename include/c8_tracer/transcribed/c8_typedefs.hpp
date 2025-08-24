#pragma once

#include "c8_tracer/vec3.hpp"

namespace c8_tracer
{
  typedef Vec3 Point;
  typedef Vec3 DirectionVector;

  typedef double LengthType;
  typedef double InverseLengthType;
  typedef double LengthTypeSq;
  typedef double TimeType;

} // namespace c8_tracer

namespace constants
{
  constexpr double c = 299792458.0;
}
