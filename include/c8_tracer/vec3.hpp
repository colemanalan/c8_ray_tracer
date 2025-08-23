#pragma once
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

namespace c8_tracer
{

  struct Vec3
  {
    double x, y, z;

    constexpr Vec3(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) {}

    // Arithmetic
    constexpr Vec3 operator+(const Vec3 &rhs) const { return {x + rhs.x, y + rhs.y, z + rhs.z}; }
    constexpr Vec3 operator-(const Vec3 &rhs) const { return {x - rhs.x, y - rhs.y, z - rhs.z}; }
    constexpr Vec3 operator*(double scalar) const { return {x * scalar, y * scalar, z * scalar}; }
    constexpr Vec3 operator/(double scalar) const { return {x / scalar, y / scalar, z / scalar}; }

    // Dot product
    constexpr double dot(const Vec3 &rhs) const { return x * rhs.x + y * rhs.y + z * rhs.z; }

    // Cross product
    constexpr Vec3 cross(const Vec3 &rhs) const
    {
      return {
          y * rhs.z - z * rhs.y,
          z * rhs.x - x * rhs.z,
          x * rhs.y - y * rhs.x};
    }

    // Magnitude
    double length() const { return std::sqrt(x * x + y * y + z * z); }
    double norm() const { return length(); }

    // Functions to match C8 calls
    double getNorm() const { return length(); } // to match C8

    // Normalization
    Vec3 normalized() const
    {
      double len = length();
      return len > 0.0 ? (*this) / len : Vec3(0.0, 0.0, 0.0);
    }

    // Convert to string
    std::string to_string() const
    {
      std::ostringstream oss;
      oss << "Vec3(" << x << ", " << y << ", " << z << ")";
      return oss.str();
    }
  };
}

// Enable std::to_string(Vec3)
namespace std
{
  inline std::string to_string(const c8_tracer::Vec3 &v)
  {
    return v.to_string();
  }
}