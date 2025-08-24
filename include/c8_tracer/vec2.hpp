#pragma once
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

namespace c8_tracer
{

  struct Vec2
  {
    double x, y;

    constexpr Vec2(double x = 0.0, double y = 0.0) : x(x), y(y) {}

    // Arithmetic
    constexpr Vec2 operator+(const Vec2 &rhs) const { return {x + rhs.x, y + rhs.y}; }
    constexpr Vec2 operator-(const Vec2 &rhs) const { return {x - rhs.x, y - rhs.y}; }
    constexpr Vec2 operator*(double scalar) const { return {x * scalar, y * scalar}; }
    constexpr Vec2 operator/(double scalar) const { return {x / scalar, y / scalar}; }

    // Dot product
    constexpr double dot(const Vec2 &rhs) const { return x * rhs.x + y * rhs.y; }

    // Magnitude
    double length() const { return std::sqrt(x * x + y * y); }
    double norm() const { return length(); }

    // Normalization
    Vec2 normalized() const
    {
      double len = length();
      return len > 0.0 ? (*this) / len : Vec2(0.0, 0.0);
    }

    // Convert to string
    std::string to_string() const
    {
      std::ostringstream oss;
      oss << "Vec2(" << x << ", " << y << ")";
      return oss.str();
    }
  };

  inline std::ostream &operator<<(std::ostream &os, const Vec2 &vec)
  {
    os << vec.to_string();
    return os;
  }
}
