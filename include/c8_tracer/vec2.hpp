#pragma once
#include <cmath>

namespace c8_tracer
{

    struct Vec2
    {
        float x, y;

        constexpr Vec2(float x = 0.0f, float y = 0.0f) : x(x), y(y) {}

        // Arithmetic
        constexpr Vec2 operator+(const Vec2 &rhs) const { return {x + rhs.x, y + rhs.y}; }
        constexpr Vec2 operator-(const Vec2 &rhs) const { return {x - rhs.x, y - rhs.y}; }
        constexpr Vec2 operator*(float scalar) const { return {x * scalar, y * scalar}; }
        constexpr Vec2 operator/(float scalar) const { return {x / scalar, y / scalar}; }

        // Dot product
        constexpr float dot(const Vec2 &rhs) const { return x * rhs.x + y * rhs.y; }

        // Magnitude
        float length() const { return std::sqrt(x * x + y * y); }

        // Normalization
        Vec2 normalized() const
        {
            float len = length();
            return len > 0.0f ? (*this) / len : Vec2(0.0f, 0.0f);
        }
    };
}
