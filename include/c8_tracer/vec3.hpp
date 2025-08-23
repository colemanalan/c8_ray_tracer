#pragma once
#include <cmath>

namespace c8_tracer
{

    struct Vec3
    {
        float x, y, z;

        constexpr Vec3(float x = 0.0f, float y = 0.0f, float z = 0.0f) : x(x), y(y), z(z) {}

        // Arithmetic
        constexpr Vec3 operator+(const Vec3 &rhs) const { return {x + rhs.x, y + rhs.y, z + rhs.z}; }
        constexpr Vec3 operator-(const Vec3 &rhs) const { return {x - rhs.x, y - rhs.y, z - rhs.z}; }
        constexpr Vec3 operator*(float scalar) const { return {x * scalar, y * scalar, z * scalar}; }
        constexpr Vec3 operator/(float scalar) const { return {x / scalar, y / scalar, z / scalar}; }

        // Dot product
        constexpr float dot(const Vec3 &rhs) const { return x * rhs.x + y * rhs.y + z * rhs.z; }

        // Cross product
        constexpr Vec3 cross(const Vec3 &rhs) const
        {
            return {
                y * rhs.z - z * rhs.y,
                z * rhs.x - x * rhs.z,
                x * rhs.y - y * rhs.x};
        }

        // Magnitude
        float length() const { return std::sqrt(x * x + y * y + z * z); }

        // Normalization
        Vec3 normalized() const
        {
            float len = length();
            return len > 0.0f ? (*this) / len : Vec3(0.0f, 0.0f, 0.0f);
        }
    };
}