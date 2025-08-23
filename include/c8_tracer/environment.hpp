#pragma once

#include "c8_tracer/vec3.hpp"

namespace c8_tracer
{
    class environment_base
    {
    public:
        environment_base() = default;
        virtual ~environment_base() = default;

        // Returns a scalar field value at a given position
        virtual float get_n(const Vec3 &position) const = 0;

        // Returns the gradient of the scalar field at a given position
        virtual Vec3 get_n_grad(const Vec3 &position) const = 0;
    };

    class isotropic_environment : public environment_base
    {
    private:
        float n_refrac_;

    public:
        isotropic_environment(float n_refrac) : n_refrac_(n_refrac) {}

        float get_n(const Vec3 &position) const override
        {
            return n_refrac_;
        }

        Vec3 get_n_grad(const Vec3 &position) const override
        {
            return Vec3(0.0f, 0.0f, 0.0f);
        }
    };

}
