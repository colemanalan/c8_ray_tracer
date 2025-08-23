#pragma once

#include "c8_tracer/logger.hpp"

namespace c8_tracer
{
    inline int add(int a, int b)
    {
        logger.info("After addition the result is " + std::to_string(a + b));
        return a + b;
    }
}