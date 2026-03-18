#pragma once

#include <iostream>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include "c8_tracer/transcribed/c8_typedefs.hpp"
#include "c8_tracer/logger.hpp"

namespace c8_tracer
{
  /**
   * Perform the Brent method for root-finding
   *
   * @param func function that takes a double and returns a signed function value.
   * This function defines the root-to-be-found, `f(x)`, must return a `LengthType`
   * @param a bounding value of the function argument corresponding to `faInit`
   * @param b bounding value of the function argument corresponding to `fbInit`
   * @param faInit function value for argument `x=a`
   * @param fbInit function value for argument `x=b`
   * @param tol will stop when the search bounds are within `tol` of each other
   * @param max_iter will fail if this many steps have been taken, regardless of found root
   *
   * @return the x-value corresponding to the root of `func`
   */
  inline double BrentMethod(
      std::function<LengthType(double)> func, double a, double b,
      LengthType faInit, LengthType fbInit, double tol = 1e-6,
      int max_iter = 100)
  {
    LengthType fa = faInit, fb = fbInit;
    if (fa * fb > 0.0)
    {
      TRACER_LOG_ERROR("Root not bracketed in initial call " +
                       std::to_string(fa) + " " + std::to_string(fb) +
                       " " + std::to_string(fa * fb));
      throw std::invalid_argument("Root not bracketed " +
                                  std::to_string(fa) + " " + std::to_string(fb));
    }

    if (std::abs(fa) < std::abs(fb))
    {
      std::swap(a, b);
      std::swap(fa, fb);
    }

    LengthType fc = fa;
    double c = a;
    double e = b - a;

    for (int iter = 0; iter < max_iter; ++iter)
    {
      if (fb == 0.0)
        return b;

      double s;
      if (fa != fc && fb != fc)
      {
        // Inverse quadratic interpolation
        s = (a * fb * fc) / ((fa - fb) * (fa - fc)) + (b * fa * fc) / ((fb - fa) * (fb - fc)) + (c * fa * fb) / ((fc - fa) * (fc - fb));
      }
      else
      {
        // Secant method
        s = b - fb * (b - a) / (fb - fa);
      }

      // Conditions to accept s or fall back to bisection
      bool cond1 = !(std::min(a, b) < s && s < std::max(a, b));
      bool cond2 = std::abs(s - b) > std::abs(e) / 2.0;
      bool cond3 = std::abs(e) < tol;
      bool cond4 = std::abs(fb) > std::abs(fc);

      if (cond1 || cond2 || cond3 || cond4)
      {
        s = (a + b) / 2.0;
        e = b - a;
      }
      else
      {
        e = b - s;
      }

      LengthType fs = func(s);
      c = b;
      fc = fb;

      if (fa * fs < 0.0)
      {
        b = s;
        fb = fs;
      }
      else
      {
        a = s;
        fa = fs;
      }

      if (std::abs(fa) < std::abs(fb))
      {
        std::swap(a, b);
        std::swap(fa, fb);
      }

      if (std::abs(b - a) < tol)
      {
        return b;
      }
    }

    return b;
  }
}
