namespace c8_tracer
{
  template <typename Func>
  double BrentMethod(
      Func func, double a, double b,
      LengthType faInit, LengthType fbInit,
      double tol = 1e-6, int max_iter = 100)
  {
    LengthType fa = faInit;
    LengthType fb = fbInit;

    if (fa * fb > 0.0)
    {
      TRACER_LOG_ERROR("Root not bracketed in initial call " +
                       std::to_string(fa) + " " + std::to_string(fb) +
                       " " + std::to_string(fa * fb));
      throw std::invalid_argument("Root not bracketed " +
                                  std::to_string(fa) + " " + std::to_string(fb));
    }

    // Ensure |fa| >= |fb|
    if (std::abs(fa) < std::abs(fb))
    {
      std::swap(a, b);
      std::swap(fa, fb);
    }

    double c = a;
    LengthType fc = fa;
    double prev_step = b - a;

    for (int iter = 0; iter < max_iter; ++iter)
    {
      if (fb == 0.0)
        return b;

      double s;

      // inverse quadratic interpolation, if possible
      if (fa != fc && fb != fc)
      {
        s = (a * fb * fc) / ((fa - fb) * (fa - fc)) +
            (b * fa * fc) / ((fb - fa) * (fb - fc)) +
            (c * fa * fb) / ((fc - fa) * (fc - fb));
      }
      else
      {
        // secant fallback
        s = b - fb * (b - a) / (fb - fa);
      }

      // conditions to reject interpolation and use bisection
      bool s_outside = (s <= std::min(a, b)) || (s >= std::max(a, b));
      bool too_far = std::abs(s - b) > std::abs(prev_step) / 2.0;
      bool small_step = std::abs(prev_step) < tol;
      bool f_increasing = std::abs(fb) > std::abs(fc);

      if (s_outside || too_far || small_step || f_increasing)
      {
        s = (a + b) * 0.5;
        prev_step = b - a;
      }
      else
      {
        prev_step = b - s;
      }

      LengthType fs = func(s);

      // Shift c → b, b → s
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

      // Ensure |fa| >= |fb|
      if (std::abs(fa) < std::abs(fb))
      {
        std::swap(a, b);
        std::swap(fa, fb);
      }

      if (std::abs(b - a) < tol)
        return b;
    }

    return b;
  }
}
