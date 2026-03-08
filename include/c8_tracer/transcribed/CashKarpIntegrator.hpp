#pragma once

#include <array>

#include "c8_tracer/environment.hpp"
#include "c8_tracer/logger.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"
#include "c8_tracer/vec3.hpp"

namespace c8_tracer
{

  struct KahanSum
  {
    double sum = 0.0;
    double c = 0.0;
    inline void add(double x)
    {
      double y = x - c;
      double t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
  };

  class CashKarpIntegrator
  {
  private:
    double minStep_;
    double maxStep_;
    double inverseTol_;

    std::array<Vec3, 6> posPerStep_;
    std::array<Vec3, 6> dirPerStep_;

    // Cash-Karp coefficients
    // clang-format off
    double const parA_[36] = {            0.,          0.,            0.,               0.,           0., 0.,
                                     1. / 5.,          0.,            0.,               0.,           0., 0.,
                                    3. / 40.,    9. / 40.,            0.,               0.,           0., 0.,
                                    3. / 10.,   -9. / 10.,       6. / 5.,               0.,           0., 0., 
                                  -11. / 54.,     5. / 2.,    -70. / 27.,        35. / 27.,           0., 0.,
                              1631. / 55296., 175. / 512., 575. / 13824., 44275. / 110592., 253. / 4096., 0.
                            };
    // clang-format on

    double const parB_[6] = {37. / 378., 0, 250. / 621., 125. / 594., 0., 512. / 1771.};
    double const parC_[6] = {2825. / 27648., 0., 18575. / 48384., 13525. / 55296.,
                             277. / 14336., 1. / 4.};

  public:
    CashKarpIntegrator(double minStep, double maxStep, double tolerance)
        : minStep_(minStep), maxStep_(maxStep), inverseTol_(1.0 / tolerance)
    {
    }

    void AdaptiveStep(Vec3 const &startPos, Vec3 const &startDir, Vec3 &endPos,
                      Vec3 &endDir, double &h0, EnvironmentBase const &env,
                      LengthType &stepLength, double &avgN,
                      bool updateStep = true)
    {

      double h = std::max(minStep_, std::min(h0, maxStep_));
      double hNew = h;
      Vec3 posError;

      // Loop to take small enough steps to keep the error below tolerance
      while (true)
      {
        Step(startPos, startDir, endPos, endDir, posError, h, env, stepLength, avgN);

        auto const ratio = posError.norm() * inverseTol_;

        TRACER_LOG_ALL("Test step: pos " + endPos.to_string() + " dir " + endDir.to_string() +
                       " err " + std::to_string(posError.norm()) + " err-ratio " + std::to_string(ratio) +
                       " h " + std::to_string(h));

        double factor = std::pow(ratio, -0.2);
        hNew = h * 0.95 * factor; // update step size to keep error close to tolerance
        hNew = std::max(0.1 * h, std::min(hNew, 5.0 * h));
        hNew = std::max(minStep_, std::min(hNew, maxStep_));

        TRACER_LOG_ALL("h0: " + std::to_string(h0) + " h: " + std::to_string(h) +
                       " hNew: " + std::to_string(hNew));

        if (ratio <= 1.0)
          h = hNew;
        break;

        if (hNew == h)
          break;

        h = hNew;

        if (ratio <= 1.0)
          break;
      }
      if (updateStep)
      {
        h0 = 0.5 * (h0 + hNew); // Update the new step to keep close to tol in future steps
      }
    }

    void Step(Vec3 const &startPos, Vec3 const &startDir, Vec3 &endPos,
              Vec3 &endDir, Vec3 &posError, double const h,
              EnvironmentBase const &env,
              LengthType &stepLength, double &avgN)
    {

      // Copy the location
      endPos = Vec3(startPos);
      endDir = Vec3(startDir);
      posError = Vec3(0.0, 0.0, 0.0);

      KahanSum lenAcc;
      KahanSum nAcc;

      Vec3 tempPos(startPos);
      Vec3 tempDir(startDir);

      for (size_t i = 0; i < 6; i++)
      {

        tempPos = startPos;
        tempDir = startDir;

        for (size_t j = 0; j < i; j++)
        {
          double a = parA_[i * 6 + j] * h;
          tempPos = tempPos + posPerStep_[j] * a;
          tempDir = tempDir + dirPerStep_[j] * a;
        }

        // Calculate optical properties
        double const n_refract = env.get_n(tempPos);
        auto const inverseRefract = 1.0 / n_refract;

        posPerStep_[i] = tempDir.normalized() * inverseRefract; // \hat{v} / n
        dirPerStep_[i] = env.get_grad_n(tempPos) * inverseRefract *
                         inverseRefract; // Grad / n^2

        double const b = parB_[i] * h;
        double const bc = (parB_[i] - parC_[i]) * h;

        endPos = endPos + posPerStep_[i] * b;
        endDir = endDir + dirPerStep_[i] * b;
        posError = posError + posPerStep_[i] * bc;

        double const posNorm = posPerStep_[i].norm();
        lenAcc.add(posNorm * b);
        nAcc.add(n_refract * posNorm * b);
      }

      stepLength = lenAcc.sum;
      avgN = nAcc.sum / stepLength;

      endDir = endDir.normalized();
    }
  };

}
