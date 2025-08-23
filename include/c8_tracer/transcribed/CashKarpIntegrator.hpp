#pragma once

#include <array>

#include "c8_tracer/environment.hpp"
#include "c8_tracer/logger.hpp"
#include "c8_tracer/vec3.hpp"

namespace c8_tracer
{
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
                      bool updateStep = true)
    {

      double h = std::max(minStep_, std::min(h0, maxStep_));
      double hNew = h;
      double ratio = 1;

      auto x0(startPos);
      auto v0(startDir);
      auto xf(startPos);
      auto vf(startDir);
      auto error(startDir);

      // Loop to take small enough steps to keep the error below tolerance
      do
      {
        h = hNew;
        Step(x0, v0, xf, vf, error, h, env);

        ratio = error.norm() * inverseTol_;
        hNew = h * 0.95f *
               pow(ratio, -0.2f); // update step size to keep error close to tolerance
        hNew = std::max(0.1f * h, std::min(hNew, 5 * h));
        hNew = std::max(minStep_, std::min(hNew, maxStep_));

        logger.debug("Test step: pos " + xf.to_string() + " dir " + vf.to_string() +
                     " err " + std::to_string(error.norm()) + " err-ratio " + std::to_string(ratio) +
                     " h " + std::to_string(h) + " hNew " + std::to_string(hNew));

        if (hNew == minStep_)
        {
          logger.debug("Hist min step size " + std::to_string(minStep_));
          break; // performed step already at the minimum
        }
      } while (ratio > 1);

      if (updateStep)
      {
        h0 =
            0.5 * (h0 + hNew); // Update the new step to keep close to tol in future steps
      }
      endPos = xf;
      endDir = vf;
    }

    void Step(Vec3 const &startPos, Vec3 const &startDir, Vec3 &endPos,
              Vec3 &endDir, Vec3 &dirError, double const h,
              EnvironmentBase const &env)
    {

      // Copy the location
      endPos = Vec3(startPos);
      endDir = Vec3(startDir);
      dirError = Vec3(0.0, 0.0, 0.0);

      for (size_t i = 0; i < 6; i++)
      {

        Vec3 tempPos(startPos);
        Vec3 tempDir(startDir);

        for (size_t j = 0; j < i; j++)
        {
          tempPos = tempPos + posPerStep_[j] * parA_[i * 6 + j] * h;
          tempDir = tempDir + dirPerStep_[j] * parA_[i * 6 + j] * h;
        }

        // Calculate optical properties
        auto const inverseRefract = 1.0 / env.get_n(tempPos);

        posPerStep_[i] = tempDir.normalized() * inverseRefract; // \hat{v} / n
        dirPerStep_[i] = env.get_grad_n(tempPos) * inverseRefract *
                         inverseRefract; // Grad / n^2

        endPos = endPos + posPerStep_[i] * parB_[i] * h;
        endDir = endDir + dirPerStep_[i] * parB_[i] * h;
        dirError = dirError + dirPerStep_[i] * (parB_[i] - parC_[i]) * h;
      }
    }
  };

}