#pragma once

#include "c8_tracer/logger.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"
#include "c8_tracer/transcribed/brent.hpp"
#include "c8_tracer/transcribed/RayTracer.hpp"

#define STOP_CLOSE_LENGTH 0.00001 // max distance for stopping criteria in Brent Loops
#define PLANE_CLOSE_TOL 0.000005  // max distance for finding intersections with planes
#define DCOS_TOL 1e-6             // limit in cosine distance before quitting opt loops
#define MAX_RAY_STEPS 10000       // max steps taken before quitting
#define INITAL_STEP_SIZE 0.1      // fist step will always been this large

namespace c8_tracer
{

#define str std::to_string

  inline RayTracer2D::RayTracer2D(DirectionVector const axis,
                                  LengthType const minStep = 0.0001,
                                  LengthType const maxStep = 10.0,
                                  double const tolerance = 1e-8,
                                  int const brentRays = 13)
      : RayTracerBase(axis), tracer_(minStep, maxStep, tolerance), nRays_(brentRays)
  {
    if (brentRays < 3)
    {
      TRACER_LOG_ERROR("Must use at least 3 Brent rays, gave" + str(brentRays));
      throw std::runtime_error("Bad construction of RayTracer2D");
    }
  }

  inline void RayTracer2D::AddReflectionLayer(Plane const layer)
  {
    // Check to ensure that this remains a 2D problem
    if (std::abs(layer.getNormal().dot(axis_)) < cos(0.001 * M_PI / 180.))
    {
      TRACER_LOG_ERROR(
          "Cannot add a reflection layer that is not perpendicular"
          "the axis. Layer normal " +
          str(layer.getNormal()) + ", ray tracer symmetry axis " + str(axis_));
      throw std::invalid_argument("Bad reflection layer");
    }

    reflectionLayers_.push_back(layer);
  }

  inline std::vector<SignalPath> RayTracer2D::GetSignalPaths(Point const &sourceXX,
                                                             Point const &targetXX,
                                                             EnvironmentBase const &env)
  {
    std::vector<SignalPath> signalPaths;
    signalPaths.reserve(2);

    bool flippedStartEnd = false;
    Point start = sourceXX;
    Point end = targetXX;

    if (_GetZ(end) > _GetZ(start))
    {
      TRACER_LOG_DEBUG("Start " + str(start) + " is below end " + str(end) + ", flipping");
      std::swap(start, end);
      flippedStartEnd = true;
    }

    // 4 iterations of 11 rays will find solutions that are at least separated by 0.000855
    auto const [foundEmit, foundReceive] = FindEmitAndReceiveBrent(start, end, env, nRays_, 4);
    TRACER_LOG_INFO("Creating paths for " + str(foundEmit.size()) + " solutions");

    if (foundEmit.size())
    {
      for (uint i = 0; i < foundEmit.size(); i++)
      {
        TRACER_LOG_DEBUG("Gettting path for x0 " + str(start) + " v0 " + str(foundEmit[i]));
        auto const path = GetSignalPath(start, foundEmit[i], end, env);
        signalPaths.push_back(flippedStartEnd ? FlipSignalPath(path) : path);
      }

      // sort results by propagation time
      std::sort(signalPaths.begin(), signalPaths.end(),
                [](auto const &a, auto const &b)
                {
                  return std::abs(a.propagation_time_) < std::abs(b.propagation_time_);
                });
    }

    return signalPaths;
  }

  inline std::vector<SignalPath> RayTracer2D::PropagateToPoint(Point const &sourceXX,
                                                               Point const &targetXX,
                                                               EnvironmentBase const &env)
  {

    std::vector<SignalPath> signalPaths;
    signalPaths.reserve(2);

    bool flippedStartEnd = false;
    Point start = sourceXX;
    Point end = targetXX;

    if (_GetZ(end) > _GetZ(start))
    {
      TRACER_LOG_DEBUG("(PropagateToPoint) Start " + str(start) + " is above end " + str(end) + ", flipping");
      std::swap(start, end);
      flippedStartEnd = true;
    }

    // Seed with straight-line direction (for direct ray)
    DirectionVector seed1 = (end - start).normalized();
    // don't want to start directly horizontal (will not propagate)
    // bump the direction slightly upward
    if (abs(_GetZ(seed1)) < 0.001)
    {
      seed1 = (seed1 + DirectionVector({0, 0, 0.01})).normalized();
    }

    DirectionVector emit1 = seed1;    // initialize
    DirectionVector receive1 = seed1; // initialize

    TRACER_LOG_INFO("(PropagateToPoint) Finding first solution with seed " + str(seed1));
    bool const solution1Found =
        FindEmitAndReceive(start, end, env, seed1, emit1, receive1);
    TRACER_LOG_INFO("(PropagateToPoint) Returned solution dir " + str(emit1));

    // If this one failed, the next one is not worth calculating (shadow zone)
    if (!solution1Found)
    {
      TRACER_LOG_WARNING("(PropagateToPoint) Direct solution finding failed, won't find second solution");
      return signalPaths;
    }

    auto const path1 = GetSignalPath(start, emit1, end, env);
    signalPaths.push_back(flippedStartEnd ? FlipSignalPath(path1) : path1);

    // Seed with weighted average of \hat{z} and the direct ray (for refracted)
    auto const phi =
        atan2(static_cast<double>(seed1.y), static_cast<double>(seed1.x));

    // in most cases, the other solution is close to the negative launch vector
    double cosine = (abs(_GetZ(emit1)) + 0.99) * 0.25;
    auto sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
    auto seed2 = DirectionVector({sine * cos(phi), sine * sin(phi), cosine});

    DirectionVector emit2 = seed2;
    DirectionVector receive2 = seed2;

    TRACER_LOG_INFO("Finding second solution with seed " + str(seed2));
    bool solution2Found = FindEmitAndReceive(start, end, env, seed2, emit2, receive2);
    TRACER_LOG_INFO("Returned solution dir " + str(emit2));

    // if similar, probably found the same solution again
    // do a scan to check for other solution
    if (abs(_GetZ(emit1) - _GetZ(emit2)) < 0.0001)
    {
      TRACER_LOG_WARNING(
          "Solution 2 is similar to solution 1, performing search for another minimum");
      int const nTest = 5;
      double const dcos = (1 - _GetZ(emit1)) / (nTest + 1);
      for (int itest = 1; itest <= nTest; itest++)
      {
        cosine += dcos;
        sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
        seed2 = DirectionVector({sine * cos(phi), sine * sin(phi), cosine});
        solution2Found = FindEmitAndReceive(start, end, env, seed2, emit2, receive2);

        if (abs(_GetZ(emit1) - _GetZ(emit2)) > 0.0001)
        {
          break;
        }
      }
      TRACER_LOG_INFO("After retry, found new solution dir " + str(emit2));
    }

    if (!solution2Found)
    {
      return signalPaths;
    }

    auto const path2 = GetSignalPath(start, emit2, end, env);
    signalPaths.push_back(flippedStartEnd ? FlipSignalPath(path2) : path2);

    return signalPaths;
  }

  inline std::tuple<std::vector<DirectionVector>, std::vector<DirectionVector>>
  RayTracer2D::FindEmitAndReceiveBrent(
      Point const &startInit, Point const &targetInit, EnvironmentBase const &env,
      uint nRaysInit, int maxIterations)
  {

    Point start = startInit;
    Point target = targetInit;
    bool flippedStartEnd = false;

    if (_GetZ(target) > _GetZ(start))
    {
      TRACER_LOG_DEBUG("Flipping start/end positions. original start " +
                       str(startInit) + ", original end " + str(targetInit));
      std::swap(start, target);
      flippedStartEnd = true;
    }

    DirectionVector testDir;
    Point testPos = target;

    auto const dr = target - start;
    auto const phi =
        atan2(static_cast<double>(dr.y), static_cast<double>(dr.x));

    auto makeDirVec = [&](double cosine)
    {
      auto const sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
      return DirectionVector({sine * cos(phi), sine * sin(phi), cosine});
    };

    // construct a distance function
    auto dist2D = [&](double cosine)
    {
      auto const startDir = makeDirVec(cosine);
      auto const success = ShootOneRayToMinimumZ(start, startDir, testPos, testDir, target, env);

      // handle rays that never hit minimum z
      if (!success)
      {
        TRACER_LOG_DEBUG("Failed ray from x0 " + str(start) + " startDir " + str(startDir) + " end at " + str(testPos));
        // throw std::runtime_error("BAD");
        return std::numeric_limits<LengthType>::infinity();
      }

      auto const dist = Get2DRadialDistance(start, target, testPos);
      TRACER_LOG_TRACE("Testing cos = " + str(cosine) + " R_err = " + str(dist));
      return dist;
    };

    std::map<double, LengthType>
        cachedSamples; // holds solutions to avoid duplicate math

    int iteration = 0;

    double cosMin = -1.0;
    double cosMax = 1.0;
    int nRays = nRaysInit;

    TRACER_LOG_INFO("Finding solution using " + std::to_string(nRays) + " rays");

    while (iteration < maxIterations)
    {
      std::vector<double> x_vals(nRays);

      for (int isol = 0; isol < nRays; isol++)
      {
        // set up initial ray direction. Pull solutions from cache if available or calculate
        x_vals[isol] = std::clamp(cosMin + isol * (cosMax - cosMin) / double(nRays - 1.0), -1.0, 1.0);
        if (cachedSamples.find(x_vals[isol]) == cachedSamples.end())
        {
          cachedSamples[x_vals[isol]] = dist2D(x_vals[isol]);
        }
      }

      // Remove entries with infinity
      x_vals.erase(
          std::remove_if(x_vals.begin(), x_vals.end(),
                         [&](double x)
                         {
                           auto it = cachedSamples.find(x);
                           if (it != cachedSamples.end() && std::isinf(it->second))
                           {
                             TRACER_LOG_DEBUG("Bad initial value found, removing cos " + str(it->first) + " from list");
                             cachedSamples.erase(it); // Remove from map
                             return true;             // Mark for removal from vector
                           }
                           return false;
                         }),
          x_vals.end());

      std::vector<DirectionVector> foundEmit;
      std::vector<DirectionVector> foundReceive;

      // collect all cached cosines within current bounds
      std::vector<double> validCosines;
      for (auto const &[cosine, error] : cachedSamples)
      {
        if (cosine >= cosMin && cosine <= cosMax)
        {
          validCosines.push_back(cosine);
        }
      }

      // sort them to ensure adjacent comparisons
      std::sort(validCosines.begin(), validCosines.end());

      // look for sign changes between neighboring consine values
      for (size_t i = 0; i + 1 < validCosines.size(); ++i)
      {
        double x0 = validCosines[i];
        double x1 = validCosines[i + 1];
        double f0 = cachedSamples[x0];
        double f1 = cachedSamples[x1];

        if (f0 * f1 <= 0.0 || std::abs(f0) < STOP_CLOSE_LENGTH || std::abs(f1) < STOP_CLOSE_LENGTH)
        {
          TRACER_LOG_DEBUG("Using brent on cos0 " + str(x0) + " f0 " +
                           str(f0) + " cos1 " + str(x1) + " f1 " + str(f1));
          auto const cosine = BrentMethod(dist2D, x0, x1, f0, f1, DCOS_TOL);
          auto const emit = makeDirVec(cosine);

          foundEmit.push_back((flippedStartEnd ? testDir * -1.0 : emit));
          foundReceive.push_back((flippedStartEnd ? testDir : emit * -1.0));
        }
      }

      if (!foundEmit.empty())
      {
        TRACER_LOG_INFO("Found " + std::to_string(foundEmit.size()) + " solutions");
        return std::make_tuple(foundEmit, foundReceive);
      }

      TRACER_LOG_DEBUG("After iteration " + str(iteration + 1) +
                       ", no solutions found after shooting " + str(cachedSamples.size()) +
                       " total rays");

      // update the cosine interval to keep the top 3 solutions
      if (iteration < maxIterations - 1)
      {
        // choose the three best errors and use that to update the limits of cosMin and cosMax
        std::vector<std::pair<double, LengthType>> sortedSamples(cachedSamples.begin(), cachedSamples.end());

        // sort by absolute error
        std::sort(sortedSamples.begin(), sortedSamples.end(),
                  [](auto const &a, auto const &b)
                  {
                    return std::abs(a.second) < std::abs(b.second);
                  });

        // take up to three best samples
        size_t nBest = std::min<size_t>(3, sortedSamples.size());
        double newCosMin = sortedSamples[0].first;
        double newCosMax = sortedSamples[0].first;

        for (size_t i = 1; i < nBest; ++i)
        {
          newCosMin = std::min(newCosMin, sortedSamples[i].first);
          newCosMax = std::max(newCosMax, sortedSamples[i].first);
        }

        cosMin = std::max(-1.0, newCosMin);
        cosMax = std::min(1.0, newCosMax);
        TRACER_LOG_DEBUG("New cosine limits are " + str(cosMin) + " and " + str(cosMax));
        TRACER_LOG_DEBUG("Closest point at " + str(sortedSamples[0].second));
      }

      if (!iteration)
      {
        nRays += 2;
      }
      iteration++;
    }

    TRACER_LOG_INFO("Solution finding failed after shooting " + str(cachedSamples.size()) +
                    " total rays. Current dcos = " + str((cosMax - cosMin) / float(nRays)));
    return std::make_tuple(std::vector<DirectionVector>{}, std::vector<DirectionVector>{});
  }

  inline bool RayTracer2D::FindEmitAndReceive(
      Point const &sourceXX, Point const &targetXX, EnvironmentBase const &env,
      DirectionVector const &seed, DirectionVector &emit, DirectionVector &receive)
  {

    auto clampCosine = [](double c)
    { return std::max(-1.0, std::min(1.0, c)); };

    auto const normSeed = seed.normalized();
    DirectionVector v1 = normSeed, v2 = normSeed, testDir = normSeed;
    Point testPos = sourceXX; // will be filled with solutions using ShootOneRay

    double dcos = -1e-5;

    Point start = sourceXX;
    Point end = targetXX;
    bool flippedStartEnd = false;

    if (_GetZ(end) > _GetZ(start))
    {
      TRACER_LOG_TRACE("Flipping start/end positions. original start " +
                       str(sourceXX) + ", original end " + str(targetXX));
      std::swap(start, end);
      flippedStartEnd = true;

      ShootOneRayToMinimumZ(start, v1, testPos, testDir, end, env);
      TRACER_LOG_TRACE("Flipping start/end positions. original " + str(testDir.normalized()) +
                       ", flipped " + str(testDir.normalized() * -1));
      v1 = v2 = testDir.normalized() * -1;
    }

    auto const phi =
        std::atan2(v1.y, v1.x);

    if (_GetZ(seed) < -0.99)
    {
      TRACER_LOG_TRACE("Direction close to vertical, setting dcos to " + str(dcos * 0.01));
      dcos *= 0.01;
    }

    auto fireRay = [&](DirectionVector const &dir)
    {
      auto const success = ShootOneRayToMinimumZ(start, dir, testPos, testDir, end, env);
      if (!success)
      {
        return Get2DRadialDistance(start, end, start) - 1.0;
      }

      return Get2DRadialDistance(start, end, testPos);
    };

    LengthType dr1, dr2;

    // calculate metric for the seed
    dr1 = fireRay(v1);
    // TRACER_LOG_TRACE("D1 --> start {}, end {}, testPos {}", start, end, testPos);

    // set launch vector #2 as a slight perturbation to calculate derivatie
    auto cosine = _GetZ(v1) + dcos;
    if (cosine < -1)
    {
      dcos *= -1;
      cosine = _GetZ(v1) + dcos;
    }
    else if (cosine > 1)
    {
      dcos *= -1;
      cosine = _GetZ(v1) + dcos;
    }

    auto sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
    v2 = DirectionVector({sine * cos(phi), sine * sin(phi), cosine});
    dr2 = fireRay(v2);
    auto derivative = dcos / (dr1 - dr2);

    // variables for binary search
    auto closeNegCos = cosine;
    auto closeNegErr = 0.0;
    auto closePosCos = cosine;
    auto closePosErr = 0.0;

    for (int istep = 1; istep <= 30;
         istep++)
    { // Don't run forever (could be in shadow zone)

      // We are well below numerical precision at this point and there is
      // no point in going forward anymore
      if (abs(dcos) < DCOS_TOL)
      {
        TRACER_LOG_TRACE("Numerical precision limit hit, dcos: " + str(dcos));
        emit = flippedStartEnd ? testDir * -1 : v1;
        receive = flippedStartEnd ? v1 * -1 : testDir;
        emit = emit.normalized();
        receive = receive.normalized();
        numericalDerivativeSteps_ += istep;

        TRACER_LOG_DEBUG("FOUND SOLUTION based on dcos machine tolerance:  v0 " + str(emit) +
                         ", pos " + str(testPos) + ", abs 3D dist " +
                         str((testPos - end).getNorm()));
        return true;
      }

      // take a step in cosine (with anealing to prevent oscillations)
      auto cosineStep = derivative * dr1 * pow(0.95, istep - 1.0);

      // keep bounded for large derivatives
      if (cosine + cosineStep > 1)
      {
        TRACER_LOG_TRACE("Step would be above 1, setting cosine to: " + str(0.5 * (1 + cosine)));
        cosine = std::min(0.5 * (1 + cosine), 1.0);
      }
      else if (cosine + cosineStep < -1)
      {
        TRACER_LOG_TRACE("Step would be below 1, setting cosine to: " + str(0.5 * (-1.0 + cosine)));
        cosine = std::max(-1.0, 0.5 * (-1.0 + cosine));
      }
      else
      {
        cosine += cosineStep;
      }

      // shoot the first ray
      cosine = clampCosine(cosine);
      sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
      v1 = DirectionVector({sine * cos(phi), sine * sin(phi), cosine});
      dr1 = fireRay(v1);

      // stop if close enough to boundary
      if (abs(dr1) < STOP_CLOSE_LENGTH)
      {
        emit = flippedStartEnd ? testDir * -1 : v1;
        receive = flippedStartEnd ? v1 * -1 : testDir;
        emit = emit.normalized();
        receive = receive.normalized();

        TRACER_LOG_DEBUG("FOUND SOLUTION within tol, v0 " + str(emit) +
                         ", pos " + str(testPos) + ", abs 3D dist " +
                         str((testPos - end).getNorm()));
        numericalDerivativeSteps_ += istep;
        return true;
      }

      // update the launch vector and shoot second ray
      cosine += dcos;
      cosine = clampCosine(cosine);
      sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
      v2 = DirectionVector({sine * cos(phi), sine * sin(phi), cosine});
      dr2 = fireRay(v2);

      // Update the binary search limits
      if (0.0 > dr1 && (dr1 > closeNegErr || 0.0 == closeNegErr))
      {
        closeNegErr = dr1;
        closeNegCos = _GetZ(v1);
        TRACER_LOG_TRACE("Step " + str(istep) +
                         " --> Updating neg binary limit " + str(closeNegCos) + " " + str(closePosErr));
      }
      else if (0.0 < dr1 && (dr1 < closePosErr || 0.0 == closePosErr))
      {
        closePosErr = dr1;
        closePosCos = _GetZ(v1);
        TRACER_LOG_TRACE("Step " + str(istep) + " --> Updating pos binary limit " +
                         str(closeNegCos) + " " + str(closePosErr));
      }

      TRACER_LOG_TRACE("Step " + str(istep) + " --> Binary search values closeNegCos: " +
                       str(closeNegCos) + " closeNegErr: " + str(closeNegErr) +
                       " closePosCos: " + str(closePosCos) + "  closePosErr: " + str(closePosErr));

      // If spanning the solution, linear interp using current solutions
      if (dr1 * dr2 < 0.0 * 0.0)
      {
        cosine = (dr1 * _GetZ(v2) - dr2 * _GetZ(v1)) / (dr1 - dr2);
        derivative *= 0.0; // don't update on next step
        dcos *= 0.1;
      }
      else if (closeNegErr * closePosErr < 0.0 * 0.0)
      {
        // cosine = 0.5 * (closeNegCos + closePosCos);
        cosine = (closePosErr * closeNegCos - closeNegErr * closePosCos) /
                 (closePosErr - closeNegErr);
        derivative *= 0.0;
        dcos *= 0.5;
        TRACER_LOG_TRACE(
            "Step " + str(istep) + " --> Taking binary step current cos: " +
            str(_GetZ(v1)) + " next cos: " + str(cosine));
      }
      else
      {

        if (dr1 == dr2)
        {
          TRACER_LOG_TRACE("Found the same points, derivative inf, returning seed value");
          v1 = seed;
          break;
        }

        // numerical derivative
        derivative = dcos / (dr1 - dr2);
        cosine = clampCosine(cosine - dcos);

        TRACER_LOG_TRACE("Step " + str(istep) + " --> phi " + str(phi) + ", cosine " +
                         str(cosine) + ", new-cosine " + str(cosine + dcos) + ", derivative " + str(derivative));
      }

      // decrease the step size if the numerical derivative is no longer a small step
      if (abs(dr1) < abs(dr1 - dr2) * 100 && dcos > 1e-10)
      {
        dcos *= 0.05;
        TRACER_LOG_TRACE("Close to the target dr1 " + str(abs(dr1)) + " diff " +
                         str(abs(dr1 - dr2)) + " dcos " + str(dcos));
      }

      TRACER_LOG_DEBUG("Step " + str(istep) + " --> v1 " + str(v1) + ", v2 " + str(v2) +
                       ", dr1 " + str(dr1) + ", dr2 " + str(dr2));
    }

    TRACER_LOG_TRACE("Could not find prop dir, v1 " + str(v1) + ", v2 " + str(v2) +
                     ", dr1 " + str(dr1) + ", dr2 " + str(dr2));

    emit = flippedStartEnd ? testDir * -1 : v1;
    receive = flippedStartEnd ? v1 * -1 : testDir;
    emit = emit.normalized();
    receive = receive.normalized();

    TRACER_LOG_DEBUG("No solution found, last check was v0 " + str(emit) +
                     ", pos " + str(testPos) + ", abs 3D dist " +
                     str((testPos - end).getNorm()));

    return false;
  }

  inline void RayTracer2D::FindIntersectionWithPlane(
      Point const &x0, DirectionVector const &v0, Point &end, DirectionVector &endDir,
      Plane const &plane, LengthType const step, EnvironmentBase const &env)
  {

    auto const currentSteps = stepsTaken_;

    TRACER_LOG_TRACE("x0 " + str(x0) + ", v0 " + str(v0));
    constexpr uint kMaxInitSteps = 10; // number of initial steps for bounding root

    LengthType initStepSize = step;
    auto const f0 = DistToPlane(plane, x0);
    auto f1 = DistToPlane(plane, end);

    uint initSteps = 0;
    while (f0 * f1 > 0)
    {
      if (initSteps >= kMaxInitSteps)
      {
        TRACER_LOG_ERROR("Failed to find bounds that surround plane! \nx0: " +
                         str(x0) + " dist " + str(f0) +
                         "\nend: " + str(end) + " dist " + str(f1) +
                         "\nv0: " + str(v0) + " step size: " + str(step) +
                         "\nplane center " + str(plane.getCenter()) + " normal " +
                         str(plane.getNormal()));
        throw std::runtime_error("Failed finding plane bounding points");
      }

      initSteps++;
      initStepSize *= 2.0;
      TakeAdaptiveStep(x0, v0, end, endDir, initStepSize, env, false);
      f1 = DistToPlane(plane, end);
    }

    LengthType const brentStep = initStepSize;
    LengthType testStep = initStepSize;
    // define stepping function to find the root of (scaling step size)
    auto root = [&](double mult)
    {
      testStep = brentStep * mult;
      TakeAdaptiveStep(x0, v0, end, endDir, testStep, env, false);
      return DistToPlane(plane, end);
    };

    auto const frac = BrentMethod(root, 0.0, 1.0, f0, f1, STOP_CLOSE_LENGTH);
    testStep = brentStep * frac;
    TakeAdaptiveStep(x0, v0, end, endDir, testStep, env, false); // take found step

    planeIntersectionSteps_ += stepsTaken_ - currentSteps;
  }

  inline std::tuple<double, double> RayTracer2D::TransmitThroughPlane(
      Point const &x0, DirectionVector const &v0, Point &end, DirectionVector &endDir,
      Plane const &plane, LengthType const step, EnvironmentBase const &env)
  {
    TRACER_LOG_TRACE("Performing transmission starting at x0: " + str(x0) + " v0: " +
                     str(v0) + " ending at xf " + str(end) + " vf: " + str(endDir));

    FindIntersectionWithPlane(x0, v0, end, endDir, plane, step, env);
    TRACER_LOG_ALL("Intersecting plane at xf " + str(end) + " vf: " + str(endDir));

    auto const n1 = env.get_n(x0);

    // get a point just on the other side of the plane based on current distance
    // to the plane plus 0.1% overshoot
    auto const farPos = x0 - plane.getNormal() * 1.001 * DistToPlane(plane, x0);
    auto const n2 = env.get_n(farPos);
    TRACER_LOG_ALL("Refrac. Index  n1 " + str(n1) + ", n2: " + str(n2));

    auto const cosineTheta1 = abs(endDir.normalized().dot(plane.getNormal()));
    auto const sineSquareTheta1 = 1.0 - cosineTheta1 * cosineTheta1;

    // init variables
    double transmitSComp = 1.0;
    double transmitPComp = 1.0;

    // check total-intertal-reflection case
    auto const ratio = n1 / n2;
    auto const sineSquareTheta2 = ratio * ratio * sineSquareTheta1; // Snell's law
    // CORSIKA_LOG_TRACE("SinSqTh2 {}", sineSquareTheta2);
    if (sineSquareTheta2 <= 1.0)
    {
      auto const cosineTheta2 = sqrt(1.0 - sineSquareTheta2);
      transmitSComp = 2 * n1 * cosineTheta1 /
                      (n1 * cosineTheta1 + n2 * cosineTheta2);
      transmitPComp = 2 * n1 * cosineTheta1 /
                      (n2 * cosineTheta1 + n1 * cosineTheta2);
    }

    auto const vv = v0.normalized();
    endDir = vv - plane.getNormal() * vv.dot(plane.getNormal()) +
             plane.getNormal() * std::sqrt(n2 * n2 - n1 * n1 + pow(vv.dot(plane.getNormal()), 2.0));
    endDir = endDir.normalized();
    TRACER_LOG_ALL("After Snell's Law dir changed from " + str(vv) + " to " + str(endDir));

    // Want to ensure that we are on the correct side of the plane after reflection
    // Ensure this is the case by forcing it to the other side
    auto const dist = DistToPlane(plane, end);
    if (dist * DistToPlane(plane, x0) > 0.0)
    { // wrong side of plane
      // CORSIKA_LOG_DEBUG("On the wrong side, current: {}, dir: {}", end, endDir);
      TRACER_LOG_ALL("On the wrong side of the plane " + str(dist) + " " + str(DistToPlane(plane, x0)));
      end = end - endDir.normalized() * dist * 1.001 / endDir.normalized().dot(plane.getNormal());
      TRACER_LOG_ALL("After correction " + str(end));
    }

    TRACER_LOG_TRACE("After transmission, xf " + str(end) + " vf: " + str(endDir));
    return {transmitSComp, transmitPComp};
  }

  inline std::tuple<double, double> RayTracer2D::ReflectOffPlane(
      Point const &x0, DirectionVector const &v0, Point &end, DirectionVector &endDir,
      Plane const &plane, LengthType const step, EnvironmentBase const &env)
  {

    TRACER_LOG_TRACE("Performing reflection starting at x0: " + str(x0) + " v0: " +
                     str(v0) + " ending at xf " + str(end) + " vf: " + str(endDir));
    FindIntersectionWithPlane(x0, v0, end, endDir, plane, step, env);

    TRACER_LOG_TRACE("Intersecting plane at xf " + str(end) + " vf: " + str(endDir));

    auto const n1 = env.get_n(x0);

    // get a point just on the other side of the plane based on current distance
    // to the plane plus 0.1% overshoot
    auto const farPos =
        x0 - plane.getNormal() * 1.001 * (x0 - plane.getCenter()).dot(plane.getNormal());
    auto const n2 = env.get_n(farPos);

    // CORSIKA_LOG_TRACE("Far pos: {} n1: {} n2: {}", farPos, n1, n2);

    auto const cosineTheta1 = abs(endDir.normalized().dot(plane.getNormal()));
    auto const sineSquareTheta1 = 1.0 - cosineTheta1 * cosineTheta1;

    // CORSIKA_LOG_TRACE("CosineTh1 {} SinSqTh1 {}", cosineTheta1, sineSquareTheta1);

    // init variables
    double reflectSComp = 1.0;
    double reflectPComp = 1.0;

    // check total-intertal-reflection case
    auto const ratio = n1 / n2;
    auto const sineSquareTheta2 = ratio * ratio * sineSquareTheta1; // Snell's law
    // CORSIKA_LOG_TRACE("SinSqTh2 {}", sineSquareTheta2);
    if (sineSquareTheta2 <= 1.0)
    {
      auto const cosineTheta2 = sqrt(1.0 - sineSquareTheta2);
      reflectSComp = (n1 * cosineTheta1 - n2 * cosineTheta2) /
                     (n1 * cosineTheta1 + n2 * cosineTheta2);
      reflectPComp = (n2 * cosineTheta1 - n1 * cosineTheta2) /
                     (n1 * cosineTheta2 + n2 * cosineTheta1);
    }

    // perform reflection
    endDir = endDir - plane.getNormal() * 2 * endDir.dot(plane.getNormal());

    // CORSIKA_LOG_TRACE("Fresnell perp: {} Fresnell parallel: {}", reflectSComp,
    //                   reflectPComp);
    // CORSIKA_LOG_TRACE("New direction: {}\n", endDir);

    // Want to ensure that we are on the correct side of the plane after reflection
    // Ensure this is the case by forcing it to the other side
    auto distanceFromPlane = (end - plane.getCenter()).dot(plane.getNormal());
    if (distanceFromPlane < 0.0)
    { // wrong side of plane
      // CORSIKA_LOG_DEBUG("On the wrong side, current: {}, dir: {}", end, endDir);
      end = end - endDir.normalized() * distanceFromPlane * 1.001 / endDir.normalized().dot(plane.getNormal());
    }

    TRACER_LOG_TRACE("New direction xf " + str(end) + " vf: " + str(endDir));
    return {reflectSComp, reflectPComp};
  }

  inline bool RayTracer2D::ShootOneRayToMinimumZ(Point const &start,
                                                 DirectionVector const &startDir,
                                                 Point &end, DirectionVector &endDir,
                                                 Point const &target,
                                                 EnvironmentBase const &env)
  {
    auto const minimumPlane =
        Plane(target, DirectionVector({0, 0, 1}));

    auto distFunc = [&](Point const &x)
    {
      return DistToPlane(minimumPlane, x);
    };

    SignalPath const sigPath = ShootOneRay(start, startDir, target, env, distFunc, false, MAX_RAY_STEPS);

    end = sigPath.getEnd();
    endDir = sigPath.receive_;

    return abs(distFunc(end)) < 2 * PLANE_CLOSE_TOL;
  }

  inline bool RayTracer2D::ShootOneRayToMaximumR(Point const &start,
                                                 DirectionVector const &startDir,
                                                 Point &end, DirectionVector &endDir,
                                                 Point const &target,
                                                 EnvironmentBase const &env)
  {
    auto distFunc = [&](Point const &x)
    {
      return Get2DRadialDistance(start, target, x);
    };

    SignalPath const sigPath = ShootOneRay(start, startDir, target, env, distFunc, false, MAX_RAY_STEPS);

    end = sigPath.getEnd();
    endDir = sigPath.receive_;

    return abs(distFunc(end)) < 2 * PLANE_CLOSE_TOL;
    ;
  }

  inline SignalPath RayTracer2D::GetSignalPath(Point const &start,
                                               DirectionVector const &startDir,
                                               Point const &target,
                                               EnvironmentBase const &env)
  {
    auto const minimumPlane =
        Plane(target, DirectionVector({0, 0, 1}));

    auto distFunc = [&](Point const &x)
    {
      return DistToPlane(minimumPlane, x);
    };

    if (distFunc(start) < 0.0)
    {
      TRACER_LOG_WARNING("Asking for path that starts at " + str(start) +
                         " which is below target " + str(target));
    }

    return ShootOneRay(start, startDir, target, env, distFunc, true, MAX_RAY_STEPS);
  }

  inline Vec3 RayTracer2D::Get2DProjection(Point const &x0, Point const &xf)
  {
    Vec3 const dr = xf - x0;
    return dr - axis_ * dr.dot(axis_);
  }

  inline LengthType RayTracer2D::Get2DRadialDistance(Point const &x0, Point const &xf,
                                                     Point const &xtest)
  {
    Vec3 const proj = Get2DProjection(x0, xf);
    return (xtest - xf).dot(proj / proj.getNorm());
  }

  inline void RayTracer2D::PrintProfiling()
  {
    LOG_INFO("               Total steps taken: " + std::to_string(stepsTaken_));
    LOG_INFO("  Steps to find correct ray path: " + std::to_string(numericalDerivativeSteps_));
    LOG_INFO("Steps finding reflection surface: " + std::to_string(planeIntersectionSteps_));
    LOG_INFO("           Total rays propagated: " + std::to_string(raysPropagated_));
  }

  inline void RayTracer2D::ResetProfiling()
  {
    numericalDerivativeSteps_ = 0;
    raysPropagated_ = 0;
    stepsTaken_ = 0;
    planeIntersectionSteps_ = 0;
  }

  inline SignalPath RayTracer2D::ShootOneRay(Point const &start, DirectionVector const &startDir,
                                             Point const &target, EnvironmentBase const &env,
                                             std::function<LengthType(Point const &)> distMethod,
                                             bool const recordPath, uint const maxSteps)
  {
    raysPropagated_++;

    TRACER_LOG_TRACE("x0 " + str(start) + ", v0 " + str(startDir) + ", target " + str(target));

    Path path(start);

    // set of stepping variables
    auto x0 = start;
    auto end = start;
    auto v0 = startDir;
    auto endDir = startDir;
    auto const normStartDir = startDir.normalized();
    LengthType stepSize = INITAL_STEP_SIZE;
    LengthType prevStep = stepSize;

    // path variables
    LengthType weightedIndexOfRefraction = 0.0;
    auto propLength = 0.0;
    double fresnelS = 1.0;
    double fresnelP = 1.0;

    auto accumulateSegment = [&](Point const &a, Point const &b)
    {
      auto const diff = b - a;
      auto const dist = diff.getNorm();
      auto const avgPos = a + diff * 0.5;
      auto const n = env.get_n(avgPos);
      weightedIndexOfRefraction += n * dist;
      propLength += dist;
    };

    uint istep = 0;

    while (distMethod(end) * distMethod(x0) >= 0 && istep < maxSteps)
    {

      if (istep > 0)
      {
        accumulateSegment(x0, end);
        if (recordPath)
        {
          path.addToEnd(end);
        }
      }

      // update current locations
      x0 = end;
      v0 = endDir;
      prevStep = stepSize;

      TakeAdaptiveStep(x0, v0, end, endDir, stepSize, env, true);

      // check for plane crossings
      for (auto const &plane : reflectionLayers_)
      {
        auto const signFinal = DistToPlane(plane, end);
        if (signFinal * DistToPlane(plane, x0) < 0.0) // step crosses the plane
        {
          if (signFinal * DistToPlane(plane, target) < 0.0)
          {
            stepSize = prevStep;
            TRACER_LOG_TRACE("Entered reflection loop between " +
                             str(x0) + " and " + str(end) + ", step size " + str(stepSize));
            auto [tempS, tempP] =
                ReflectOffPlane(x0, v0, end, endDir, plane, stepSize, env);
            fresnelS *= tempS;
            fresnelP *= tempP;
          }
          else
          {
            TRACER_LOG_TRACE("Entered transmission loop between " +
                             str(x0) + " and " + str(end) + ", step size " + str(stepSize));
            auto [tempS, tempP] =
                TransmitThroughPlane(x0, v0, end, endDir, plane, stepSize, env);
            fresnelS *= tempS;
            fresnelP *= tempP;
          }
        }

        istep++;
      }
    }

    // if we actually found a potential solution (i.e. didn't run out of steps)
    // find the zero-crossing of the distMethod function
    if (distMethod(end) * distMethod(x0) <= 0)
    {
      TRACER_LOG_TRACE("Stopping criteria bounded, x0: " + str(x0) +
                       ", v0: " + str(v0) + ", xf: " + str(end) + ", vf: " + str(endDir));

      LengthType currentStepSize;

      // define stepping function to find the root of (scaling step size)
      auto root = [&](double mult)
      {
        currentStepSize = stepSize * mult;
        TakeAdaptiveStep(x0, v0, end, endDir, currentStepSize, env, false);
        return distMethod(end);
      };

      auto const frac = BrentMethod(root, 0.0, 1.0, distMethod(x0), distMethod(end), STOP_CLOSE_LENGTH);
      currentStepSize = stepSize * frac;
      TakeAdaptiveStep(x0, v0, end, endDir, currentStepSize, env, false); // take found step
    }
    else
    {
      TRACER_LOG_TRACE("Stopping criteria not found! Number of steps taken " + str(istep));
    }

    path.addToEnd(end);
    accumulateSegment(x0, end);

    auto const avgIndex = weightedIndexOfRefraction / propLength;
    auto const propTime = weightedIndexOfRefraction / constants::c;

    auto const n0 = env.get_n(start);
    auto const nf = env.get_n(end);

    TRACER_LOG_TRACE("Solution found ending at " + str(end));

    return SignalPath(propTime, avgIndex, n0, nf, normStartDir, endDir.normalized(),
                      propLength, path, fresnelS, fresnelP);
  }

  void RayTracer2D::TakeAdaptiveStep(Vec3 const &startPos, Vec3 const &startDir, Vec3 &endPos,
                                     Vec3 &endDir, double &h0, EnvironmentBase const &env,
                                     bool updateStep)
  {
    tracer_.AdaptiveStep(startPos, startDir, endPos, endDir, h0, env, updateStep);
    stepsTaken_++;
  }

} // namespace c8_tracer
