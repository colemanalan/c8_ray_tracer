#pragma once

#include "c8_tracer/logger.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"
#include "c8_tracer/transcribed/RayTracer.hpp"

#define CLOSE_TOL 0.00001
#define DCOS_TOL 1e-9
#define INITAL_STEP_SIZE 0.1
#define PLANE_CLOSE_TOL 0.000005

namespace c8_tracer
{

#define str std::to_string

  SignalPath FlipSignalPath(SignalPath const &inPath)
  {
    Path flippedPath(inPath.getEnd());

    for (int ibin = inPath.getNSegments() - 2; ibin >= 0; ibin--)
    {
      flippedPath.addToEnd(inPath.getPoint(ibin));
    }

    return SignalPath(inPath.propagation_time_, inPath.average_refractive_index_,
                      inPath.refractive_index_destination_,
                      inPath.refractive_index_source_, inPath.receive_ * -1,
                      inPath.emit_ * -1, inPath.R_distance_, flippedPath,
                      inPath.fresnelS_, inPath.fresnelP_);
  }

  inline RayTracer2D::RayTracer2D(DirectionVector const axis,
                                  LengthType const minStep = 0.0001,
                                  LengthType const maxStep = 10.0,
                                  double const tolerance = 1e-8)
      : tracer_(minStep, maxStep, tolerance), axis_(axis.normalized()) {}

  inline void RayTracer2D::AddReflectionLayer(Plane const layer)
  {
    // Check to ensure that this remains a 2D problem
    if (std::abs(layer.getNormal().dot(axis_)) < cos(0.001 * M_PI / 180.))
    {
      logger.error(
          "Cannot add a reflection layer that is not perpendicular"
          "the axis. Layer normal " +
          str(layer.getNormal()) + ", ray tracer symmetry axis " + str(axis_));
      throw std::invalid_argument("Bad reflection layer");
    }

    reflectionLayers_.push_back(layer);
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

    if (end.z > start.z)
    {
      logger.debug("(PropagateToPoint) Start " + str(start) + " is above end " + str(end) + ", flipping");
      std::swap(start, end);
      flippedStartEnd = true;
    }

    // Seed with straight-line direction (for direct ray)
    DirectionVector seed1 = (end - start).normalized();
    // don't want to start directly horizontal (will not propagate)
    // bump the direction slightly upward
    if (abs(seed1.z) < 0.001)
    {
      seed1 = (seed1 + DirectionVector({0, 0, 0.01})).normalized();
    }

    DirectionVector emit1 = seed1;    // initialize
    DirectionVector receive1 = seed1; // initialize

    logger.info("(PropagateToPoint) Finding first solution with seed " + str(seed1));
    bool const solution1Found =
        FindEmitAndReceive(start, end, env, seed1, emit1, receive1);
    logger.info("(PropagateToPoint) Returned solution dir " + str(emit1));

    // If this one failed, the next one is not worth calculating (shadow zone)
    if (!solution1Found)
    {
      logger.warning("(PropagateToPoint) Direct solution finding failed, won't find second solution");
      return signalPaths;
    }

    auto const path1 = GetSignalPath(start, emit1, end, env);
    signalPaths.push_back(flippedStartEnd ? FlipSignalPath(path1) : path1);

    // Seed with weighted average of \hat{z} and the direct ray (for refracted)
    auto const phi =
        atan2(static_cast<double>(seed1.y), static_cast<double>(seed1.x));

    // in most cases, the other solution is close to the negative launch vector
    double cosine = (abs(emit1.z) + 0.99) * 0.25;
    auto sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
    auto seed2 = DirectionVector({sine * cos(phi), sine * sin(phi), cosine});

    DirectionVector emit2 = seed2;
    DirectionVector receive2 = seed2;

    logger.info("(PropagateToPoint) Finding second solution with seed " + str(seed2));
    bool solution2Found = FindEmitAndReceive(start, end, env, seed2, emit2, receive2);
    logger.info("(PropagateToPoint) Returned solution dir " + str(emit2));

    // if similar, probably found the same solution again
    // do a scan to check for other solution
    if (abs(emit1.z - emit2.z) < 0.0001)
    {
      logger.warning(
          "(PropagateToPoint) Solution 2 is similar to solution 1, performing search for another minimum");
      int const nTest = 5;
      double const dcos = (1 - emit1.z) / (nTest + 1);
      for (int itest = 1; itest <= nTest; itest++)
      {
        cosine += dcos;
        sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
        seed2 = DirectionVector({sine * cos(phi), sine * sin(phi), cosine});
        solution2Found = FindEmitAndReceive(start, end, env, seed2, emit2, receive2);

        if (abs(emit1.z - emit2.z) > 0.0001)
        {
          break;
        }
      }
      logger.info("(PropagateToPoint) After retry, found new solution dir " + str(emit2));
    }

    if (!solution2Found)
    {
      return signalPaths;
    }

    auto const path2 = GetSignalPath(start, emit2, end, env);
    signalPaths.push_back(flippedStartEnd ? FlipSignalPath(path2) : path2);

    return signalPaths;
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

    if (end.z > start.z)
    {
      logger.trace("(PropagateToPoint) Flipping start/end positions. original start " +
                   str(sourceXX) + ", original end " + str(targetXX));
      std::swap(start, end);
      flippedStartEnd = true;

      ShootOneRayToMinimumZ(start, v1, testPos, testDir, end, env);
      logger.trace("(PropagateToPoint) Flipping start/end positions. original " + str(testDir.normalized()) +
                   ", flipped " + str(testDir.normalized() * -1));
      v1 = v2 = testDir.normalized() * -1;
    }

    auto const phi =
        std::atan2(v1.y, v1.x);

    if (seed.z < -0.99)
    {
      logger.trace("(PropagateToPoint) Direction close to vertical, setting dcos to " + str(dcos * 0.01));
      dcos *= 0.01;
    }

    uint nSteps = 0;
    auto fireRay = [&](DirectionVector const &dir)
    {
      nSteps += ShootOneRayToMinimumZ(start, dir, testPos, testDir, end, env);
      return Get2DRadialDistance(start, end, testPos);
    };

    LengthType dr1, dr2;

    // calculate metric for the seed
    dr1 = fireRay(v1);
    // logger.trace("D1 --> start {}, end {}, testPos {}", start, end, testPos);

    // set launch vector #2 as a slight perturbation to calculate derivatie
    auto cosine = v1.z + dcos;
    if (cosine < -1)
    {
      dcos *= -1;
      cosine = v1.z + dcos;
    }
    else if (cosine > 1)
    {
      dcos *= -1;
      cosine = v1.z + dcos;
    }

    auto sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
    v2 = DirectionVector({sine * cos(phi), sine * sin(phi), cosine});
    dr2 = fireRay(v2);
    auto derivative = dcos / (dr1 - dr2);
    // logger.trace("D2 --> start " + str(start) + ", end " + str(end) + ", testPos " + str(testPos));

    // logger.trace("Start of propagate v1 {}, v2 {}, dr1 {}, dr2 {}", v1, v2, dr1,
    //               dr2);
    // logger.trace("Start of propagate derivative {}", derivative);
    // logger.trace("Start cosine value will be {}", cosine + derivative * dr1);

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
        logger.trace("(PropagateToPoint) Numerical precision limit hit, dcos: " + str(dcos));
        emit = flippedStartEnd ? testDir * -1 : v1;
        receive = flippedStartEnd ? v1 * -1 : testDir;
        emit = emit.normalized();
        receive = receive.normalized();
        numericalDerivativeSteps_ += istep;

        logger.debug("(PropagateToPoint) FOUND SOLUTION based on dcos machine tolerance:  v0 " + str(emit) +
                     ", pos " + str(testPos) + ", abs 3D dist " +
                     str((testPos - end).getNorm()));
        return true;
      }

      // take a step in cosine (with anealing to prevent oscillations)
      auto cosineStep = derivative * dr1 * pow(0.95, istep - 1.0);

      // keep bounded for large derivatives
      if (cosine + cosineStep > 1)
      {
        logger.trace("(PropagateToPoint) Step would be above 1, setting cosine to: " + str(0.5 * (1 + cosine)));
        cosine = std::min(0.5 * (1 + cosine), 1.0);
      }
      else if (cosine + cosineStep < -1)
      {
        logger.trace("(PropagateToPoint) Step would be below 1, setting cosine to: " + str(0.5 * (-1.0 + cosine)));
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
      if (abs(dr1) < CLOSE_TOL)
      {
        emit = flippedStartEnd ? testDir * -1 : v1;
        receive = flippedStartEnd ? v1 * -1 : testDir;
        emit = emit.normalized();
        receive = receive.normalized();

        logger.debug("(PropagateToPoint) FOUND SOLUTION within tol, v0 " + str(emit) +
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
        closeNegCos = v1.z;
        logger.trace("(PropagateToPoint) Step " + str(istep) +
                     " --> Updating neg binary limit " + str(closeNegCos) + " " + str(closePosErr));
      }
      else if (0.0 < dr1 && (dr1 < closePosErr || 0.0 == closePosErr))
      {
        closePosErr = dr1;
        closePosCos = v1.z;
        logger.trace("(PropagateToPoint) Step " + str(istep) + " --> Updating pos binary limit " +
                     str(closeNegCos) + " " + str(closePosErr));
      }

      logger.trace("(PropagateToPoint) Step " + str(istep) + " --> Binary search values closeNegCos: " +
                   str(closeNegCos) + " closeNegErr: " + str(closeNegErr) +
                   " closePosCos: " + str(closePosCos) + "  closePosErr: " + str(closePosErr));

      // If spanning the solution, linear interp using current solutions
      if (dr1 * dr2 < 0.0 * 0.0)
      {
        cosine = (dr1 * v2.z - dr2 * v1.z) / (dr1 - dr2);
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
        logger.trace(
            "(PropagateToPoint) Step " + str(istep) + " --> Taking binary step current cos: " +
            str(v1.z) + " next cos: " + str(cosine));
      }
      else
      {

        if (dr1 == dr2)
        {
          logger.trace("(PropagateToPoint) Found the same points, derivative inf, returning seed value");
          v1 = seed;
          break;
        }

        // numerical derivative
        derivative = dcos / (dr1 - dr2);
        cosine = clampCosine(cosine - dcos);

        logger.trace("(PropagateToPoint) Step " + str(istep) + " --> phi " + str(phi) + ", cosine " +
                     str(cosine) + ", new-cosine " + str(cosine + dcos) + ", derivative " + str(derivative));
      }

      // decrease the step size if the numerical derivative is no longer a small step
      if (abs(dr1) < abs(dr1 - dr2) * 100 && dcos > 1e-10)
      {
        dcos *= 0.05;
        logger.trace("(PropagateToPoint) Close to the target dr1 " + str(abs(dr1)) + " diff " +
                     str(abs(dr1 - dr2)) + " dcos " + str(dcos));
      }

      logger.debug("(PropagateToPoint) Step " + str(istep) + " --> v1 " + str(v1) + ", v2 " + str(v2) +
                   ", dr1 " + str(dr1) + ", dr2 " + str(dr2));
    }

    logger.trace("(PropagateToPoint) Could not find prop dir, v1 " + str(v1) + ", v2 " + str(v2) +
                 ", dr1 " + str(dr1) + ", dr2 " + str(dr2));

    numericalDerivativeSteps_ += nSteps;
    emit = flippedStartEnd ? testDir * -1 : v1;
    receive = flippedStartEnd ? v1 * -1 : testDir;
    emit = emit.normalized();
    receive = receive.normalized();

    logger.debug("(PropagateToPoint) No solution found, last check was v0 " + str(emit) +
                 ", pos " + str(testPos) + ", abs 3D dist " +
                 str((testPos - end).getNorm()));

    return false;
  }
  inline void RayTracer2D::FindIntersectionWithPlane(
      Point const &x0, DirectionVector const &v0, Point &end, DirectionVector &endDir,
      Plane const &plane, LengthType const step, EnvironmentBase const &env)
  {

    logger.trace("(FindIntersectionWithPlane) x0 " + str(x0) + ", v0 " + str(v0));

    constexpr auto closeThreshold = PLANE_CLOSE_TOL; // defines when solution is found
    constexpr uint kMaxInitSteps = 10;               // number of initial steps for bounding root
    constexpr uint kMaxBinarySteps = 50;             // max steps to find solutions
    constexpr auto epsilon = 1e-10;                  // tolerance for division

    // initializing values that will be passed by reference
    double fracHigh = 1.0;
    double fracMid = 0.0;
    double fracLow = 0.0;
    LengthType currentStepSize;

    // lambda for calculating distance from the plane
    auto computeDist = [&](Point const &pt)
    {
      return (pt - plane.getCenter()).dot(plane.getNormal());
    };

    // lambda for taking one step forward and bumping the counter
    auto takeStep = [&](double frac)
    {
      currentStepSize = frac * step;
      stepsTaken_++;
      tracer_.AdaptiveStep(x0, v0, end, endDir, currentStepSize, env, false);
    };

    auto distLow = computeDist(x0);
    takeStep(fracHigh);
    auto distHigh = computeDist(end);

    uint counter = 0;

    // Sanity check to ensure that the solution is bounded in the first place
    // Take increasingly larger steps if need be
    while (distLow * distHigh >= 0.0 * 0.0 && counter++ < kMaxInitSteps)
    {
      // CORSIKA_LOG_DEBUG(
      //     "Failed not-bounded check\n x0 {}\n v0 {}\n fracLow {}"
      //     "\n distLow {}\n fracHigh {}\n distHigh {}",
      //     x0, v0, fracLow, distLow, fracHigh, distHigh);

      fracHigh += 0.5;
      takeStep(fracHigh);
      distHigh = computeDist(end);
    }

    if (counter >= kMaxInitSteps)
    { // Can happen where ray tracer gets stuck due to
      // large jumps in n
      // CORSIKA_LOG_WARN(
      //     "Could not properly initialize binary search to find intersection with plane"
      //     "step frac: {}, distance to plane: {}. Will use straight track approximation "
      //     "to find intersection",
      //     fracHigh, distHigh);
      // CORSIKA_LOG_DEBUG(
      //     "\nx0 {}\nv0 {}\nstep {}\nend {}\nendDir {}\nfracLow {}\n"
      //     "distLow {}\n fracHigh {}\n distHigh {}",
      //     x0, v0, step, end, endDir, fracLow, distLow, fracHigh, distHigh);

      auto denom = endDir.dot(plane.getNormal());
      if (abs(denom) > epsilon)
      {
        auto const t = (plane.getCenter() - end).dot(plane.getNormal()) / denom;
        end = end + endDir * t;
      }

      return;
    }

    // set to a value in the middle
    fracMid = 0.5 * fracHigh;
    counter = 0;

    // Begin binary search
    while (counter++ < kMaxBinarySteps)
    {
      binarySearchReflection_++;
      takeStep(fracMid);
      auto distMid = computeDist(end);

      if (abs(distLow) < closeThreshold)
      {
        logger.trace("(FindIntersectionWithPlane) Finishing at low point! " + str(distLow));
        takeStep(fracLow);
        break;
      }
      else if (abs(distMid) < closeThreshold)
      {
        logger.trace("(FindIntersectionWithPlane) Finishing at mid point! " + str(distMid));
        takeStep(fracMid);
        break;
      }
      else if (abs(distHigh) < closeThreshold)
      {
        logger.trace("(FindIntersectionWithPlane) Finishing at high point! " + str(distHigh));
        takeStep(fracHigh);
        break;
      }

      // CORSIKA_LOG_DEBUG("Low: {}, Med: {}, High: {}", distLow, distMid, distHigh);

      if (distLow * distMid < 0.0 * 0.0)
      {
        logger.trace("(FindIntersectionWithPlane) Point is btw low/mid " + str(distLow * distMid));
        auto delta = distLow - distMid;
        fracHigh = fracMid;
        distHigh = distMid;
        fracMid = (abs(delta) > epsilon * 1.0f)
                      ? (distLow * fracMid - distMid * fracLow) / delta
                      : 0.5 * (fracLow + fracMid);
      }
      else if (distHigh * distMid < 0.0 * 0.0)
      {
        logger.trace("(FindIntersectionWithPlane) Point is btw mid/high " + str(distHigh * distMid));
        auto delta = distMid - distHigh;
        fracLow = fracMid;
        distLow = distMid;
        fracMid = (abs(delta) > epsilon * 1.0f)
                      ? (distMid * fracHigh - distHigh * fracMid) / delta
                      : 0.5 * (fracMid + fracHigh);
      }
      else
      {
        // happens when you exactly hit the plane
        logger.trace("(FindIntersectionWithPlane) Finishing at mid point! " + str(distMid));
        takeStep(fracMid);
        break;
      }

    } // End while loop

    if (counter >= kMaxBinarySteps)
    {
      // CORSIKA_LOG_DEBUG(
      //     "Could not find reflection after {} steps, deviation {}, tolerance {}", counter,
      // computeDist(end), closeThreshold);
    }
  }

  inline std::tuple<double, double> RayTracer2D::ReflectOffPlane(
      Point const &x0, DirectionVector const &v0, Point &end, DirectionVector &endDir,
      Plane const &plane, LengthType const step, EnvironmentBase const &env)
  {

    logger.trace("(ReflectOffPlane) Performing reflection starting at x0: " + str(x0) + " v0: " +
                 str(v0) + " ending at xf " + str(end) + " vf: " + str(endDir));
    FindIntersectionWithPlane(x0, v0, end, endDir, plane, step, env);

    logger.trace("(ReflectOffPlane) Intersecting plane at xf " + str(end) + " vf: " + str(endDir));

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

    logger.trace("(ReflectOffPlane) New direction xf " + str(end) + " vf: " + str(endDir));
    return {reflectSComp, reflectPComp};
  }

  void RayTracer2D::FindRadius(Point const &x0, DirectionVector const &v0, Point &end,
                               DirectionVector &endDir, LengthTypeSq const rMaxSq,
                               LengthType const step, EnvironmentBase const &env)
  {
    auto fracHigh = 2.0; // Extend a bit beyond for numerical precision safety
    auto fracMid = 1.0;
    auto fracLow = 0.0;

    LengthType testStep; // Need this to pass by reference

    testStep = fracLow * step;
    stepsTaken_++;
    tracer_.AdaptiveStep(x0, v0, end, endDir, testStep, env, false);
    auto distLow = rMaxSq - Get2DProjection(x0, end).getSquaredNorm();

    testStep = fracHigh * step;
    stepsTaken_++;
    tracer_.AdaptiveStep(x0, v0, end, endDir, testStep, env, false);
    auto distHigh = rMaxSq - Get2DProjection(x0, end).getSquaredNorm();

    // Begin binary search
    while (true)
    {
      binarySearchBoundary_++;

      // only ever need to calculate the mid point
      testStep = fracMid * step;
      stepsTaken_++;
      tracer_.AdaptiveStep(x0, v0, end, endDir, testStep, env, false);
      auto distMid = rMaxSq - Get2DProjection(x0, end).getSquaredNorm();

      // CORSIKA_LOG_TRACE("Low: {}, Med: {}, High: {}", distLow, distMid, distHigh);

      if (distHigh * distMid < 0.0 * 0.0 * 0.0 * 0.0)
      { // between mid and high
        // CORSIKA_LOG_TRACE("Point is btw mid/high {}", distHigh * distMid);

        // check stopping criteria
        if (abs(distHigh) < 0.001f * 0.001f)
        {
          // CORSIKA_LOG_TRACE("Finishing at high point! {}", distHigh);
          testStep = fracHigh * step;
          stepsTaken_++;
          tracer_.AdaptiveStep(x0, v0, end, endDir, testStep, env, false);
          break;
        }
        else if (abs(distMid) < 0.001f * 0.001f)
        {
          // CORSIKA_LOG_TRACE("Finishing at mid point! {}", distMid);
          testStep = fracMid * step;
          stepsTaken_++;
          tracer_.AdaptiveStep(x0, v0, end, endDir, testStep, env, false);
          break;
        }

        // update mid and low points using linear interpolation
        auto const temp = fracMid;
        fracMid = (distMid * fracHigh - distHigh * fracMid) / (distMid - distHigh);
        fracMid = 0.5 * (fracHigh + fracMid);
        fracLow = temp;
        distLow = distMid;
        // CORSIKA_LOG_TRACE("New limits will be low {} mid {} high {}", fracLow, fracMid,
        // fracHigh);
      }
      else if (distLow * distMid < 0.0 * 0.0 * 0.0 * 0.0)
      { // between low and mid
        // CORSIKA_LOG_TRACE("Point is btw low/mid {}", distLow * distMid);

        // check stopping criteria
        if (abs(distLow) < 0.0001 * 0.0001)
        {
          // CORSIKA_LOG_TRACE("Finishing at low point! {}", distLow);
          testStep = fracLow * step;
          stepsTaken_++;
          tracer_.AdaptiveStep(x0, v0, end, endDir, testStep, env, false);
          break;
        }
        else if (abs(distMid) < 0.0001 * 0.0001)
        {
          // CORSIKA_LOG_TRACE("Finishing at mid point! {}", distMid);
          testStep = fracMid * step;
          stepsTaken_++;
          tracer_.AdaptiveStep(x0, v0, end, endDir, testStep, env, false);
          break;
        }

        // update high and mid points using linear interpolation
        auto const temp = fracMid;
        fracMid = (distLow * fracMid - distMid * fracLow) / (distLow - distMid);
        fracHigh = temp;
        distHigh = distMid;
        // CORSIKA_LOG_TRACE("New limits will be low {} mid {} high {}", fracLow, fracMid,
        // fracHigh);
      }
      else
      {
        // shouldn't ever happen, but here for safety
        // CORSIKA_LOG_DEBUG("Point is no longer bounded in binary search!");
        // CORSIKA_LOG_DEBUG("Low h={} d={}, Mid h={} d={}, High h={} d={},", fracLow * step,
        // distLow, fracMid * step, distMid, fracHigh * step, distHigh);
        assert(false);
        break;
      }
    } // End while loop
  }

  inline uint RayTracer2D::ShootOneRayToMinimumZ(Point const &start,
                                                 DirectionVector const &startDir,
                                                 Point &end, DirectionVector &endDir,
                                                 Point const &minZ,
                                                 EnvironmentBase const &env)
  {

    raysPropagated_++;

    // Initialize state
    Point x0 = end = start;
    DirectionVector v0 = endDir = startDir;

    auto const minimumPlane =
        Plane(minZ, DirectionVector({0, 0, 1}));

    LengthType stepSize = INITAL_STEP_SIZE;
    LengthType prevStep = stepSize;

    logger.trace("(ShootOneRayToMinimumZ) x0 " + str(start) + ", v0 " + str(startDir) + ", minZ " + str(minZ));

    uint constexpr maxSteps = 10000;
    uint istep = 0;

    do
    {
      x0 = end;
      v0 = endDir;
      prevStep = stepSize;

      tracer_.AdaptiveStep(x0, v0, end, endDir, stepSize, env);
      istep++;

      // check for plane crossings
      for (auto const &plane : reflectionLayers_)
      {
        if ((end - plane.getCenter()).dot(plane.getNormal()) < 0.0)
        {
          stepSize = prevStep;

          logger.trace("(ShootOneRayToMinimumZ) Entered reflection loop between " +
                       str(x0) + " and " + str(end) + ", step size " + str(stepSize));
          if ((x0 - plane.getCenter()).dot(plane.getNormal()) < 0.0)
          {
            logger.error("(ShootOneRayToMinimumZ) First point should not be on this side of plane " +
                         str(x0) + "  " + str((x0 - plane.getCenter()).dot(plane.getNormal())));
            throw std::invalid_argument("Wrong!");
          }
          ReflectOffPlane(x0, v0, end, endDir, plane, stepSize, env);
        }
      }

      if (istep == maxSteps)
      {
        logger.debug("(ShootOneRayToMinimumZ) max steps have been used, stopping");
        break;
      }
    } while ((end - minimumPlane.getCenter()).dot(minimumPlane.getNormal()) >
             0.0); // Stop when you passed the antennas

    // RAY_TRACER_PRINT("Stopped after {} steps, step size {}", istep, stepSize);
    // RAY_TRACER_PRINT("Prev, x0 {}, v0 {}", x0, v0);
    // RAY_TRACER_PRINT("Shot, xf {}, vf {}", end, endDir);
    // RAY_TRACER_PRINT("Will find crossing point minZ at {}, distance {}", minZ,
    //                  (end - minimumPlane.getCenter()).dot(minimumPlane.getNormal()));

    // Find the crossing point
    if (maxSteps > istep)
    {
      logger.trace("(ShootOneRayToMinimumZ) point found below the plane x0: " + str(x0) +
                   ", v0: " + str(v0) + ", xf: " + str(end) + ", vf: " + str(endDir));
      FindIntersectionWithPlane(x0, v0, end, endDir, minimumPlane, stepSize, env);
    }
    // RAY_TRACER_PRINT("Corrected, xf {}, vf {}, dist {}", end, endDir,
    //                  (end - minimumPlane.getCenter()).dot(minimumPlane.getNormal()));

    endDir = endDir.normalized();
    return istep;
  }

  inline uint RayTracer2D::ShootOneRayToMaximumR(Point const &start,
                                                 DirectionVector const &startDir,
                                                 Point &end, DirectionVector &endDir,
                                                 LengthTypeSq const rMaxSq,
                                                 EnvironmentBase const &env)
  {

    raysPropagated_++;

    // set up some variables
    auto x0 = end = start;
    auto v0 = endDir = startDir;

    LengthType stepSize = INITAL_STEP_SIZE; // This will be updated to keep good tolerance

    // CORSIKA_LOG_DEBUG("Shooting ray, x0 {}, v0 {}, r2 {}", start, startDir, rMaxSq);

    uint const maxSteps = 10000;
    uint istep = 0;
    do
    {
      x0 = end;
      v0 = endDir;
      tracer_.AdaptiveStep(x0, v0, end, endDir, stepSize, env);

      istep++;

      // check for plane crossings
      for (auto const &plane : reflectionLayers_)
      {
        if ((end - plane.getCenter()).dot(plane.getNormal()) < 0.0)
        {
          // CORSIKA_LOG_DEBUG("Entered reflection loop between {} and {}, step size {}", x0,
          //                   end, stepSize);
          ReflectOffPlane(x0, v0, end, endDir, plane, stepSize, env);
        }
      }

      // std::cout << x0 << std::endl;

      if (istep == maxSteps)
      {
        // CORSIKA_LOG_DEBUG("{} steps have been used, stopping, dir {}", istep, startDir);
        break;
      }
    } while (Get2DProjection(start, end).getSquaredNorm() <
             rMaxSq); // Stop when you passed the antennas

    // CORSIKA_LOG_TRACE("Stopped after {} steps, step size {}", istep, stepSize);
    // CORSIKA_LOG_TRACE("Prev, x0 {}, v0 {}, r2 {}", x0, v0,
    //                   Get2DProjection(start, x0).getSquaredNorm());
    // CORSIKA_LOG_TRACE("Shot, xf {}, vf {}, r2 {}", end, endDir,
    //                   Get2DProjection(start, end).getSquaredNorm());
    // CORSIKA_LOG_TRACE("Will find the edge from {} with distance {} {} {}", x0, rMaxSq,
    //                   Get2DProjection(start, x0).getSquaredNorm(),
    //                   rMaxSq - Get2DProjection(start, x0).getSquaredNorm());

    // Take one last partial step to get close to the propagation limit
    if (maxSteps > istep)
    {
      auto dr = sqrt(rMaxSq) - sqrt(Get2DProjection(start, x0).getSquaredNorm());
      FindRadius(x0, v0, end, endDir, dr * dr, stepSize, env);
    }
    // CORSIKA_LOG_DEBUG("Corrected, xf {}, vf {}, r2 {}", end, endDir,
    //                   Get2DProjection(start, end).getSquaredNorm());

    endDir = endDir.normalized();

    return istep;
  }

  inline SignalPath RayTracer2D::GetSignalPath(Point const &start,
                                               DirectionVector const &startDir,
                                               Point const &target,
                                               EnvironmentBase const &env)
  {

    Path path(start);
    logger.debug("(GetSignalPath) Making signal path using emit: " + str(startDir) + " from " + str(start) + " to " + str(target));

    auto const targetZ = target.z;

    // set up some variables
    auto x0 = start;
    auto end = start;
    auto v0 = startDir;

    auto const normStartDir = startDir.normalized();
    DirectionVector endDir = normStartDir;

    if (start.z < targetZ)
    {
      logger.warning("(GetSignalPath) Asking to propagate with an initial position below the target");
    }

    LengthType constexpr initialStepSize = INITAL_STEP_SIZE;
    LengthType stepSize = initialStepSize; // This will be updated to keep good tolerance

    LengthType weightedIndexOfRefraction = 0.0;
    auto propLength = 0.0;

    double fresnelS = 1.0;
    double fresnelP = 1.0;

    uint constexpr maxSteps = 500000;
    uint istep = 0;

    auto accumulateSegment = [&](Point const &a, Point const &b)
    {
      auto const diff = b - a;
      auto const dist = diff.getNorm();
      auto const avgPos = a + diff * 0.5;
      auto const n = env.get_n(avgPos);
      weightedIndexOfRefraction += n * dist;
      propLength += dist;
    };

    logger.debug("(GetSignalPath) Going from " + str(end.z) + " to " + str(targetZ));


    while (end.z > targetZ && istep < maxSteps)
    {
      if (istep > 0)
      {
        path.addToEnd(end);
        accumulateSegment(x0, end);
      }

      x0 = end;
      v0 = endDir;

      tracer_.AdaptiveStep(x0, v0, end, endDir, stepSize, env);
      ++istep;

      for (auto const &plane : reflectionLayers_)
      {
        if ((end - plane.getCenter()).dot(plane.getNormal()) < 0.0)
        {
          auto [tempS, tempP] =
              ReflectOffPlane(x0, v0, end, endDir, plane, stepSize, env);
          fresnelS *= tempS;
          fresnelP *= tempP;
        }
      }
    }

    logger.debug("(GetSignalPath) Finding final intersection with plane using x0:" + str(x0) + ", v0: " + str(v0));
    // Final intersection if not capped
    if (istep < maxSteps)
    {
      auto const minimumPlane = Plane(target, DirectionVector({0, 0, 1}));
      FindIntersectionWithPlane(x0, v0, end, endDir, minimumPlane, stepSize, env);
    }

    path.addToEnd(end);
    accumulateSegment(x0, end);

    auto const avgIndex = weightedIndexOfRefraction / propLength;
    auto const propTime = weightedIndexOfRefraction / constants::c;

    auto const n0 = env.get_n(start);
    auto const nf = env.get_n(end);

    return SignalPath(propTime, avgIndex, n0, nf, normStartDir, endDir.normalized(),
                      propLength, path, fresnelS, fresnelP);
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
    logger.info("  Steps to find correct ray path " + std::to_string(numericalDerivativeSteps_));
    logger.info("           Total rays propagated " + std::to_string(raysPropagated_));
    logger.info("               Total steps taken " + std::to_string(stepsTaken_));
    logger.info("Steps finding reflection surface " + std::to_string(binarySearchReflection_));
    logger.info("    Steps finding stopping point " + std::to_string(binarySearchBoundary_));
  }

  inline void RayTracer2D::ResetProfiling()
  {
    numericalDerivativeSteps_ = 0;
    raysPropagated_ = 0;
    stepsTaken_ = 0;
    binarySearchReflection_ = 0;
    binarySearchBoundary_ = 0;
  }

} // namespace c8_tracer
