#pragma once

#include <functional>

#include "c8_tracer/environment.hpp"
#include "c8_tracer/plane.hpp"
#include "c8_tracer/ray_tracer_base.hpp"
#include "c8_tracer/transcribed/CashKarpIntegrator.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"

namespace c8_tracer
{

  class RayTracer2D : public RayTracerBase
  {
  public:
    /*
    Finds the ray solutions between sets of points in space.
    This is implemented to solve problems where the gradient of the index of refraction
    only changes along one dimension.

    Arguments:
      axis:
        axis along which the gradent changes
      minStep:
        minimum step size that can be take by the propagator
      maxStep:
        maximum step size that can be take by the propagator
      tolerance:
        relative uncertainty on the numerical integration for each step. The step size
    will be adjusted to keep it within this tolerance
      brentRays:
        how many initial rays to cast when creating a Brent-Method-based search

    */
    RayTracer2D(DirectionVector const axis, LengthType const minStep,
                LengthType const maxStep, double const tolerance, int const brentRays);
    ~RayTracer2D() {}

    /*
    Adds a reflection layer for the propagation. Rays only reflect such that they
    stay on the side of the plane where the normal vector is pointing
    */
    void AddReflectionLayer(Plane const layer);

    /*
    This is the main function which finds the solutions from `start` to `end`

    Arguments:
      start:
        starting location of the launch point
      end:
        the target location for the propagation
      env:
        description of the refractive index and gradient
    */
    std::vector<SignalPath> PropagateToPoint(Point const &start, Point const &end,
                                             EnvironmentBase const &env) override;

    std::vector<SignalPath> GetSignalPaths(Point const &start, Point const &end,
                                           EnvironmentBase const &env) override;

    /*
    Finds the initial propagation direction such that a ray propagates from `start` to
    `end`. This does a simple numerical gradient descent starting from `seed` until
    convergence is found. If a solution is not found after several steps, this will return
    0

    Arguments:
      start:
        starting location of the launch point
      end:
        the target location for the propagation
      env:
        description of the refractive index and gradient
      seed:
        initial guess of the launch vector
      emit:
        this will be updated with the found emit vector
      receive:
        this will be updated with the found receive vector

    Returns:
      if a valid solution was found
    */
    bool FindEmitAndReceive(Point const &start, Point const &end, EnvironmentBase const &env,
                            DirectionVector const &seed, DirectionVector &emit,
                            DirectionVector &receive) override;

    /*
    Finds the initial propagation direction such that a ray propagates from `start` to
    `end`. This shoots several initial rays to find rays that bound the solutions and
    then used the Brent Method to ultimately find the optimal path(s). If no solutions
    are found after the initial rays are cast, the angular bounds will be updated
    to keep the interval that includes the top 3 solutions and casts another set of nRays

    Arguments:
      start:
        starting location of the launch point
      end:
        the target location for the propagation
      env:
        description of the refractive index and gradient
      nRays:
        number of rays to shoot each iteration
      maxIterations:
        number of times that rays will be shot again if no solutions are found

    Returns:
      tuple of a list of emit vectors and a list of receive vectors.
    */
    std::tuple<std::vector<DirectionVector>, std::vector<DirectionVector>>
    FindEmitAndReceiveBrent(
        Point const &start, Point const &end, EnvironmentBase const &env, uint nRays,
        int maxIterations);

    /*
    Propagates a ray from `start` until the lateral distance (see `Get2DProjection`)
    exceeds a specified value. This function handles the reflections off the known planes.

    Arguments:
      start:
        starting location of the launch point
      startDir:
        initial direction of propagation at `start`
      end:
        this value will be updated with the final location
      endDir:
        this value will be update with direction at `end`
      target:
        the point that defines the maximum radius
      env:
        description of the refractive index and gradient

    Returns:
      the number of steps in the propagation from `start` to `end`
    */
    bool ShootOneRayToMaximumR(Point const &start, DirectionVector const &startDir,
                               Point &end, DirectionVector &endDir, Point const &target,
                               EnvironmentBase const &env);

    /*
    Propagates a ray from `start` until reaching a minimum z value
    This function handles the reflections off the known planes.

    Arguments:
      start:
        starting location of the launch point
      startDir:
        initial direction of propagation at `start`
      end:
        this value will be updated with the final location
      endDir:
        this value will be update with direction at `end`
      target:
        point at which the minimum z value will be considered
      env:
        description of the refractive index and gradient

    Returns:
      the number of steps in the propagation from `start` to `end`
    */
    bool ShootOneRayToMinimumZ(Point const &start, DirectionVector const &startDir,
                               Point &end, DirectionVector &endDir, Point const &target,
                               EnvironmentBase const &env);

    /*
    This does a binary search to find the point on the plane and that the ray
    intersects. This function expects to only need to perform a single step
    and optimizes the step size such that the ray arrives at the plane

    Arguments:
      x0:
        starting location of the launch point (should be only 1 step away from the
        boundary)
      startDir:
        initial direction of propagation at `x0`
      end:
        this value will be updated with the final location (will be on the reflection
        plane)
      endDir:
        this value will be update with direction at `end` after reflecting it
      plane:
        plane that is doing the reflecting
      step:
        maximum step size in the binary search
      env:
        description of the refractive index and gradient
    */
    void FindIntersectionWithPlane(Point const &x0, DirectionVector const &v0, Point &end,
                                   DirectionVector &endDir, Plane const &plane,
                                   LengthType const step, EnvironmentBase const &env);

    /*
    Finds the intersection with the plane using `FindIntersectionWithPlane` and gives
    the final direction after reflection

    Arguments:
      x0:
        starting location of the launch point (should be only 1 step away from the
        boundary)
      startDir:
        initial direction of propagation at `x0`
      end:
        this value will be updated with the final location (will be on the reflection
        plane)
      endDir:
        this value will be update with direction at `end` after reflecting it
      plane:
        plane that is doing the reflecting
      step:
        maximum step size in the binary search
      env:
        description of the refractive index and gradient

    Returns:
      tuple of (Fresnel-S, Fresnel-P)
    */
    std::tuple<double, double> ReflectOffPlane(Point const &x0, DirectionVector const &v0,
                                               Point &end, DirectionVector &endDir,
                                               Plane const &plane, LengthType const step,
                                               EnvironmentBase const &env);

    /*
    Finds the intersection with the plane using `FindIntersectionWithPlane` and gives
    the final direction after transmission

    Arguments:
      x0:
        starting location of the launch point (should be only 1 step away from the
        boundary)
      startDir:
        initial direction of propagation at `x0`
      end:
        this value will be updated with the final location (will be on the reflection
        plane)
      endDir:
        this value will be update with direction at `end` after reflecting it
      plane:
        plane that is doing the reflecting
      step:
        maximum step size in the binary search
      env:
        description of the refractive index and gradient

    Returns:
      tuple of (Fresnel-S, Fresnel-P)
    */
    std::tuple<double, double> TransmitThroughPlane(Point const &x0, DirectionVector const &v0,
                                                    Point &end, DirectionVector &endDir,
                                                    Plane const &plane, LengthType const step,
                                                    EnvironmentBase const &env);

    /*
    This helps find the point that is exactly a specified lateral distance (see
    `Get2DProjection`) away from `x0`. The function performs a binary search to find
    the step size required such that the propagation is at the correct distance. This
    function expects to only need to perform a single step and optimizes the step size
    such that the ray arrives at the boundary

    Arguments:
      x0:
        starting location of the launch point (should be only 1 step away from the
        boundary)
      v0:
        initial direction of propagation at `x0`
      end:
        this value will be updated with the final location (will be sqrt(`rMaxSq`) from
      endDir:
        this value will be update with direction at `end` point
      rMaxSq:
        plane that is doing the reflecting step: maximum step size in the binary
        search

      env: description of the refractive index and gradient
    */
    void FindRadius(Point const &x0, DirectionVector const &v0, Point &end,
                    DirectionVector &endDir, LengthTypeSq const rMaxSq,
                    LengthType const step, EnvironmentBase const &env);

    /*
    Gets the projection in the x-y plane of the vector from `x0` to `x1`
    The x-y plane is defined in the coordinate system of `x0`
    */
    Vec3 Get2DProjection(Point const &x0, Point const &xf);

    /*
    Gets the distance in the x-y plane along the line connecting x0 and xf
    The radial distance from xf and xtest is returned
    */
    LengthType Get2DRadialDistance(Point const &x0, Point const &xf, Point const &xtest);

    SignalPath GetSignalPath(Point const &start, DirectionVector const &startDir,
                             Point const &target, EnvironmentBase const &env) override;

    void PrintProfiling();
    void ResetProfiling();

    LengthType _GetZ(Point const &p) const { return p.dot(axis_); }

  private:
    CashKarpIntegrator tracer_;
    std::vector<Plane> reflectionLayers_;

    int nRays_;

    // profiling parameters
    unsigned long int numericalDerivativeSteps_ = 0;
    unsigned long int stepsTaken_ = 0;
    unsigned long int planeIntersectionSteps_ = 0;
    unsigned long int raysPropagated_ = 0;

    LengthType DistToPlane(Plane const &p, Point const &x) { return (x - p.getCenter()).dot(p.getNormal()); }

    // Main function for propagating rays forward and handling effects during propagation
    SignalPath ShootOneRay(Point const &start, DirectionVector const &startDir,
                           Point const &target, EnvironmentBase const &env,
                           std::function<LengthType(Point const &)> stopMethod,
                           bool const accumulate, uint const maxSteps);

    // Wrapper for the Adaptive step of the ray tracer
    void TakeAdaptiveStep(Vec3 const &startPos, Vec3 const &startDir, Vec3 &endPos,
                      Vec3 &endDir, double &h0, EnvironmentBase const &env,
                      bool updateStep = true);

  };
}

#include "c8_tracer/transcribed/RayTracer.inl"
