#pragma once

#include "c8_tracer/environment.hpp"
#include "c8_tracer/transcribed/c8_typedefs.hpp"
#include "c8_tracer/signal_path.hpp"

namespace c8_tracer
{

  enum SolutionMethod
  {
    NGD = 0,
    Brent = 1
  };

  class RayTracerBase
  {
  protected:
    DirectionVector axis_;

  public:
    explicit RayTracerBase(const DirectionVector &axis)
        : axis_(axis) {}
    ~RayTracerBase() = default;

    DirectionVector const &GetAxis() const { return axis_; }

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
    virtual std::vector<SignalPath> PropagateToPoint(Point const &start, Point const &end,
                                                     EnvironmentBase const &env) = 0;

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
    virtual bool FindEmitAndReceive(Point const &start, Point const &end, EnvironmentBase const &env,
                                    DirectionVector const &seed, DirectionVector &emit,
                                    DirectionVector &receive) = 0;
    /*
    Gets the distance in the x-y plane along the line connecting x0 and xf
    The radial distance from xf and xtest is returned
    */
    virtual SignalPath GetSignalPath(Point const &start, DirectionVector const &startDir,
                                     Point const &target, EnvironmentBase const &env) = 0;

    /*
    This is the replacement of PropagateToPoint to get ray paths from `start` to `end`

    Arguments:
      start:
        starting location of the launch point
      end:
        the target location for the propagation
      env:
        description of the refractive index and gradient
    */
    virtual std::vector<SignalPath> GetSignalPaths(Point const &start, Point const &end,
                                                   EnvironmentBase const &env) = 0;
  };

}