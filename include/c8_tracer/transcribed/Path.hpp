#pragma once

#include <deque>

#include "c8_tracer/transcribed/c8_typedefs.hpp"

namespace c8_tracer
{

  /**
   * This class represents a (potentially) curved path between two
   * points using N >= 1 straight-line segments.
   */
  class Path
  {

  protected:
    std::deque<Point> points_; ///< The points that make up this path.
    LengthType length_ = 0.0;  ///< The length of the path.

    using iterator = std::deque<Point>::iterator;
    using const_iterator = std::deque<Point>::const_iterator;

  public:
    /**
     * Create a Path with a given starting Point.
     */
    Path(Point const &point);

    /**
     * Initialize a Path from an existing collection of Points.
     */
    Path(std::deque<Point> const &points);

    /**
     * Add a new Point to the end of the path.
     */
    inline void addToEnd(Point const &point);

    /**
     * Remove a point from the end of the path.
     */
    inline void removeFromEnd();

    /**
     * Get the total length of the path.
     */
    inline LengthType getLength() const;

    /**
     * Get the starting point of the path.
     */
    inline Point const &getStart() const;

    /**
     * Get the end point of the path.
     */
    inline Point const &getEnd() const;

    /**
     * Get a specific point of the path.
     */
    inline Point const &getPoint(std::size_t const index) const;

    /**
     * Return an iterator to the start of the Path.
     */
    inline const_iterator begin() const;

    /**
     * Return an iterator to the end of the Path.
     */
    inline const_iterator end() const;

    /**
     * Return an iterator to the start of the Path.
     */
    inline iterator begin();

    /**
     * Return an iterator to the end of the Path.
     */
    inline iterator end();

    /**
     * Get the number of steps in the path.
     * This is one less than the number of points that
     * defines the path.
     */
    inline int getNSegments() const;

    inline std::string to_string() const
    {
      return "Path[" + std::to_string(getNSegments()) +
             " steps, start: " + std::to_string(getStart()) +
             ", end: " + std::to_string(getEnd()) + "]";
    }

  }; // class Path

  inline std::ostream &operator<<(std::ostream &os, const Path &path)
  {
    os << path.to_string();
    return os;
  }

} // namespace c8_tracer

#include "c8_tracer/transcribed/Path.inl"
