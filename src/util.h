/*
 * -------------------------------------------
 * Copyright (c) 2021 - 2025 Prashant K. Jha
 * -------------------------------------------
 * https://github.com/CEADpx/multiphysics-peridynamics
 *
 * Distributed under the Boost Software License, Version 1.0. (See accompanying
 * file LICENSE)
 */
#pragma once

#include <libmesh/point.h>
#include <vector>
#include <cmath>

// tolerance for float comparison
#define COMPARE_EPS 1e-5

namespace util {

/*!
 * @brief Returns true if a > b
 * @param a Value a
 * @param b Value b
 * @return True if a is definitely greater than b
 */
inline bool isGreater(const double &a, const double &b) {
    return (a - b) > ((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) *
      COMPARE_EPS);
};

/*!
 * @brief Returns true if a < b
 * @param a Value a
 * @param b Value b
 * @return True if a is definitely less than b
 */
inline bool isLess(const double &a, const double &b) {
    return (b - a) > ((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) *
      COMPARE_EPS);
};

/*!
 * @brief Computes hat function at given point
 *
 * Hat function:
 *      f ^
 *        |
 *        |
 *     1  o
 *        |           /|\
 *        |         /  |  \
 *        |       /    |    \
 *        |     /      |      \
 *        |   /        |        \
 *        | /          |          \
 *        o____________o____________o______\ x
 *                                         /
 *      x_min                      x_max
 *
 * @param x Point in real line
 * @param x_min Left side point in real line
 * @param x_max Right side point in real line
 * @return value Evaluation of hat function at x
 */
inline double hatFunction(const double &x, const double &x_min, const double &x_max) {

    if (util::isGreater(x, x_min - 1.0E-12) and
        util::isLess(x, x_max + 1.0E-12)) {
  
      double x_mid = 0.5 * (x_min + x_max);
      double l = x_mid - x_min;
  
      // check if this is essentially a point load (dirac)
      if (l < 1.0E-12)
        return 1.0;
  
      if (util::isLess(x, x_mid))
        return (x - x_min) / l;
      else
        return (x_max - x) / l;
    } else
      return 0.0;
  };

/*!
 * @brief Computes hat function at given point
 *
 * This version does not test if point x is in valid interval.
 *
 * Hat function:
 *      f ^
 *        |
 *        |
 *     1  o
 *        |           /|\
 *        |         /  |  \
 *        |       /    |    \
 *        |     /      |      \
 *        |   /        |        \
 *        | /          |          \
 *        o____________o____________o______\ x
 *                                         /
 *      x_min                      x_max
 *
 * @param x Point in real line
 * @param x_min Left side point in real line
 * @param x_max Right side point in real line
 * @return value Evaluation of hat function at x
 */
inline double hatFunctionQuick(const double &x, const double &x_min,
                        const double &x_max) {

                            double x_mid = 0.5 * (x_min + x_max);
                            double l = x_mid - x_min;
                          
                            // check if this is essentially a point load (dirac)
                            if (l < 1.0E-12)
                              return 1.0;
                          
                            if (util::isLess(x, x_mid))
                              return (x - x_min) / l;
                            else
                              return (x_max - x) / l;
                          };

/*!
 * @brief Compute linear step function
 *
 * Step function:
 *
 * f ^
 *   |             __________
 *   |            /
 *   |           /
 *   |   _______/
 *   |  /
 *   | /
 *   |/_________________________ t
 *       x1   x1+x2
 *
 *  - Linear (with slope 1) in [0,l1), constant in [l1,l1+l2)
 *  - Periodic with periodicity l1+l2
 *
 * @param x  Point in real line
 * @param x1 Point such that function is linear with slope 1 in [0, x1)
 * @param x2 Point such that function is constant in [x1, x1 + x2)
 * @return value Evaluation of step function at x
 */
inline double linearStepFunc(const double &x, const double &x1, const double &x2) {

    //
    // a = floor(x/(x1+x2))
    // xl = a * (x1 + x2), xm = xl + x1, xr = xm + x2
    // fl = a * x1
    //
    // At xl, value of the function is = period number \times x1
    //
    // From xl to xm, the function grows linear with slope 1
    // so the value in between [xl, xm) will be
    // fl + (x - xl) = a*x1 + (x - a*(x1+x2)) = x - a*x2
    //
    // In [xm, xr) function is constant and the value is
    // fl + (xm - xl) = a*x1 + (xl + x1 - xl) = (a+1)*x1
  
    double period = std::floor(x / (x1 + x2));
  
    if (util::isLess(x, period * (x1 + x2) + x1))
      return x - period * x2;
    else
      return (period + 1.) * x1;
  };

  inline double linearStepSpecifiedConstFunc(const double &x, const double &x1, const double &x2, const double &const_val) {

    //
    // a = floor(x/(x1+x2))
    // xl = a * (x1 + x2), xm = xl + x1, xr = xm + x2
    // fl = a * x1
    //
    // At xl, value of the function is = period number \times x1
    //
    // From xl to xm, the function grows linear with slope 1
    // so the value in between [xl, xm) will be
    // fl + (x - xl) = a*x1 + (x - a*(x1+x2)) = x - a*x2
    //
    // In [xm, xr) function is constant and the value is const_val
  
    double period = std::floor(x / (x1 + x2));
  
    if (util::isLess(x, period * (x1 + x2) + x1))
      return x - period * x2;
    else
      return const_val;
  };

/*!
 * @brief Compute gaussian function in 1-d
 *
 * Guassian (1-d) function: \f$ f(r) = a \exp(-\frac{r^2}{\beta}). \f$
 *
 * Here \f$ a\f$ is the amplitude and \f$ \beta \f$ is the exponential factor.
 *
 * @param r Distance from origin
 * @param a Amplitude
 * @param beta Factor in exponential function
 * @return value Component of guassian 1-d function
 */
inline double gaussian(const double &r, const double &a, const double &beta) {
    return a * std::exp(-r * r / beta);
};

/*!
 * @brief Compute gaussian function in 2-d
 *
 * Guassian (2-d) function:
 * \f[ f(x,y) = (f_1(x,y), f_2(x,y)), \f]
 * where
 * \f[ f_1(x,y) =  a \exp(-\frac{(x-x_c)^2 + (y-y_c)^2}{\beta}) d_1,
 * \quad f_1(x,y) =  a \exp(-\frac{(x-x_c)^2 + (y-y_c)^2}{\beta}) d_2.
 * \f]
 * Here \f$ (x_c,y_c) \f$ is the center of the pulse, \f$ a\f$ is the
 * amplitude, \f$ \beta \f$ is the exponential factor, and \f$ (d_1,d_2)\f$
 * is the direction of the pulse.
 *
 * @param x  Coordinates of point
 * @param params List of parameters
 * @param dof Component of guassian function
 * @return value Component of guassian 2-d vector function along dof
 */
inline double gaussian2d(const libMesh::Point &x, const size_t &dof,
                  const std::vector<double> &params) {
    if (params.size() < 6) {
        std::cerr << "Error: Not enough parameters to compute guassian 2-d "
                        "function.\n";
        exit(1);
        }
    
        return util::gaussian(
                    (x - libMesh::Point(params[0], params[1], 0.)).norm(), params[5],
                    params[4]) *
                params[2 + dof];
    
};

/*!
 * @brief Compute sum of two gaussian function in 2-d
 *
 * Double guassian (2-d) function:
 * \f[ f(x,y) = (f_1(x,y), f_2(x,y)) + (g_1(x,y), g_2(x,y)), \f]
 * where \f$ (f_1,f_2)\f$ and \f$(g_1, g_2)\f$ are two guassian 2-d function
 * as described in guassian2d() with different values of \f$ (x_c, y_c), a,
 * (d_1, d_2)\f$.
 *
 * @param x  Coordinates of point
 * @param params List of parameters
 * @param dof Component of guassian function
 * @return value Component of guassian 2-d vector function along dof
 */
inline double doubleGaussian2d(const libMesh::Point &x, const size_t &dof,
                        const std::vector<double> &params)  {

                            if (params.size() < 10) {
                              std::cerr << "Error: Not enough parameters to compute guassian 2-d "
                                           "function.\n";
                              exit(1);
                            }
                          
                            return util::gaussian(
                                       (x - libMesh::Point(params[0], params[1], 0.)).norm(), params[9],
                                       params[8]) *
                                       params[4 + dof] +
                                   util::gaussian(
                                       (x - libMesh::Point(params[2], params[3], 0.)).norm(), params[9],
                                       params[8]) *
                                       params[6 + dof];
                          };

/*!
 * @brief Compute harmonic mean of m1 and m2
 *
 * @param m1 Mass 1
 * @param m2 Mass 2
 * @return m Harmonic mean
 */
inline double equivalentMass(const double &m1, const double &m2) {
    return 2. * m1 * m2 / (m1 + m2);
};

inline double harmonicMean(const double &m1, const double &m2) {
    return 2. * m1 * m2 / (m1 + m2);
};

/*!
 * @brief Adds a unique element to a STL container (e.g., std::vector, std::list) of elements of class T.
 *
 * If the element already exists in the container (using operator==), it is not added.
 *
 * @tparam Container STL container type (e.g., std::vector<T>)
 * @tparam T Element type
 * @param container The container to add to
 * @param elem The element to add
 * @return true if the element was added, false if it already existed
 */
template <typename Container, typename T>
inline bool addUnique(Container& container, const T& elem) {
    if (std::find(container.begin(), container.end(), elem) == container.end()) {
        container.push_back(elem);
        return true;
    }
    return false;
}

inline void setFixity(const size_t &i, const unsigned int &dof,
  const bool &flag, std::vector<uint8_t> &fixed) {

// to set i^th bit as true of integer a,
// a |= 1UL << (i % 8)

// to set i^th bit as false of integer a,
// a &= ~(1UL << (i % 8))

flag ? (fixed[i] |= 1UL << dof) : (fixed[i] &= ~(1UL << dof));
}


inline bool isDofFree(const size_t &i, const unsigned int &dof, const std::vector<uint8_t> &fixed) {

  // below checks if d_fix has 1st bit (if dof=0), 2nd bit (if dof=1), 3rd
  // bit (if dof=2) is set to 1 or 0. If set to 1, then it means it is fixed,
  // and therefore it returns false
  return !(fixed[i] >> dof & 1UL);
};

/**
 * @name Rotation and Transformation Functions
 */
/**@{*/

/*!
 * @brief Rotates a vector in xy-plane in clockwise direction
 * @param x Point
 * @param theta Angle
 * @return Point after rotation
 */
inline std::vector<double> rotateCW2D(const std::vector<double> &x,
                               const double &theta) {
  return std::vector<double>{x[0] * std::cos(theta) + x[1] * std::sin(theta),
                            -x[0] * std::sin(theta) + x[1] * std::cos(theta),
                            0.0};
}

/*!
 * @brief Rotates a vector in xy-plane in clockwise direction
 * @param x Point
 * @param theta Angle
 * @return Point after rotation
 */
inline libMesh::Point rotateCW2D(const libMesh::Point &x, const double &theta) {
  return libMesh::Point(x(0) * std::cos(theta) + x(1) * std::sin(theta),
                       -x(0) * std::sin(theta) + x(1) * std::cos(theta), 
                       0.0);
}

/*!
 * @brief Rotates a vector in xy-plane in anti-clockwise direction
 * @param x Point
 * @param theta Angle
 * @return Point after rotation
 */
inline std::vector<double> rotateACW2D(const std::vector<double> &x,
                                const double &theta) {
  return rotateCW2D(x, -theta);
}

/*!
 * @brief Rotates a vector in xy-plane in anti-clockwise direction
 * @param x Point
 * @param theta Angle
 * @return Point after rotation
 */
inline libMesh::Point rotateACW2D(const libMesh::Point &x, const double &theta) {
  return rotateCW2D(x, -theta);
}

/*!
 * @brief Rotates a vector in xy-plane assuming ACW convention
 * @param x Point
 * @param theta Angle
 * @return Point after rotation
 */
inline std::vector<double> rotate2D(const std::vector<double> &x, 
                                  const double &theta) {
  return std::vector<double>{x[0] * std::cos(theta) - x[1] * std::sin(theta),
                            x[0] * std::sin(theta) + x[1] * std::cos(theta),
                            0.0};
}

/*!
 * @brief Rotates a vector in xy-plane assuming ACW convention
 * @param x Point
 * @param theta Angle
 * @return Point after rotation
 */
inline libMesh::Point rotate2D(const libMesh::Point &x, const double &theta) {
  return libMesh::Point(x(0) * std::cos(theta) - x(1) * std::sin(theta),
                       x(0) * std::sin(theta) + x(1) * std::cos(theta), 
                       0.0);
}

/*!
 * @brief Computes derivative of rotation wrt to time
 *
 * If \f$ R(x,t) = Q(at)x \f$ then \f$ dR/dt = a Q' x \f$. This function returns \f$ Q' x \f$.
 *
 * @param x Point
 * @param theta Angle
 * @return Point after rotation
 */
inline libMesh::Point derRotate2D(const libMesh::Point &x, const double &theta) {
  return libMesh::Point(-x(0) * std::sin(theta) - x(1) * std::cos(theta),
                       x(0) * std::cos(theta) - x(1) * std::sin(theta), 
                       0.0);
}

/*!
 * @brief Returns the vector after rotating by desired angle
 *
 * @param p Vector
 * @param theta Angle of rotation
 * @param axis Axis of rotation
 * @return x Vector after rotation
 */
inline libMesh::Point rotate(const libMesh::Point &p, const double &theta, 
                           const libMesh::Point &axis) {
  auto ct = std::cos(theta);
  auto st = std::sin(theta);

  // dot
  double p_dot_n = p * axis;

  // cross
  libMesh::Point n_cross_p = axis.cross(p);

  return (1. - ct) * p_dot_n * axis + ct * p + st * n_cross_p;
}

/*!
 * @brief Computes angle between two vectors
 * @param a Vector 1
 * @param b Vector 2
 * @return angle Angle between vector a and b
 */
inline double angle(const libMesh::Point &a, const libMesh::Point &b) {
  if ((a - b).norm_sq() < 1.0E-12)
    return 0.;

  // since we do not know which side of plane given by normal
  // a x b / |a x b| is +ve, we compute the angle using cosine
  return std::acos((b * a) / (b.norm() * a.norm()));
}

/*!
 * @brief Computes angle between two vectors
 * @param a Vector 1
 * @param b Vector 2
 * @param axis Axis of rotation
 * @param is_axis If true then axis is the axis of orientation, otherwise
 * axis specifies the +ve side of the plane in which a and b are
 * @return angle Angle between vector a and b
 */
inline double angle(const libMesh::Point &a, const libMesh::Point &b, 
                   const libMesh::Point &axis, bool is_axis = true) {
  if ((a - b).norm_sq() < 1.0E-12)
    return 0.;

  if (is_axis) {
    // normal to plane of rotation
    libMesh::Point n = axis / axis.norm();
    libMesh::Point na = n.cross(a);

    double theta = std::atan((b * na) / (a * b - (b * n) * (a * n)));
    if (theta < 0.)
      theta += M_PI;

    if (b * na < 0.)
      theta = M_PI + theta;

    return theta;
  } else {
    auto theta = angle(a, b);

    // normal to a and b
    libMesh::Point n_ab = a.cross(b);

    double orient = axis * n_ab;
    if (orient < 0.)
      return 2. * M_PI - theta;
    else
      return theta;
  }
}

/*!
 * @brief Get vector with plus or minus 1 depending on the sign of component
 * of another vector
 *
 * @param v Vector from which a new vector is created
 * @return vec Vector
 */
inline libMesh::Point signVector(const libMesh::Point &v) {
  return libMesh::Point(v(0) > 0.0 ? 1.0 : -1.0, v(1) > 0.0 ? 1.0 : -1.0, v(2) > 0.0 ? 1.0 : -1.0);
}

inline double getDeterminant(const std::vector<libMesh::Point> &rows) {
  double a =
      rows[0](0) * (rows[1](1) * rows[2](2) - rows[1](2) * rows[2](1));
  double b =
      rows[0](1) * (rows[1](0) * rows[2](2) - rows[1](2) * rows[2](0));
  double c =
      rows[0](2) * (rows[1](0) * rows[2](1) - rows[1](1) * rows[2](0));

  return a - b + c;
}

/*!
 * @brief Returns all corner points in the box
 * @param dim Dimension of the box
 * @param box Pair of points representing cuboid (rectangle in 2d)
 * @return Vector Vector of corner points
 */
inline std::vector<libMesh::Point> getCornerPoints(
  const size_t &dim, const std::pair<libMesh::Point, libMesh::Point> &box) {
  if (dim == 1)
    return {box.first, box.second};
  else if (dim == 2)
    return {box.first, libMesh::Point(box.second(0), box.first(1), 0.),
            box.second, libMesh::Point(box.first(0), box.second(1), 0.)};
  else if (dim == 3) {
    double a = box.second(0) - box.first(0);
    double b = box.second(1) - box.first(1);
    double c = box.second(2) - box.first(2);
    return {box.first,
            box.first + libMesh::Point(a, 0., 0.),
            box.first + libMesh::Point(a, b, 0.),
            box.first + libMesh::Point(0., b, 0.),
            box.first + libMesh::Point(0., 0., c),
            box.first + libMesh::Point(a, 0., c),
            box.second,
            box.first + libMesh::Point(0., b, c)};
  } else {
    std::cerr << "Error: Check dimension = " << dim << ".\n";
    exit(1);
  }
}

/*!
 * @brief Checks if point is inside a rectangle
 * @param x Point
 * @param x_min X coordinate of left-bottom corner point
 * @param x_max X coordinate of right-top corner point
 * @param y_min Y coordinate of left-bottom corner point
 * @param y_max Y coordinate of right-top corner point
 * @return bool True if point inside rectangle, else false
 */
inline bool isPointInsideRectangle(const libMesh::Point &x, const double &x_min,
                                          const double &x_max, const double &y_min,
                                          const double &y_max) {
  return !(util::isLess(x(0), x_min - 1.0E-12) or
         util::isLess(x(1), y_min - 1.0E-12) or
         util::isGreater(x(0), x_max + 1.0E-12) or
         util::isGreater(x(1), y_max + 1.0E-12));
}

/*!
 * @brief Checks if point is inside a rectangle
 * @param x Point
 * @param x_lb Coordinate of left-bottom corner point
 * @param x_rt Coordinate of right-top corner point
 * @return bool True if point inside rectangle, else false
 */
inline bool isPointInsideRectangle(const libMesh::Point &x, const libMesh::Point &x_lb,
                                          const libMesh::Point &x_rt) {
  return !(util::isLess(x(0), x_lb(0) - 1.0E-12) or
         util::isLess(x(1), x_lb(1) - 1.0E-12) or
         util::isGreater(x(0), x_rt(0) + 1.0E-12) or
         util::isGreater(x(1), x_rt(1) + 1.0E-12));
}

/*!
 * @brief Checks if point is inside an angled rectangle
 * @param x Point
 * @param x_min X coordinate of left-bottom corner point
 * @param x_max X coordinate of right-top corner point
 * @param y_min Y coordinate of left-bottom corner point
 * @param y_max Y coordinate of right-top corner point
 * @param theta Angle of orientation of rectangle from x-axis
 * @return bool True if point inside rectangle, else false
 */
inline bool isPointInsideAngledRectangle(const libMesh::Point &x, const double &x1,
  const double &x2, const double &y1,
  const double &y2, const double &theta) {
  // we assume that the rectangle has passed the test

  //
  //                             (x2,y2)
  //                            o
  //
  //
  //
  //
  //
  //        o
  //      (x1,y1)

  // get divisors
  libMesh::Point lam = rotateCW2D(libMesh::Point(x2 - x1, y2 - y1, 0.0), theta);

  // double lam1 = (x2-x1) * std::cos(theta) + (y2-y1) * std::sin(theta);
  // double lam2 = -(x2-x1) * std::sin(theta) + (y2-y1) * std::cos(theta);

  // get mapped coordinate of x
  libMesh::Point xmap = rotateCW2D(libMesh::Point(x(0) - x1, x(1) - y1, 0.0), theta);

  // double xmap = (x[0]-x1) * std::cos(theta) + (x[1]-y1) * std::sin(theta);
  // double ymap = -(x[0]-x1) * std::sin(theta) + (x[1]-y1) * std::cos(theta);

  // check if mapped coordinate are out of range [0, lam1] and [0, lam2]
  return !(util::isLess(xmap(0), -1.0E-12) or
      util::isLess(xmap(1), -1.0E-12) or
      util::isGreater(xmap(0), lam(0) + 1.0E-12) or
      util::isGreater(xmap(1), lam(1) + 1.0E-12));
}

/*!
 * @brief Checks if point is inside a cuboid (rectangle in 2-d, line in 1-d)
 * @param dim Dimension
 * @param x Point
 * @param x_lbb Coordinate of left-bottom-back corner point
 * @param x_rtf Coordinate of right-top-front corner point
 * @return bool True if point inside cuboid, else false
 */
inline bool isPointInsideCuboid(const size_t &dim, const libMesh::Point &x,
    const libMesh::Point &x_lbb,
    const libMesh::Point &x_rtf) {
  if (dim == 1)
    return !(isLess(x(0), x_lbb(0) - 1.0E-12) or isGreater(x(0), x_rtf(0) + 1.0E-12));
  else if (dim == 2)
    return !(isLess(x(0), x_lbb(0) - 1.0E-12) or
        isLess(x(1), x_lbb(1) - 1.0E-12) or
        isGreater(x(0), x_rtf(0) + 1.0E-12) or
        isGreater(x(1), x_rtf(1) + 1.0E-12));
  else
    return !(isLess(x(0), x_lbb(0) - 1.0E-12) or
          isLess(x(1), x_lbb(1) - 1.0E-12) or
          isLess(x(2), x_lbb(2) - 1.0E-12) or
          isLess(x(2), x_lbb(2) - 1.0E-12) or
          isGreater(x(0), x_rtf(0) + 1.0E-12) or
          isGreater(x(1), x_rtf(1) + 1.0E-12) or
          isGreater(x(2), x_rtf(2) + 1.0E-12));
}

/*!
 * @brief Computes the area of triangle
 * @param nodes Vertices of the triangle
 * @return area Area of triangle
 */
inline double getTriangleArea(const std::vector<libMesh::Point> &nodes) {
  return 0.5 * ((nodes[1](0) - nodes[0](0)) * (nodes[2](1) - nodes[0](1)) -
                (nodes[2](0) - nodes[0](0)) * (nodes[1](1) - nodes[0](1)));
}

/*!
 * @brief Computes the volume of tetrahedron
 * @param nodes Vertices of the tetrahedron
 * @return area Area of tetrahedron
 */
inline double getTetVolume(const std::vector<libMesh::Point> &nodes) {
  auto a = nodes[1] - nodes[0];
  auto b = nodes[2] - nodes[0];
  auto c = nodes[3] - nodes[0];

  return (1. / 6.) * getDeterminant({a, b, c});
}

/*!
 * @brief Checks if the point C in on the line between A and B 
 * @param A Point 1
 * @param B Point 2
 * @param C Point 3 
 * @return Direction
 */
inline int direction(const libMesh::Point &A , const libMesh::Point &B, const libMesh::Point &C) {
  int val =
      (B(1) - A(1)) * (C(0) - B(0)) - (B(0) - A(0)) * (C(1) - B(1));
  if (val == 0)
    return 0;  // colinear
  else if (val < 0)
    return 2;  // anti-clockwise direction
  return 1;    // clockwise direction
}

inline bool onLine(const libMesh::Point &A, const libMesh::Point &B,
  const libMesh::Point &C) {
  if (C(0) <= std::max(A(0), B(0)) && C(0) <= std::min(A(0), B(0)) &&
      (C(1) <= std::max(A(1), B(1)) && C(1) <= std::min(A(1), B(1))))
    return true;

  return false;
}

/*!
 * @brief Checks if the two lines between A and B and between C and C intersect
 * @param A Start point of the first line
 * @param B End point of the first line
 * @param C Start point of the second line
 * @param D End point of the second line
 * @return intersection 
 */
inline bool doLinesIntersect (const libMesh::Point &A , const libMesh::Point &B, const libMesh::Point &C, const libMesh::Point & D) {
  // four direction for two lines and points of other line
  int dir1 = direction(A, B, C);
  int dir2 = direction(A, B, D);
  int dir3 = direction(C, D, A);
  int dir4 = direction(C, D, B);

  if (dir1 != dir2 && dir3 != dir4) return true;  // they are intersecting

  if (dir1 == 0 && onLine(A, B, C))  // when p2 of line2 are on the line1
    return true;

  if (dir2 == 0 && onLine(A, B, D))  // when p1 of line2 are on the line1
    return true;

  if (dir3 == 0 && onLine(C, D, A))  // when p2 of line1 are on the line2
    return true;

  if (dir4 == 0 && onLine(C, D, B))  // when p1 of line1 are on the line2
    return true;

  return false;
 }

/** @}*/

} // namespace util
