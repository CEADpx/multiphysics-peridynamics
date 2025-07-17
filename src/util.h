/*
 * -------------------------------------------
 * Copyright (c) 2021 - 2024 Prashant K. Jha
 * -------------------------------------------
 * PeriDEM https://github.com/prashjha/PeriDEM
 *
 * Distributed under the Boost Software License, Version 1.0. (See accompanying
 * file LICENSE)
 */

#ifndef UTIL_FUNCTION_H
#define UTIL_FUNCTION_H

#include <libmesh/point.h>
#include "io.h"
#include "json.h"
#include "geometry.h"
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

/** @}*/

} // namespace util

#endif // UTIL_FUNCTION_H
