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

#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include "geometry.h"
#include "util.h"
#include <libmesh/point.h>

namespace inp {

/**
 * @brief Simple boundary condition class for single domain simulation
 */
struct BCBase {
  // Type of BC (Force_BC, Displacement_BC, IC)
  std::string d_type;

  // Region for applying BC
  std::shared_ptr<geom::GeomObject> d_region_p;

  // Time function parameters
    /*!
    * @brief Name of the formula with respect to time
    *
    * List of allowed values are:
    * - "" (none)
    * - \a constant
    * - \a linear
    * - \a linear_step
    * - \a linear_slow_fast
    * - \a rotation
    */
  std::string d_time_fn_type;  // constant, linear
  std::vector<double> d_time_fn_params;

      /*!
     * @brief Name of the formula of with respect to spatial coordinate
     *
     * List of allowed values are:
     * - "" (none)
     * - \a constant
     * - \a hat_x
     * - \a hat_y
     * - \a sin
     * - \a rotation
     */
     std::string d_spatial_fn_type;
     std::vector<double> d_spatial_fn_params;

  // Direction to apply BC (x=0, y=1, z=2)
  std::vector<size_t> d_direction;

  // For displacement BC
  bool d_is_zero;

  // For initial conditions
  libMesh::Point d_ic_vector;

  // Constructor
  BCBase(const std::string& bcType = "") 
    : d_type(bcType), d_is_zero(false) {}

  // Print information
  std::string printStr(int nt = 0, int lvl = 0) const {
    std::ostringstream oss;
    auto tabS = util::io::getTabS(nt);
    oss << tabS << "------- inp::BCBase --------" << std::endl
        << std::endl;
    oss << tabS << "BC Type: " << d_type << std::endl;
    if (d_region_p)
      oss << tabS << "Region: " << d_region_p->getName() << std::endl;
      oss << tabS<< "Region data: " << std::endl;
      oss << d_region_p->printStr(nt + 1, lvl) << std::endl;
    if (!d_time_fn_type.empty()) {
      oss << tabS << "Time Function: " << d_time_fn_type << std::endl;
      oss << tabS << "Time Parameters: ";
      for (auto p : d_time_fn_params) oss << p << " ";
      oss << std::endl;
    }
    if (!d_spatial_fn_type.empty()) {
      oss << tabS << "Spatial Function: " << d_spatial_fn_type << std::endl;
      oss << tabS << "Spatial Parameters: ";
      for (auto p : d_spatial_fn_params) oss << p << " ";
      oss << std::endl;
    }
    if (!d_direction.empty()) {
      oss << tabS << "Direction: ";
      for (auto d : d_direction) oss << d << " ";
      oss << std::endl;
    }
    if (d_is_zero)
      oss << tabS << "Zero Displacement" << std::endl;
    if (d_type == "IC")
      oss << tabS << "IC Vector: (" << d_ic_vector(0) << ", " 
          << d_ic_vector(1) << ", " << d_ic_vector(2) << ")" << std::endl;
    oss << tabS << std::endl;
    return oss.str();
  }

  void print(int nt = 0, int lvl = 0) const { std::cout << printStr(nt, lvl); }

  /**
    * @brief Check if a point is in the specified region
    * @param x Point to check
    * @param bc Boundary condition containing region info
    * @return true if point is in region
    */
    bool isInRegion(const libMesh::Point &x) const {
        if (!d_region_p)
          return true;  // no region specified means apply everywhere
        return d_region_p->isInside(x);
    }

    void getBCDispVel(const double time, const libMesh::Point &x, libMesh::Point &ux, libMesh::Point &vx) const {
        
                double umax = d_time_fn_params[0];
                double du = 0.;
                double dv = 0.;

    
                  if (d_spatial_fn_type == "sin_x") {
                    double a = M_PI * d_spatial_fn_params[0];
                    umax = umax * std::sin(a * x(0));
                  } else if (d_spatial_fn_type == "sin_y") {
                    double a = M_PI * d_spatial_fn_params[0];
                    umax = umax * std::sin(a * x(1));
                  } else if (d_spatial_fn_type == "sin_z") {
                    double a = M_PI * d_spatial_fn_params[0];
                    umax = umax * std::sin(a * x(2));
                  } else if (d_spatial_fn_type == "linear_x") {
                    double a = d_spatial_fn_params[0];
                    umax = umax * a * x(0);
                  } else if (d_spatial_fn_type == "linear_y") {
                    double a = d_spatial_fn_params[0];
                    umax = umax * a * x(1);
                  } else if (d_spatial_fn_type == "linear_z") {
                    double a = d_spatial_fn_params[0];
                    umax = umax * a * x(2);
                  }

                  // apply time function
                  if (d_time_fn_type == "constant")
                    du = umax;
                  else if (d_time_fn_type == "linear") {
                    du = umax * time;
                    dv = umax;
                  } else if (d_time_fn_type == "quadratic") {
                    du = umax * time + d_time_fn_params[1] * time * time;
                    dv = umax + d_time_fn_params[1] * time;
                  } else if (d_time_fn_type == "sin") {
                    double a = M_PI * d_time_fn_params[1];
                    du = umax * std::sin(a * time);
                    dv = umax * a * std::cos(a * time);
                  }

                  ux = libMesh::Point();
                  vx = libMesh::Point();
                  for (auto d : d_direction) {
                    ux(d) = du;
                    vx(d) = dv;
                  }

                  if (d_time_fn_type == "rotation") {
                    auto x0 = libMesh::Point(d_time_fn_params[1], d_time_fn_params[2],
                                          d_time_fn_params[3]);
                    auto dx = x - x0;
                    auto r_x = util::rotate2D(
                            dx, d_time_fn_params[0] * time);
                    auto dr_x = util::derRotate2D(
                            dx, d_time_fn_params[0] * time);

                    ux += r_x - dx;
                    vx += d_time_fn_params[0] * dr_x;
                  }
    }

    void getBCForce(const double time, const libMesh::Point &x, libMesh::Point &fx) const {
                double fmax = 1.0;


                  if (d_spatial_fn_type == "sin_x") {
                    double a = M_PI * d_spatial_fn_params[0];
                    fmax = d_spatial_fn_params[0] * std::sin(a * x(0));
                  } else if (d_spatial_fn_type == "sin_y") {
                    double a = M_PI * d_spatial_fn_params[0];
                    fmax = d_spatial_fn_params[0] * std::sin(a * x(1));
                  } else if (d_spatial_fn_type == "sin_z") {
                    double a = M_PI * d_spatial_fn_params[0];
                    fmax = d_spatial_fn_params[0] * std::sin(a * x(2));
                  } else if (d_spatial_fn_type == "linear_x") {
                    double a = d_spatial_fn_params[0];
                    fmax = d_spatial_fn_params[0] * a * x(0);
                  } else if (d_spatial_fn_type == "linear_y") {
                    double a = d_spatial_fn_params[0];
                    fmax = d_spatial_fn_params[0] * a * x(1);
                  } else if (d_spatial_fn_type == "linear_z") {
                    double a = d_spatial_fn_params[0];
                    fmax = d_spatial_fn_params[0] * a * x(2);
                  }

                  // apply time function
                  if (d_time_fn_type == "linear")
                    fmax *= time;
                  else if (d_time_fn_type == "linear_step")
                    fmax *= util::linearStepFunc(time, d_time_fn_params[1],
                                                 d_time_fn_params[2]);
                  else if (d_time_fn_type == "linear_slow_fast") {
                    if (util::isGreater(time, d_time_fn_params[1]))
                      fmax *= d_time_fn_params[3] * time;
                    else
                      fmax *= d_time_fn_params[2] * time;
                  } else if (d_time_fn_type == "sin") {
                    double a = M_PI * d_time_fn_params[1];
                    fmax *= std::sin(a * time);
                  }

                  // multiply by the slope
                  fmax *= d_time_fn_params[0];

                  fx = libMesh::Point();
                  for (auto d : d_direction) {
                    fx(d) = fmax;
                  }
              }
};

} // namespace inp 