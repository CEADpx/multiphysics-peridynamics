/*
 * -------------------------------------------
 * Copyright (c) 2021 - 2024 Prashant K. Jha
 * -------------------------------------------
 * PeriDEM https://github.com/prashjha/PeriDEM
 *
 * Distributed under the Boost Software License, Version 1.0. (See accompanying
 * file LICENSE)
 */

#ifndef MATERIAL_PD_INFLUENCEFN_H
#define MATERIAL_PD_INFLUENCEFN_H

#include "io.h"
#include <cstring>
#include <iostream>
#include <vector>

namespace material {

/*! @brief A base class for computing influence function */
class BaseInfluenceFn {

public:
  /*! @brief Constructor */
  BaseInfluenceFn() = default;

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  virtual double get(const double &r) const = 0;

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th} \f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return moment Moment
   */
  virtual double getMoment(const size_t &i) const = 0;

  /*!
   * @brief Returns the string containing printable information about the object
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   * @return string String containing printable information about the object
   */
  virtual std::string printStr(int nt = 0, int lvl = 0) const {

    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- material::BaseInfluenceFn --------" << std::endl << std::endl;
    oss << tabS << "Provides abstraction for different influence function "
                   "types" << std::endl;
    oss << tabS << std::endl;

    return oss.str();
  }

  /*!
   * @brief Prints the information about the object
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   */
  virtual void print(int nt = 0, int lvl = 0) const { std::cout << printStr(nt, lvl); }

};

/*! @brief A class to implement constant influence function */
class ConstInfluenceFn : public BaseInfluenceFn {

public:
  /*!
   * @brief Constructor
   * @param params List of parameters
   * @param dim Dimension
   */
  ConstInfluenceFn(const std::vector<double> &params, const size_t &dim) {
    d_a0 = params.empty() ? double(dim + 1) : params[0];
  };

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  double get(const double &r) const override { return d_a0; }

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th}\f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return moment Moment
   */
  double getMoment(const size_t &i) const override { return d_a0 / double(i + 1); }

  /*!
   * @brief Returns the string containing printable information about the object
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   * @return string String containing printable information about the object
   */
  std::string printStr(int nt = 0, int lvl = 0) const override {

    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- material::ConstInfluenceFn --------" << std::endl << std::endl;
    oss << tabS << "Constant function with constant = " << d_a0 << std::endl;
    oss << tabS << "First moment = " << getMoment(1) << std::endl;
    oss << tabS << "Second moment = " << getMoment(2) << std::endl;
    oss << tabS << "Third moment = " << getMoment(3) << std::endl;
    oss << tabS << std::endl;

    return oss.str();
  }

  /*!
   * @brief Prints the information about the object
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   */
  void print(int nt = 0, int lvl = 0) const override {
    std::cout << printStr(nt, lvl);
  }

private:
  /*! @brief Constant such that J(r) = Constant */
  double d_a0;
};

/*! @brief A class to implement linear influence function
 *
 * \f$ J(r) = a0 + a1 r \f$
 */
class LinearInfluenceFn : public BaseInfluenceFn {

public:
  /*!
   * @brief Constructor
   * @param params List of parameters
   * @param dim Dimension
   */
  LinearInfluenceFn(const std::vector<double> &params, const size_t &dim)
  : BaseInfluenceFn(), d_a0(0.), d_a1(0.) {

    if (params.empty()) {
      // choose a0, a1 = -a0 such that \int_0^1 J(r) r^d dr = 1
      // and J(r) = a0 (1 - r)
      if (dim == 1) {
        d_a0 = 6.;
        d_a1 = -d_a0;
      } else if (dim == 2) {
        d_a0 = 12.;
        d_a1 = -d_a0;
      } else if (dim == 3) {
        d_a0 = 20.;
        d_a1 = -d_a0;
      }
    } else {
      d_a0 = params[0];
      if (params.size() < 2)
        d_a1 = -d_a0;
      else
        d_a1 = params[1];
    }
  };

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  double get(const double &r) const override { return d_a0 + d_a1 * r; }

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th}\f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return moment Moment
   */
  double getMoment(const size_t &i) const override {
    return (d_a0 / double(i + 1)) + (d_a1 / double(i + 2));
  }

  /*!
   * @brief Returns the string containing printable information about the object
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   * @return string String containing printable information about the object
   */
  std::string printStr(int nt = 0, int lvl = 0) const override {

    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- material::LinearInfluenceFn --------" << std::endl << std::endl;
    oss << tabS << "Linear function a0 + a1*r with constants: a0 = "
                << d_a0 << ", a1 = " << d_a1 << std::endl;
    oss << tabS << "First moment = " << getMoment(1) << std::endl;
    oss << tabS << "Second moment = " << getMoment(2) << std::endl;
    oss << tabS << "Third moment = " << getMoment(3) << std::endl;
    oss << tabS << std::endl;

    return oss.str();
  }

  /*!
   * @brief Prints the information about the object
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   */
  void print(int nt = 0, int lvl = 0) const override {
    std::cout << printStr(nt, lvl);
  }

private:
  /*! @brief Constants such that J(r) = d_a0 + d_a1 * r */
  double d_a0;

  /*! @brief Constants such that J(r) = d_a0 + d_a1 * r */
  double d_a1;
};

/*! @brief A class to implement Gaussian influence function
 *
 * \f$ J(r) = \alpha \exp(-r^2/\beta) \f$
 */
class GaussianInfluenceFn : public BaseInfluenceFn {

public:
  /*!
   * @brief Constructor
   * @param params List of parameters
   * @param dim Dimension
   */
  GaussianInfluenceFn(const std::vector<double> &params, const size_t &dim)
  : BaseInfluenceFn(), d_alpha(0.), d_beta(0.) {

    if (params.empty()) {
      // beta = 0.2 (default value)
      // choose alpha such that \int_0^1 J(r) r^d dr = 1
      d_beta = 0.2;
      if (dim == 1)
        d_alpha = 2. / (d_beta * (1. - std::exp(-1. / d_beta)));
      else if (dim == 2)
        d_alpha = (4.0 / d_beta) * 1.0 /
                  (std::sqrt(M_PI * d_beta) * std::erf(1.0 / std::sqrt(d_beta)) -
                  2.0 * std::exp(-1.0 / d_beta));
      else if (dim == 3)
        d_alpha = (2.0 / d_beta) * 1.0 /
                  (d_beta - (d_beta + 1.) * std::exp(-1.0 / d_beta));
    } else {
      d_alpha = params[0];
      d_beta = params[1];
    }
  };

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  double get(const double &r) const override {
    return d_alpha * std::exp(-r * r / d_beta);
  }

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th}\f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return moment Moment
   */
  double getMoment(const size_t &i) const override  {

    double sq1 = std::sqrt(d_beta);
    double sq2 = std::sqrt(M_PI);
    // M_i = \int_0^1 alpha exp(-r^2/beta) r^i dr
  
    if (i == 0) {
      // M0 = 0.5 * \alpha (\beta)^(1/2) * (pi)^(1/2) * erf((1/beta)^(1/2))
  
      return 0.5 * d_alpha * sq1 * sq2 * std::erf(1. / sq1);
    } else if (i == 1) {
      // M1 = 0.5 * \alpha \beta (1 - exp(-1/beta))
  
      return 0.5 * d_alpha * d_beta * (1. - std::exp(-1. / d_beta));
    } else if (i == 2) {
      // M2 = 0.5 * \alpha (\beta)^(3/2) * [0.5 * (pi)^(1/2) erf((1/beta)^(1/2)
      // ) - (1/beta)^(1/2) * exp(-1/beta) ]
  
      return 0.5 * d_alpha * d_beta * sq1 *
             (0.5 * sq2 * std::erf(1. / sq1) -
              (1. / sq1) * std::exp(-1. / d_beta));
    } else if (i == 3) {
      // M3 = 0.5 * \alpha (\beta)^(2) * [1 - ((1/beta) + 1) * exp(-1/beta)]
  
      return 0.5 * d_alpha * d_beta * d_beta *
             (1. - (1. + 1. / d_beta) * std::exp(-1. / d_beta));
    } else {
      std::cerr << "Error: getMoment() accepts argument i from 0 to 3.\n";
      exit(1);
    }
  
  
  };

  /*!
   * @brief Returns the string containing printable information about the object
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   * @return string String containing printable information about the object
   */
  std::string printStr(int nt = 0, int lvl = 0) const override {

    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- material::GaussianInfluenceFn --------" << std::endl << std::endl;
    oss << tabS << "Gaussian function a0 * exp(-r*r / a1) with constants: a0 = "
                << d_alpha << ", a1 = " << d_beta << std::endl;
    oss << tabS << "First moment = " << getMoment(1) << std::endl;
    oss << tabS << "Second moment = " << getMoment(2) << std::endl;
    oss << tabS << "Third moment = " << getMoment(3) << std::endl;
    oss << tabS << std::endl;

    return oss.str();
  }

  /*!
   * @brief Prints the information about the object
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   */
  void print(int nt = 0, int lvl = 0) const override {
    std::cout << printStr(nt, lvl);
  }

private:
  /*! @brief Constants */
  double d_alpha;

  /*! @brief Constants */
  double d_beta;
};

} // namespace material

#endif // MATERIAL_PD_INFLUENCEFN_H
