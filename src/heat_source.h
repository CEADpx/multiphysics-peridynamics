#ifndef HEAT_SOURCE_H
#define HEAT_SOURCE_H

#include "geometry.h"
#include "libmesh/point.h"
#include <memory>
#include "util.h"

namespace loading{
/**
 * @brief Abstract base class for heat sources in the thermomechanical model
 * 
 * This class defines the interface for heat sources. Derived classes must implement
 * the get() method to provide the volumetric heat generation rate at a given
 * position and time.
 */
class BaseHeatSource {
public:
  /**
   * @brief Virtual destructor
   */
  virtual ~BaseHeatSource() = default;

  /**
   * @brief Get heat source value at given position and time
   * @param x Position vector
   * @param time Current time
   * @return Heat source value in W/m³
   */
  virtual double get(const libMesh::Point& x, const double& time) const {
    return 0.0;  // Base class returns zero heat source
  }

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
    oss << tabS << "------- loading::BaseHeatSource --------" << std::endl
        << std::endl;
    oss << tabS << std::endl;

    return oss.str();
  };

  /*!
   * @brief Prints the information about the object
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   */
  virtual void print(int nt = 0, int lvl = 0) const {
    std::cout << printStr(nt, lvl);
  };
};

/**
 * @brief Example implementation of a Gaussian heat source
 */
class GaussianHeatSource : public BaseHeatSource {
public:
  /**
   * @brief Constructor
   * @param center Center position of the Gaussian
   * @param sigma Standard deviation of the Gaussian
   * @param amplitude Maximum heat source intensity (W/m³)
   */
  GaussianHeatSource(const libMesh::Point& center, const double& sigma, const double& amplitude)
    : d_center(center), d_sigma(sigma), d_amplitude(amplitude) {}


  double get(const libMesh::Point& x, const double& time) const override {
    double r = (x - d_center).norm();
    return d_amplitude * std::exp(-r*r/(2.0*d_sigma*d_sigma));
  }

  libMesh::Point d_center;   ///< Center of the Gaussian
  double d_sigma;            ///< Standard deviation
  double d_amplitude;        ///< Maximum intensity

  std::string printStr(int nt = 0, int lvl = 0) const override {
    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- loading::GaussianHeatSource --------" << std::endl
        << std::endl;
    oss << tabS << "Center: " << d_center << std::endl;
    oss << tabS << "Sigma: " << d_sigma << std::endl;
    oss << tabS << "Amplitude: " << d_amplitude << std::endl;
    oss << tabS << std::endl;
    return oss.str();
  }

  void print(int nt = 0, int lvl = 0) const override {
    std::cout << printStr(nt, lvl);
  }
};

/**
 * @brief Example implementation of a time-dependent pulsed heat source
 */
class PulsedHeatSource : public BaseHeatSource {
public:
  /**
   * @brief Constructor
   * @param center Center position of the heat source
   * @param radius Radius of influence
   * @param amplitude Maximum heat source intensity (W/m³)
   * @param frequency Pulse frequency (Hz)
   */
  PulsedHeatSource(const libMesh::Point& center, const double& radius, 
                   const double& amplitude, const double& frequency)
    : d_center(center), d_radius(radius), d_amplitude(amplitude), d_frequency(frequency) {}


  double get(const libMesh::Point& x, const double& time) const override {
    double r = (x - d_center).norm();
    if (r > d_radius) return 0.0;
    
    // Sinusoidal pulse
    return d_amplitude * std::sin(2.0 * M_PI * d_frequency * time) * 
           (1.0 - r/d_radius);  // Linear decay with distance
  }

  libMesh::Point d_center;   ///< Center of the heat source
  double d_radius;           ///< Radius of influence
  double d_amplitude;        ///< Maximum intensity
  double d_frequency;        ///< Pulse frequency

  std::string printStr(int nt = 0, int lvl = 0) const override {
    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- loading::PulsedHeatSource --------" << std::endl
        << std::endl;
    oss << tabS << "Center: " << d_center << std::endl;
    oss << tabS << "Radius: " << d_radius << std::endl;
    oss << tabS << "Amplitude: " << d_amplitude << std::endl;
    oss << tabS << "Frequency: " << d_frequency << std::endl;
    oss << tabS << std::endl;
    return oss.str();
  }

  void print(int nt = 0, int lvl = 0) const override {
    std::cout << printStr(nt, lvl);
  } 
};

/**
 * @brief Example implementation of a Gaussian heat source
 */
 class HeatSource : public BaseHeatSource {
  public:

     HeatSource(std::string spatial_fn_type = "", 
          std::vector<double> spatial_fn_params = {}, 
          std::string time_fn_type = "", 
          std::vector<double> time_fn_params = {1},
          std::shared_ptr<geom::GeomObject> spatial_geom_p = nullptr)
      : d_spatial_fn_type(spatial_fn_type), 
        d_spatial_fn_params(spatial_fn_params), 
        d_time_fn_type(time_fn_type), 
        d_time_fn_params(time_fn_params),
        d_spatial_geom_p(std::move(spatial_geom_p)) {}
  
    double get(const libMesh::Point& x, const double& time) const override {
      double fmax = 1.0;

      bool is_inside = d_spatial_geom_p == nullptr ? true : d_spatial_geom_p->isInside(x);
      if (!is_inside) return 0.0;

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
      } else if (d_spatial_fn_type == "gaussian") {
        double sigma = d_spatial_fn_params[0];
        double beta = 2.0*sigma*sigma;
        double a = d_spatial_fn_params[1];
        libMesh::Point center(0, 0, 0);
        if (d_spatial_fn_params.size() > 2) {
          for (size_t i = 2; i < d_spatial_fn_params.size(); ++i) {
            center(i-2) = d_spatial_fn_params[i];
          }
        }
        double r = (x - center).norm();
        fmax = util::gaussian(r, a, beta);
      }

      // apply time function
      if (d_time_fn_type == "linear")
        fmax *= time;
      else if (d_time_fn_type == "linear_step")
        fmax *= util::linearStepFunc(time, d_time_fn_params[1],
                                      d_time_fn_params[2]);
      else if (d_time_fn_type == "linear_step_const_value")
        fmax *= util::linearStepSpecifiedConstFunc(time, d_time_fn_params[2],
                                    d_time_fn_params[3], d_time_fn_params[1]);
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

      return fmax;
    }
  
    std::string d_spatial_fn_type;   ///< Type of spatial function
    std::vector<double> d_spatial_fn_params; ///< Parameters for spatial function
    std::string d_time_fn_type;   ///< Type of time function
    std::vector<double> d_time_fn_params; ///< Parameters for time function
    std::shared_ptr<geom::GeomObject> d_spatial_geom_p;
  
    std::string printStr(int nt = 0, int lvl = 0) const override {
      auto tabS = util::io::getTabS(nt);
      std::ostringstream oss;
      oss << tabS << "------- loading::HeatSource --------" << std::endl
          << std::endl;
      oss << tabS << "Spatial function type: " << util::io::printStr(d_spatial_fn_type) << std::endl;
      oss << tabS << "Spatial function parameters: " << util::io::printStr(d_spatial_fn_params) << std::endl;
      oss << tabS << "Time function type: " << util::io::printStr(d_time_fn_type) << std::endl;
      oss << tabS << "Time function parameters: " << util::io::printStr(d_time_fn_params) << std::endl;
      oss << tabS << std::endl;
      return oss.str();
    }
  
    void print(int nt = 0, int lvl = 0) const override {
      std::cout << printStr(nt, lvl);
    }
  };

  class HeatSourceCollection {
    public:
      std::vector<HeatSource> d_heat_sources;
      
      HeatSourceCollection() {};

      void addHeatSource(HeatSource heat_source) {
        d_heat_sources.push_back(heat_source);
      }
    
      /**
       * @brief Get heat source value at given position and time
       * @param x Position vector
       * @param time Current time
       * @return Heat source value in W/m³
       */
      double get(const libMesh::Point& x, const double& time) const {
        double fmax = 0.0;
        for (const auto& heat_source : d_heat_sources) {
          fmax += heat_source.get(x, time);
        }
        return fmax;
      }
    
        /*!
       * @brief Returns the string containing printable information about the object
       *
       * @param nt Number of tabs to append before printing
       * @param lvl Information level (higher means more information)
       * @return string String containing printable information about the object
       */
      std::string printStr(int nt = 0, int lvl = 0) const {
    
        auto tabS = util::io::getTabS(nt);
        std::ostringstream oss;
        oss << tabS << "------- loading::HeatSourceCollection --------" << std::endl
            << std::endl;
        oss << tabS << "Number of heat sources: " << d_heat_sources.size() << std::endl;
        for (size_t i = 0; i < d_heat_sources.size(); i++) {
          oss << tabS << "Heat source " << i << ": " << std::endl;
          oss << d_heat_sources[i].printStr(nt + 1, lvl) << std::endl;
        }
        oss << tabS << std::endl;
    
        return oss.str();
      };
    
      /*!
       * @brief Prints the information about the object
       *
       * @param nt Number of tabs to append before printing
       * @param lvl Information level (higher means more information)
       */
      void print(int nt = 0, int lvl = 0) const {
        std::cout << printStr(nt, lvl);
      };
    };

} // namespace loading

#endif // HEAT_SOURCE_H 