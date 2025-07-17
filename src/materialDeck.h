#ifndef MATERIAL_DECK_H
#define MATERIAL_DECK_H

#include "io.h"

namespace inp {

enum class ModelType { StateBased, BondBased };

struct MaterialDeck {
  
  // Mechanical properties
  ModelType d_model_type;
  double d_E;         // Young's modulus (Pa)
  double d_K;         // Bulk modulus (Pa)
  double d_Kadjust;         // Bulk modulus (Pa)
  double d_G;         // Shear modulus (Pa)
  double d_nu;        // Poisson's ratio
  double d_rho;       // Density (kg/m^3)
  double d_Gc;         // Critical energy release rate (J/m^2)

  // Plane stress/strain flag
  bool d_isPlaneStrain;

  // Problem dimension
  int d_dim;

  // Peridynamic horizon
  double d_horizon; // (m)

  // Thermal properties
  double d_Ktherm;    // Thermal conductivity (W/m/K)
  double d_Tref;      // Reference temperature (K)
  double d_Cv;        // Specific heat (J/kg/K)
  double d_alpha;     // Thermal expansion coefficient (1/K)
  bool d_robinBC;     // Flag to use Robin BC
  double d_hconvect;  // Convective heat transfer coefficient (W/m^2/K)

  // Influence function type: 0 = Constant, 1 = Linear, 2 = Gaussian
  int d_influenceFnType = 1;

  // Influence function parameters
  std::vector<double> d_influenceFnParams;

  // constructor
  MaterialDeck() {}

  void setDefaults() {
    std::string s = "default";

    if (s == "default") {
      d_model_type = ModelType::StateBased;
      d_E = 1.0;
      d_nu = 0.25;
      d_K = d_E/(3*(1-2*d_nu));
      d_Kadjust = d_K;
      d_G = d_E/(2*(1+d_nu)); 
      d_rho = 1.0;
      d_Gc = 1.0;

      d_isPlaneStrain = false;

      d_dim = 3;

      d_horizon = 1.0;

      d_Ktherm = 0.1;
      d_Tref = 273.0;
      d_Cv = 10.0;
      d_alpha = 1e-6;
      d_robinBC = true;
      d_hconvect = 50.;

      d_influenceFnType = 2;
      d_influenceFnParams = {1.0, 0.25};
    }

    if (s == "granite") { 
      d_model_type = ModelType::StateBased;
      d_E = 67.0e9;
      d_nu = 0.33;
      d_K = d_E/(3*(1-2*d_nu));
      d_Kadjust = d_K;
      d_G = d_E/(2*(1+d_nu)); 
      d_rho = 2650.0;
      d_Gc = 70.0;

      d_isPlaneStrain = false;

      d_dim = 3;

      d_horizon = 1.0;

      d_Ktherm = 3.5;
      d_Tref = 293.0;
      d_Cv = 1015.0;
      d_alpha = 3.5e-6;
      d_robinBC = true;
      d_hconvect = 100.;

      d_influenceFnType = 1;
      d_influenceFnParams = {1.0, -1.0};
    }
  }

  std::string printStr(int nt = 0, int lvl = 0) const {
    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- inp::MaterialDeck --------" << std::endl
        << std::endl;
    oss << tabS << "E = " << d_E << std::endl;
    oss << tabS << "model_type = " << (d_model_type == ModelType::StateBased ? "StateBased" : "BondBased") << std::endl;
    oss << tabS << "K = " << d_K << std::endl;
    oss << tabS << "Kadjust = " << d_Kadjust << std::endl;
    oss << tabS << "G = " << d_G << std::endl;
    oss << tabS << "nu = " << d_nu << std::endl;
    oss << tabS << "rho = " << d_rho << std::endl;
    oss << tabS << "Gc = " << d_Gc << std::endl;
    oss << tabS << "isPlaneStrain = " << d_isPlaneStrain << std::endl;
    oss << tabS << "dim = " << d_dim << std::endl;
    oss << tabS << "horizon = " << d_horizon << std::endl;
    oss << tabS << "Ktherm = " << d_Ktherm << std::endl;
    oss << tabS << "Tref = " << d_Tref << std::endl;
    oss << tabS << "Cv = " << d_Cv << std::endl;
    oss << tabS << "alpha = " << d_alpha << std::endl;
    oss << tabS << "robinBC = " << d_robinBC << std::endl;
    oss << tabS << "hconvective = " << d_hconvect << std::endl;
    oss << tabS << "influenceFnType = " << d_influenceFnType << std::endl;
    oss << tabS << "influenceFnParams = ";
    for (auto param : d_influenceFnParams) oss << param << " ";
    oss << std::endl; 
    oss << tabS << std::endl;
    return oss.str();
  }

  void print(int nt = 0, int lvl = 0) const {
    std::cout << printStr(nt, lvl);
  }
};

} // namespace inp

#endif // MATERIAL_DECK_H