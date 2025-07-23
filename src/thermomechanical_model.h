/*
 * -------------------------------------------
 * Copyright (c) 2021 - 2024 Prashant K. Jha
 * -------------------------------------------
 */

#ifndef THERMOMECHANICAL_MODEL_H
#define THERMOMECHANICAL_MODEL_H

// libMesh includes
#include "heat_source.h"
#include "libmesh_includes.h"

// Local includes
#include "material.h"
#include "mpi_communicator.h"
#include "io.h"
#include "crack.h"

namespace loading {
  class HeatSource;
  class Loading;
}

namespace geom {
  class Fracture;
}


namespace model {

/*!
 * @brief Class implementing coupled thermomechanical model using libMesh and nodal peridynamics
 */
class ThermomechanicalModel {
public:
  /*!
   * @brief Constructor
   * @param mesh libMesh mesh object
   * @param deck Material deck containing parameters
   * @param dt Time step size
   * @param heat_source_p Pointer to heat source object
   */
  ThermomechanicalModel(libMesh::ReplicatedMesh& mesh, 
    libMesh::EquationSystems& equation_systems,
    libMesh::TransientLinearImplicitSystem& temperature_system,
    libMesh::ExplicitSystem& theta_dot_system,
    libMesh::ExplicitSystem& mechanical_system,
    inp::MaterialDeck& deck, 
    double dt);

  /*!
   * @brief Initialize the systems and compute static quantities
   */
  void initialize();

  void secondaryInitialize();

  /*!
   * @brief Set the fixity to free (0) or fixed (1)
   * @param i Id of node
   * @param dof Dof which is affected
   * @param flag Set fixity to fixed if true or free
   */
  void setFixity(const size_t &i, const unsigned int &dof, const bool &flag);

  void setupNeighborList();

  void computeNeighborVolume();

  void computeMx();

  void setupGhostNodesAndCommunicator();

  void syncGhostData();

  void assembleHeatMatrix();

  void assembleHeatRHS(const double& time);

  void updateThetaAndDamage();

  void advance();

  void solveHeatEquation();

  void copyTemperature();

  void copyThetaDot();

  void computeForces();
  
  void computePeriForces();
  
  void updateKinematics();

  void updateMechanicalSystem();

  void write(unsigned int file_number, bool print_debug = false);

  void setObservationPoints(const std::vector<libMesh::Point> &obs_points);

  void updateObsData();

  void updateQoi();

  double vectorNorm(const std::vector<libMesh::Point> &vec);

  double vectorNorm(const std::vector<double> &vec);

  void updateCoupledData();

  void debugVector(const std::vector<double> &vec, const std::string &name);

  void debugVector(const std::vector<libMesh::Point> &vec, const std::string &name);

  /*!
   * @brief Update dilation at nodes
   */
  

  /*!
   * @brief Compute peridynamic forces at nodes
   */
  

  /*!
   * @brief Update displacement and velocity using central difference
   * @param dt Time step size
   */
  

  // libMesh objects
  libMesh::ReplicatedMesh& d_mesh;                      //!< Mesh
  libMesh::EquationSystems& d_equation_systems; //!< Equation systems

  double d_dt;
  double d_time;
  std::string d_displacement_update_method;
  bool d_use_nodal_fem;

  // Systems for temperature field
  libMesh::TransientLinearImplicitSystem& d_temperature_system;  //!< Temperature system
  libMesh::ExplicitSystem& d_theta_dot_system;         //!< Volumetric strain rate system
  libMesh::ExplicitSystem& d_mechanical_system;         //!< Mechanical system

  // Local vectors for peridynamics
  std::vector<libMesh::Point> d_displacement;  //!< Displacement at nodes
  std::vector<libMesh::Point> d_displacement_old;  //!< Displacement at nodes
  std::vector<libMesh::Point> d_velocity;      //!< Velocity at nodes
  std::vector<double> d_temperature;           //!< Temperature at nodes
  std::vector<double> d_mx;                 //!< Weighted volume at nodes
  std::vector<double> d_theta;             //!< Volumetric strain at nodes
  std::vector<double> d_theta_old;         //!< Volumetric strain at nodes
  std::vector<double> d_theta_dot;         //!< Volumetric strain rate at nodes
  std::vector<double> d_damage;         //!< Damage at nodes
  std::vector<libMesh::Point> d_force;     //!< Force at nodes
  std::vector<std::vector<libMesh::dof_id_type>> d_neighbor_list; //!< Neighbor list at nodes
  std::vector<std::vector<double>> d_neighbor_volume;       //!< Neighbor volume for integration
  std::vector<double> d_nodal_volume; //!< Nodal volume for integration
  std::vector<uint8_t> d_displacement_fixed;  //!< Fixed nodes (three bits for x, y, z)
  std::vector<uint8_t> d_force_fixed;  //!< Specified force nodes (three bits for x, y, z)

  // fracture
  std::vector<geom::EdgeCrack> d_edge_cracks; //!< Edge cracks
  std::unique_ptr<geom::Fracture> d_fracture_p; //!< Fracture model

  // Material properties
  std::shared_ptr<material::Material> d_material_p;  //!< Material model

  // MPI communicator for parallel data exchange
  std::unique_ptr<MPICommunicator> d_cm_p;  //!< MPI communicator

  // Heat source
  std::shared_ptr<loading::HeatSourceCollection> d_heat_sources_p;

  // Loading
  std::unique_ptr<loading::Loading> d_loading_p;

  // observation points
  std::vector<libMesh::Point> d_obs_points;
  std::vector<libMesh::dof_id_type> d_obs_nodes;
  std::vector<bool> d_obs_points_owned;
  std::vector<double> d_obs_T;
  std::vector<std::vector<double>> d_obs_u;
  std::vector<double> d_obs_damage;

  std::map<std::string, double> d_qoi;

  std::string printStr(int nt = 0, int lvl = 0) const; 
  void print(int nt = 0, int lvl = 0) const;
};

} // namespace model

#endif // THERMOMECHANICAL_MODEL_H 