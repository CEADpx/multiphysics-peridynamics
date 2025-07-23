#ifndef LOADING_H
#define LOADING_H

#include "bc.h"
#include "thermomechanical_model.h"
#include <libmesh/point.h>
#include <vector>
#include <cmath>
#include <map>
#include <utility>
#include <memory>

namespace loading {

/**
 * @brief Class to handle all loading conditions (displacement, force, initial)
 */
class Loading {
public:

  /**
   * @brief Default constructor
   */
  Loading(model::ThermomechanicalModel &model)
      : d_model(model) {};

  /**
   * @brief Constructor
   * @param model_p Pointer to thermomechanical model
   * @param disp_bcs Displacement boundary conditions
   * @param force_bcs Force boundary conditions
   * @param init_conds Initial conditions
   */
  Loading(model::ThermomechanicalModel &model,
          const std::vector<inp::BCBase> &disp_bcs,
          const std::vector<inp::BCBase> &force_bcs,
          const std::vector<inp::BCBase> &init_conds)
      : d_model(model), d_disp_bcs(disp_bcs), d_force_bcs(force_bcs), 
        d_init_conds(init_conds) {
    initialize();
  }

  void initialize() {
    buildNodeSets();
    setFixity(true);
    setFixity(false);
  }

  /**
   * @brief Build sets of nodes where BCs are active
   */
  void buildNodeSets() {
    bool debug = false;

    // Clear existing node sets
    d_disp_nodes.clear();
    d_force_nodes.clear();
    d_init_nodes.clear();
    d_disp_nodes.resize(d_disp_bcs.size());
    d_force_nodes.resize(d_force_bcs.size());
    d_init_nodes.resize(d_init_conds.size());

    // Get mesh from model
    const auto& mesh = d_model.d_mesh;

    size_t n_disp_nodes = 0;
    size_t n_force_nodes = 0;
    size_t n_init_nodes = 0;
    for (size_t loc_i=0; loc_i < d_model.d_cm_p->d_owned_size; loc_i++) {

      const auto& i = d_model.d_cm_p->d_owned_and_ghost_ids[loc_i];
      const libMesh::Point& xi = *mesh.node_ptr(i);

      for (size_t j=0; j < d_disp_bcs.size(); j++) {
        const auto& bc = d_disp_bcs[j];
        if (bc.isInRegion(xi)) {
          d_disp_nodes[j].push_back(i);
          n_disp_nodes++;
          if (debug)
            printf("Disp BC: BC id = %zu, Node = %u, Coordinates = (%f, %f, %f), Owner = %d, Rank = %d\n", j, i, xi(0), xi(1), xi(2), mesh.node_ptr(i)->processor_id(), d_model.d_cm_p->d_rank);
        }
      }

      for (size_t j=0; j < d_force_bcs.size(); j++) {
        const auto& bc = d_force_bcs[j];
        if (bc.isInRegion(xi)) {
          d_force_nodes[j].push_back(i);
          n_force_nodes++;
          if (debug)
            printf("Force BC: BC id = %zu, Node = %u, Coordinates = (%f, %f, %f), Owner = %d, Rank = %d\n", j, i, xi(0), xi(1), xi(2), mesh.node_ptr(i)->processor_id(), d_model.d_cm_p->d_rank);
        }
      }

      for (size_t j=0; j < d_init_conds.size(); j++) {
        const auto& bc = d_init_conds[j];
        if (bc.isInRegion(xi)) {
          d_init_nodes[j].push_back(i);
          n_init_nodes++;
          if (debug)
            printf("Init BC: BC id = %zu, Node = %u, Coordinates = (%f, %f, %f), Owner = %d, Rank = %d\n", j, i, xi(0), xi(1), xi(2), mesh.node_ptr(i)->processor_id(), d_model.d_cm_p->d_rank);
        }
      }
    }
    if (debug) {
      printf("n_disp_nodes = %zu\n", n_disp_nodes);
      printf("n_force_nodes = %zu\n", n_force_nodes);
      printf("n_init_nodes = %zu\n", n_init_nodes);
    }
  }

  /**
   * @brief Apply displacement boundary conditions
   * @param time Current time
   */
  void applyDisplacement(const double time) const {
    for (size_t i = 0; i < d_disp_bcs.size(); ++i) {
      const auto& bc = d_disp_bcs[i];
      const auto& nodes = d_disp_nodes[i];
      
      for (const auto& node_id : nodes) {
        if (bc.d_is_zero) continue;

        libMesh::Point ux, vx;
        bc.getBCDispVel(time, *d_model.d_mesh.node_ptr(node_id), ux, vx);
        for (auto d : bc.d_direction) {
          d_model.d_displacement_old[node_id](d) = d_model.d_displacement[node_id](d);
          d_model.d_displacement[node_id](d) = ux(d);
          d_model.d_velocity[node_id](d) = vx(d);
        }
      }
    }
  }

  /**
   * @brief Apply force boundary conditions
   * @param time Current time
   */
  void applyForce(const double time) const {
    for (size_t i = 0; i < d_force_bcs.size(); ++i) {
      const auto& bc = d_force_bcs[i];
      const auto& nodes = d_force_nodes[i];
      
      for (const auto& node_id : nodes) {
        libMesh::Point fx;
        bc.getBCForce(time, *d_model.d_mesh.node_ptr(node_id), fx);
        for (auto d : bc.d_direction) {
          d_model.d_force[node_id](d) += fx(d);
        }
      }
    }
  }

  /**
   * @brief Apply initial conditions
   */
  void applyInitialCondition() const {
    for (size_t i = 0; i < d_init_conds.size(); ++i) {
      const auto& ic = d_init_conds[i];
      const auto& nodes = d_init_nodes[i];

      for (const auto& node_id : nodes) {
        d_model.d_velocity[node_id] = ic.d_ic_vector;
      }
    }
  }

  /**
   * @brief Set fixity for nodes based on boundary conditions
   * @param fixity_for_disp If true, set fixity for displacement BCs, if false for force BCs
   */
  void setFixity(bool fixity_for_disp = true) const {

    bool debug = false;

    if (fixity_for_disp) {
      // Set fixity for displacement BCs
      for (size_t i = 0; i < d_disp_bcs.size(); ++i) {
        const auto& bc = d_disp_bcs[i];
        const auto& nodes = d_disp_nodes[i];

        for (const auto& node_id : nodes) {
          for (const auto& dir : bc.d_direction) {
            bool is_free;
            if (debug) is_free = util::isDofFree(node_id, dir, d_model.d_displacement_fixed); 
            
            util::setFixity(node_id, dir, true, d_model.d_displacement_fixed);
            
            if (debug) {
              bool is_free_after = util::isDofFree(node_id, dir, d_model.d_displacement_fixed);
              printf("Disp BC: BC id = %zu, Node = %zu, dir = %zu, is_free = %s, is_free_after = %s\n", 
                  i, node_id, dir, is_free ? "true" : "false", is_free_after ? "true" : "false" );
            }
          }
        }
      }
    } else {
      // Set fixity for force BCs
      for (size_t i = 0; i < d_force_bcs.size(); ++i) {
        const auto& bc = d_force_bcs[i];
        const auto& nodes = d_force_nodes[i];

        for (const auto& node_id : nodes) {
          for (const auto& dir : bc.d_direction) {
            util::setFixity(node_id, dir, true, d_model.d_force_fixed);
          }
        }
      }
    }
  }


  model::ThermomechanicalModel &d_model;  ///< thermomechanical model

  std::vector<inp::BCBase> d_disp_bcs;   ///< Displacement boundary conditions
  std::vector<inp::BCBase> d_force_bcs;   ///< Force boundary conditions
  std::vector<inp::BCBase> d_init_conds;  ///< Initial conditions

  std::vector<std::vector<size_t>> d_disp_nodes;   ///< Nodes for each displacement BC
  std::vector<std::vector<size_t>> d_force_nodes;  ///< Nodes for each force BC
  std::vector<std::vector<size_t>> d_init_nodes;   ///< Nodes for each initial condition

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
    oss << tabS << "------- loading::Loading --------" << std::endl
        << std::endl;
    oss << tabS << "Disp BCs: " << std::endl;
    for (size_t i = 0; i < d_disp_bcs.size(); ++i) {
      oss << d_disp_bcs[i].printStr(nt + 1, lvl) << std::endl;
      oss << tabS << tabS << "Number of nodes: " << d_disp_nodes[i].size() << std::endl;
    }
    oss << tabS << "Force BCs: " << std::endl;
    for (size_t i = 0; i < d_force_bcs.size(); ++i) {
      oss << d_force_bcs[i].printStr(nt + 1, lvl) << std::endl;
      oss << tabS << tabS << "Number of nodes: " << d_force_nodes[i].size() << std::endl;
    }
    oss << tabS << "Init BCs: " << std::endl;
    for (size_t i = 0; i < d_init_conds.size(); ++i) {
      oss << d_init_conds[i].printStr(nt + 1, lvl) << std::endl;
      oss << tabS << tabS << "Number of nodes: " << d_init_nodes[i].size() << std::endl;
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

#endif // LOADING_H 