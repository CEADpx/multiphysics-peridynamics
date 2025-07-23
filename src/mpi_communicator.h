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

#include "libmesh_includes.h"
#include <vector>
#include "io.h"
#include "util.h"

namespace model {

struct ProcessorComm {
  std::vector<libMesh::dof_id_type> send_list;    //!< List of nodes to send
  std::vector<libMesh::dof_id_type> recv_list;    //!< List of nodes to receive
  std::vector<libMesh::Point> point_send_buffer; //!< Reusable buffer for Point data
  std::vector<libMesh::Point> point_recv_buffer; //!< Reusable buffer for Point data
  std::vector<double> scalar_send_buffer; //!< Reusable buffer for scalar data
  std::vector<double> scalar_recv_buffer; //!< Reusable buffer for scalar data
};

/*!
 * @brief Class for handling MPI communication of nodal data
 */
class MPICommunicator {
public:
  libMesh::ReplicatedMesh& d_mesh;                    //!< Reference to mesh
  const libMesh::Parallel::Communicator& d_comm; //!< MPI communicator
  const unsigned int d_rank;                //!< Current processor rank
  const unsigned int d_n_procs;             //!< Number of processors
  std::vector<std::vector<libMesh::dof_id_type>> d_ghost_ids;     //!< Global IDs of ghost nodes for each owned node
  std::vector<libMesh::dof_id_type> d_owned_and_ghost_ids; //!< Global IDs of owned and ghost nodes
  libMesh::dof_id_type d_owned_size; //!< Number of owned nodes

  // Pre-computed communication data for each processor
  std::vector<ProcessorComm> d_proc_comm;   //!< Communication data for each processor

  /*!
   * @brief Constructor
   * @param mesh Reference to libMesh mesh object
   * @param ghost_ids Vector of ghost node IDs
   */
  MPICommunicator(libMesh::ReplicatedMesh& mesh)
  : d_mesh(mesh),
    d_comm(mesh.comm()),
    d_rank(d_comm.rank()),
    d_n_procs(d_comm.size()) {};

  /*!
   * @brief Synchronize displacement data
   * @param displacement Vector of nodal displacements
   */
  void syncDisplacement(std::vector<libMesh::Point>& displacement) {
    // For each processor, exchange data
    for (unsigned int p = 0; p < d_n_procs; ++p) {
      if (p == d_rank) continue;
  
      auto& proc_data = d_proc_comm[p];
      if (proc_data.send_list.empty() && proc_data.recv_list.empty()) continue;
  
      // Pack data into pre-allocated buffer
      packPointData(displacement, p);
  
      // Exchange data
      d_comm.send_receive(p, proc_data.point_send_buffer,
                         p, proc_data.point_recv_buffer);
  
      // Unpack received data
      unpackPointData(displacement, p);
    }
  };

  /*!
   * @brief Synchronize scalar data (temperature, theta)
   * @param data Vector of nodal scalar data
   */
  void syncScalarData(std::vector<double>& data) {
    // For each processor, exchange data
    for (unsigned int p = 0; p < d_n_procs; ++p) {
      if (p == d_rank) continue;
  
      auto& proc_data = d_proc_comm[p];
      if (proc_data.send_list.empty() && proc_data.recv_list.empty()) continue;
  
      // Pack data into pre-allocated buffer
      packScalarData(data, p);
  
      // Exchange data
      d_comm.send_receive(p, proc_data.scalar_send_buffer,
                         p, proc_data.scalar_recv_buffer);
  
      // Unpack received data
      unpackScalarData(data, p);
    }
  };

  /*!
   * @brief Initialize communication data structures
   * @param ghost_ids Vector of ghost node IDs
   */
  void initCommunication() {
    // Resize processor communication data
    d_proc_comm.resize(d_n_procs);
  
    // For each processor, determine which nodes to send/receive
    for (unsigned int p = 0; p < d_n_procs; ++p) {
      if (p == d_rank) continue;
  
      auto& proc_data = d_proc_comm[p];
  
      // We use ghost list to find nodes to send (node to ghost id belongs to) and receive (ghost id of node)
      for (libMesh::dof_id_type loc_i = 0; loc_i < d_owned_size; loc_i++) {
        const auto i = d_owned_and_ghost_ids[loc_i];  // Global node ID
        const libMesh::Point& xi = *d_mesh.node_ptr(i);
        // check that d_rank owns this node if not throw error
        if (d_mesh.node_ptr(i)->processor_id() != d_rank) {
          throw std::runtime_error("Node " + std::to_string(i) + " must be owned by processor " + std::to_string(d_rank));
        }
  
        if (d_ghost_ids[i].size() == 0) continue;
  
        for (const auto& ghost_id : d_ghost_ids[i]) {
          // find the owner of the ghost id
          const bool is_p_owner = d_mesh.node_ptr(ghost_id)->processor_id() == p;
          if (is_p_owner) {
            util::addUnique(proc_data.send_list, i);
            util::addUnique(proc_data.recv_list, ghost_id);
          }
        }
      }
  
      // Pre-allocate buffers
      proc_data.point_send_buffer.resize(proc_data.send_list.size());
      proc_data.point_recv_buffer.resize(proc_data.recv_list.size());
      proc_data.scalar_send_buffer.resize(proc_data.send_list.size());
      proc_data.scalar_recv_buffer.resize(proc_data.recv_list.size());
    }
  };

  /*!
   * @brief Pack Point data into pre-allocated buffer
   * @param data Source data
   * @param proc_id Processor ID
   */
  void packPointData(const std::vector<libMesh::Point>& data, unsigned int proc_id) {
    auto& proc_data = d_proc_comm[proc_id];
    for (size_t i = 0; i < proc_data.send_list.size(); ++i) {
      proc_data.point_send_buffer[i] = data[proc_data.send_list[i]];
    }
  };

  /*!
   * @brief Pack scalar data into pre-allocated buffer
   * @param data Source data
   * @param proc_id Processor ID
   */
  void packScalarData(const std::vector<double>& data, unsigned int proc_id) {
    auto& proc_data = d_proc_comm[proc_id];
    for (size_t i = 0; i < proc_data.send_list.size(); ++i) {
      proc_data.scalar_send_buffer[i] = data[proc_data.send_list[i]];
    }
  };

  /*!
   * @brief Unpack Point data from pre-allocated buffer
   * @param data Target data to update
   * @param proc_id Processor ID
   */
  void unpackPointData(std::vector<libMesh::Point>& data, unsigned int proc_id) {
    auto& proc_data = d_proc_comm[proc_id];
    for (size_t i = 0; i < proc_data.recv_list.size(); ++i) {
      data[proc_data.recv_list[i]] = proc_data.point_recv_buffer[i];
    }
  };

  /*!
   * @brief Unpack scalar data from pre-allocated buffer
   * @param data Target data to update
   * @param proc_id Processor ID
   */
  void unpackScalarData(std::vector<double>& data, unsigned int proc_id) {
    auto& proc_data = d_proc_comm[proc_id];
    for (size_t i = 0; i < proc_data.recv_list.size(); ++i) {
      data[proc_data.recv_list[i]] = proc_data.scalar_recv_buffer[i];
    }
  };

  std::string printStr(int nt = 0, int lvl = 0) const {
    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- model::MPICommunicator --------" << std::endl
        << std::endl;
    oss << tabS << "rank = " << d_rank << std::endl;
    oss << tabS << "n_procs = " << d_n_procs << std::endl;
    oss << tabS << "owned_size = " << d_owned_size << std::endl;
    oss << tabS << "ghost_size = " << d_owned_and_ghost_ids.size() - d_owned_size << std::endl;
    return oss.str(); 
  }

  void print(int nt = 0, int lvl = 0) const {
    std::cout << printStr(nt, lvl);
  }
};

} // namespace model