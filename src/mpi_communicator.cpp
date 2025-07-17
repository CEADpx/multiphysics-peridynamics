#include "mpi_communicator.h"
#include "util.h"

namespace model {

MPICommunicator::MPICommunicator(libMesh::ReplicatedMesh& mesh)
  : d_mesh(mesh),
    d_comm(mesh.comm()),
    d_rank(d_comm.rank()),
    d_n_procs(d_comm.size()) {}

void MPICommunicator::initCommunication() {
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
}

void MPICommunicator::syncDisplacement(std::vector<libMesh::Point>& displacement) {
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
}

void MPICommunicator::syncScalarData(std::vector<double>& data) {
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
}

void MPICommunicator::packPointData(const std::vector<libMesh::Point>& data, 
                                  unsigned int proc_id) {
  auto& proc_data = d_proc_comm[proc_id];
  for (size_t i = 0; i < proc_data.send_list.size(); ++i) {
    proc_data.point_send_buffer[i] = data[proc_data.send_list[i]];
  }
}

void MPICommunicator::packScalarData(const std::vector<double>& data, 
                                   unsigned int proc_id) {
  auto& proc_data = d_proc_comm[proc_id];
  for (size_t i = 0; i < proc_data.send_list.size(); ++i) {
    proc_data.scalar_send_buffer[i] = data[proc_data.send_list[i]];
  }
}

void MPICommunicator::unpackPointData(std::vector<libMesh::Point>& data, 
                                    unsigned int proc_id) {
  auto& proc_data = d_proc_comm[proc_id];
  for (size_t i = 0; i < proc_data.recv_list.size(); ++i) {
    data[proc_data.recv_list[i]] = proc_data.point_recv_buffer[i];
  }
}

void MPICommunicator::unpackScalarData(std::vector<double>& data, 
                                     unsigned int proc_id) {
  auto& proc_data = d_proc_comm[proc_id];
  for (size_t i = 0; i < proc_data.recv_list.size(); ++i) {
    data[proc_data.recv_list[i]] = proc_data.scalar_recv_buffer[i];
  }
}

} // namespace model 