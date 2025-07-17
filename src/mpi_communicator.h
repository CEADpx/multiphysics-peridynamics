#ifndef MPI_COMMUNICATOR_H
#define MPI_COMMUNICATOR_H

#include "libmesh_includes.h"
#include <vector>
#include <map>
#include "io.h"

namespace model {

/*!
 * @brief Class for handling MPI communication of nodal data
 */
class MPICommunicator {
public:
  /*!
   * @brief Constructor
   * @param mesh Reference to libMesh mesh object
   * @param ghost_ids Vector of ghost node IDs
   */
  MPICommunicator(libMesh::ReplicatedMesh& mesh);

  /*!
   * @brief Synchronize displacement data
   * @param displacement Vector of nodal displacements
   */
  void syncDisplacement(std::vector<libMesh::Point>& displacement);

  /*!
   * @brief Synchronize scalar data (temperature, theta)
   * @param data Vector of nodal scalar data
   */
  void syncScalarData(std::vector<double>& data);

  
  libMesh::ReplicatedMesh& d_mesh;                    //!< Reference to mesh
  const libMesh::Parallel::Communicator& d_comm; //!< MPI communicator
  const unsigned int d_rank;                //!< Current processor rank
  const unsigned int d_n_procs;             //!< Number of processors
  std::vector<std::vector<libMesh::dof_id_type>> d_ghost_ids;     //!< Global IDs of ghost nodes for each owned node
  std::vector<libMesh::dof_id_type> d_owned_and_ghost_ids; //!< Global IDs of owned and ghost nodes
  libMesh::dof_id_type d_owned_size; //!< Number of owned nodes

  // Pre-computed communication data for each processor
  struct ProcessorComm {
    std::vector<libMesh::dof_id_type> send_list;    //!< List of nodes to send
    std::vector<libMesh::dof_id_type> recv_list;    //!< List of nodes to receive
    std::vector<libMesh::Point> point_send_buffer; //!< Reusable buffer for Point data
    std::vector<libMesh::Point> point_recv_buffer; //!< Reusable buffer for Point data
    std::vector<double> scalar_send_buffer; //!< Reusable buffer for scalar data
    std::vector<double> scalar_recv_buffer; //!< Reusable buffer for scalar data
  };
  std::vector<ProcessorComm> d_proc_comm;   //!< Communication data for each processor

  /*!
   * @brief Initialize communication data structures
   * @param ghost_ids Vector of ghost node IDs
   */
  void initCommunication();

  /*!
   * @brief Pack Point data into pre-allocated buffer
   * @param data Source data
   * @param proc_id Processor ID
   */
  void packPointData(const std::vector<libMesh::Point>& data, unsigned int proc_id);

  /*!
   * @brief Pack scalar data into pre-allocated buffer
   * @param data Source data
   * @param proc_id Processor ID
   */
  void packScalarData(const std::vector<double>& data, unsigned int proc_id);

  /*!
   * @brief Unpack Point data from pre-allocated buffer
   * @param data Target data to update
   * @param proc_id Processor ID
   */
  void unpackPointData(std::vector<libMesh::Point>& data, unsigned int proc_id);

  /*!
   * @brief Unpack scalar data from pre-allocated buffer
   * @param data Target data to update
   * @param proc_id Processor ID
   */
  void unpackScalarData(std::vector<double>& data, unsigned int proc_id);

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

#endif // MPI_COMMUNICATOR_H 