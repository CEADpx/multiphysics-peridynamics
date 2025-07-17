
#include "geometry.h"
#include "thermomechanical_model.h"
#include "libmesh_includes.h"

#include "bc.h"
#include "heat_source.h"
#include "io.h"
#include "loading.h"
#include "material.h"
#include "mpi.h"

void checkMesh(libMesh::ReplicatedMesh &mesh);

void applyTemperatureBoundaryConditions(libMesh::EquationSystems &equation_systems, const inp::MaterialDeck &material_deck, const double &T0, const double & T1);

void setDisplacementAndForceConditions(std::shared_ptr<model::ThermomechanicalModel> &model_p, const geom::Cuboid &domain, const double &tFinal, const double &load_rate);

void checkMPICommunicatorAndExchange(const std::shared_ptr<model::ThermomechanicalModel>& model_p);

void checkMPIReduce(const std::shared_ptr<model::ThermomechanicalModel>& model_p);

void debugLoading(const std::shared_ptr<model::ThermomechanicalModel>& model_p);

void debugPeriForce(const std::shared_ptr<model::ThermomechanicalModel>& model_p);

void setObservationPoints(std::shared_ptr<model::ThermomechanicalModel> &model_p, const geom::Cuboid &domain);

int main(int argc, char** argv) {
  // Initialize libMesh
  libMesh::LibMeshInit init(argc, argv);
  util::io::setCommunicator(init.comm().rank());

  // Problem parameters
  int dim = 2;
  double Lx = 4.0, Ly = 1.0, Lz = 0.0;
  if (dim == 3) {
    Lz = 1.0;
  }
  libMesh::Point center(0.5*Lx, 0.5*Ly, 0.5*Lz);
  geom::Cuboid domain(Lx, Ly, Lz, center);
  const unsigned int nx = 80, ny = 20, nz = 20; // Number of elements in each direction
  const double T0 = 273.0; // Initial temperature in K
  const double T1 = 273.0; // Dirichlet BC at x = L (if no robinBC)
  const double horizon = 4*(domain.d_Lx/nx); // Example horizon, should match material deck
  const double load_rate = 0.0; // Load rate in N/s

  // Create mesh
  libMesh::Mesh mesh(init.comm());
  if (dim == 2) {
    libMesh::MeshTools::Generation::build_square(mesh, nx, ny, 0., domain.d_Lx, 0., domain.d_Ly, libMesh::QUAD4);
  } else {
    libMesh::MeshTools::Generation::build_cube(mesh, nx, ny, nz, 0., domain.d_Lx, 0., domain.d_Ly, 0., domain.d_Lz, libMesh::HEX8);
  }

  //checkMesh(mesh);

  // Setup input deck
  inp::MaterialDeck material_deck;
  material_deck.setDefaults();
  material_deck.d_dim = dim;
  material_deck.d_horizon = horizon;
  double dt = 0.001;
  int nsteps = 1000;
  double tFinal = dt*nsteps;
  int write_interval = 10;

  // Create a Gaussian heat source at the center of the domain
  std::string sfn_type = "gaussian";
  std::vector<double> sfn_params = {domain.d_Ly/10, 1.0, domain.d_x(0), domain.d_x(1), domain.d_x(2)};
  std::string tfn_type = "linear_step_const_value";
  std::vector<double> tfn_params = {-50.0/tFinal, 0.0, 0.1*tFinal, tFinal};
  auto heat_source_p = std::make_unique<loading::HeatSource>(sfn_type, sfn_params, tfn_type, tfn_params);

  // Create equation systems
  libMesh::EquationSystems equation_systems(mesh);
  auto& temperature_system = equation_systems.add_system<libMesh::TransientLinearImplicitSystem>
    ("Temperature");
  temperature_system.add_variable("temp", libMesh::FIRST);

  auto& theta_dot_system = equation_systems.add_system<libMesh::ExplicitSystem>
    ("ThetaDot");
  theta_dot_system.add_variable("theta_dot", libMesh::FIRST);

  auto& mechanical_system = equation_systems.add_system<libMesh::ExplicitSystem>
    ("Mechanical");
  mechanical_system.add_variable("ux", libMesh::FIRST);
  mechanical_system.add_variable("uy", libMesh::FIRST);
  if (mesh.mesh_dimension() == 3) {
    mechanical_system.add_variable("uz", libMesh::FIRST);
  }
  mechanical_system.add_variable("vx", libMesh::FIRST);
  mechanical_system.add_variable("vy", libMesh::FIRST);
  if (mesh.mesh_dimension() == 3) {
    mechanical_system.add_variable("vz", libMesh::FIRST);
  }
  mechanical_system.add_variable("fx", libMesh::FIRST);
  mechanical_system.add_variable("fy", libMesh::FIRST);
  if (mesh.mesh_dimension() == 3) {
    mechanical_system.add_variable("fz", libMesh::FIRST);
  }
  mechanical_system.add_variable("theta", libMesh::FIRST);

  // we need to avoid libmesh resetting the matrix and rhs to zero
  temperature_system.zero_out_matrix_and_rhs = false;

  // Apply temperature conditions
  applyTemperatureBoundaryConditions(equation_systems, material_deck, T0, T1);

  // initialize equation systems
  equation_systems.init();

  // Create model
  auto model_p = std::make_shared<model::ThermomechanicalModel>(mesh, equation_systems, 
          temperature_system, theta_dot_system, mechanical_system, 
          material_deck, dt, std::move(heat_source_p));

  // Initialize model
  model_p->initialize();

  setDisplacementAndForceConditions(model_p, domain, tFinal, load_rate);

  model_p->secondaryInitialize();

  // get nodes at obs points
  setObservationPoints(model_p, domain);

  if (false) {
    std::cout << "\n\n\ndebugging loading" << std::endl;
    debugLoading(model_p);
    std::cout << "debugging loading done\n\n" << std::endl;
  }

  if (false){
    std::cout << "\n\n\ndebugging MPICommunicator and exchange" << std::endl;
    checkMPICommunicatorAndExchange(model_p);
    checkMPIReduce(model_p);
    std::cout << "debugging MPICommunicator and exchange done\n\n" << std::endl;
  }

  if (false) {
    std::cout << "Material\n\n" << std::endl;
    model_p->d_material_p->print();
    std::cout << "\n\n\n" << std::endl;
  }

  std::cout << "\n\n++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "ThermomechanicalModel information" << std::endl;
  model_p->print();
  std::cout << "++++++++++++++++++++++++++++++++++++\n\n" << std::endl;

  // Time stepping loop
  if (init.comm().rank() == 0) {
    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "Time stepping loop" << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
  }

  if (init.comm().rank() == 0) {
    std::cout << "\n\n\n++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "Step " << 0 << ", time = " << model_p->d_time << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
  }
  model_p->write(0, true);

  for (unsigned int step = 1; step <= nsteps; ++step) {

    if (step% write_interval == 0) {
      if (init.comm().rank() == 0) {
        std::cout << "\n\n\n++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "Step " << step << ", time = " << model_p->d_time << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
      }
    }

    // advance the model
    model_p->advance();

    if (step % write_interval == 0) {
      model_p->write(step/write_interval + 1, true);

      // print heat source at the center of the domain
      std::cout << "Heat source at the center of the domain: " << model_p->d_heat_source_p->get(domain.d_x, model_p->d_time) << std::endl;
    }

    if (step % write_interval == 0 && false) {
      debugPeriForce(model_p);
    }
  }

  return 0;
}

void checkMesh(libMesh::ReplicatedMesh &mesh) {
  // Comprehensive node numbering check
  printf("=== Node Numbering Consistency Check ===\n");
  printf("n_nodes = %d\n", mesh.n_nodes());
  printf("max_node_id = %d\n", mesh.max_node_id());
  
  // Collect all node IDs
  std::set<libMesh::dof_id_type> node_ids;
  std::vector<libMesh::dof_id_type> node_id_vector;
  for (const auto& node : mesh.node_ptr_range()) {
    node_ids.insert(node->id());
    node_id_vector.push_back(node->id());
  }
  
  printf("Unique node IDs found: %zu\n", node_ids.size());
  printf("Total nodes in mesh: %d\n", mesh.n_nodes());
  
  // Check if node IDs are consecutive starting from 0
  bool consecutive = true;
  libMesh::dof_id_type expected_id = 0;
  for (const auto& id : node_ids) {
    if (id != expected_id) {
      printf("ERROR: Expected node ID %d, but found %d\n", expected_id, id);
      consecutive = false;
      break;
    }
    expected_id++;
  }
  
  if (consecutive) {
    printf("✓ Node IDs are consecutive starting from 0\n");
  } else {
    printf("✗ Node IDs are NOT consecutive starting from 0\n");
  }
  
  // Check if max_node_id matches the highest node ID
  if (node_ids.empty()) {
    printf("ERROR: No nodes found in mesh!\n");
    exit(EXIT_FAILURE);
  }
  
  // Check if we can access all nodes by their IDs
  bool all_nodes_accessible = true;
  for (const auto& id : node_ids) {
    if (id > mesh.max_node_id()) {
      printf("ERROR: Node ID %d > max_node_id %d\n", id, mesh.max_node_id());
      all_nodes_accessible = false;
    }
    const libMesh::Node* node = mesh.node_ptr(id);
    if (!node) {
      printf("ERROR: Cannot access node with ID %d\n", id);
      all_nodes_accessible = false;
    }
  }
  
  if (all_nodes_accessible) {
    printf("✓ All nodes are accessible by their IDs\n");
  } else {
    printf("✗ Some nodes are not accessible by their IDs\n");
  }
  
  printf("=== End Node Numbering Check ===\n\n");
  
  // If node numbering is not consistent, exit
  if (!consecutive || !all_nodes_accessible) {
    printf("ERROR: Node numbering is inconsistent. Exiting.\n");
    exit(EXIT_FAILURE);
  }
}

libMesh::Number initialTemperature(const libMesh::Point &p, const libMesh::Parameters &es,
  const std::string &system_name, const std::string &var_name) {
    libmesh_assert_equal_to(system_name, "Temperature");
    if (var_name == "temp") {
      return es.get<double>("T0");
    }
    return 0.;
}

void initial_condition(libMesh::EquationSystems &es, const std::string &system_name) {
  if (system_name == "Temperature") {
    auto &sys = es.get_system<libMesh::TransientLinearImplicitSystem>(system_name);
    sys.project_solution(initialTemperature, nullptr, es.parameters);
  }
}

void applyTemperatureBoundaryConditions(libMesh::EquationSystems &equation_systems, const inp::MaterialDeck &material_deck, const double &T0, const double &T1) {
  // x=0 face: T=273K, x=L face: T=373K

  // Get the temperature system and dof map
  auto& temperature_system = equation_systems.get_system<libMesh::TransientLinearImplicitSystem>("Temperature");
  const auto& dof_map = temperature_system.get_dof_map();
  std::vector<libMesh::dof_id_type> dof_indices;

  // Set initial temperature everywhere
  equation_systems.parameters.set<double>("T0") = T0;
  temperature_system.attach_init_function(initial_condition);

  if (material_deck.d_robinBC) {
    return;
  }

  // Boundary conditions
  // create constant function for bc and apply
  std::set<libMesh::boundary_id_type> ids;
  {
    ids.insert(0);
    libMesh::ConstFunction<libMesh::Number> T0_fn(T0);
    libMesh::DirichletBoundary diri_bc(ids, {0}, &T0_fn);
    temperature_system.get_dof_map().add_dirichlet_boundary(diri_bc);
  }

  // create constant function for bc and apply
  if (false) {
    ids.clear();
    if (equation_systems.get_mesh().mesh_dimension() == 3) {
      ids.insert(2);
    } else if (equation_systems.get_mesh().mesh_dimension() == 2) {
      ids.insert(1);
    }
    libMesh::ConstFunction<libMesh::Number> T1_fn(T1);
    libMesh::DirichletBoundary diri_bc(ids, {0}, &T1_fn);
    temperature_system.get_dof_map().add_dirichlet_boundary(diri_bc);
  }
}

void setDisplacementAndForceConditions(std::shared_ptr<model::ThermomechanicalModel> &model_p, const geom::Cuboid &domain, const double &tFinal, const double &load_rate) {
  auto& loading_p = model_p->d_loading_p;

  // displacement BC: volume of thickness horizon from x = 0 to x = horizon
  const auto& horizon = model_p->d_material_p->d_deck.d_horizon;
  auto dim = model_p->d_mesh.mesh_dimension();
  {
    loading_p->d_disp_bcs.resize(1);
    auto& bc = loading_p->d_disp_bcs[0];
    bc.d_type = "Displacement_BC";
    auto left_face_center = domain.d_x - libMesh::Point(0.5*domain.d_Lx, 0.0, 0.0);
    bc.d_region_p = std::make_shared<geom::Cuboid>(horizon, domain.d_Ly, domain.d_Lz, left_face_center + libMesh::Point(0.5*horizon, 0.0, 0.0));
    bc.d_direction = {0, 1, 2};
    bc.d_is_zero = true;
  }
  
  // force BC: volume of thickness horizon from x = L - horizon to x = L
  if (false) {
    loading_p->d_force_bcs.resize(1);
    auto& bc = loading_p->d_force_bcs[0];
    bc.d_type = "Force_BC";
    auto right_face_center = domain.d_x + libMesh::Point(0.5*domain.d_Lx, 0.0, 0.0);
    bc.d_region_p = std::make_shared<geom::Cuboid>(horizon, domain.d_Ly, domain.d_Lz, right_face_center - libMesh::Point(0.5*horizon, 0.0, 0.0));
    bc.d_direction = {1};
    bc.d_time_fn_type = "linear_step";
    bc.d_time_fn_params = {load_rate, 0.25*tFinal, tFinal};
  }
}

// Function to check out MPICommunicator d_cm_p of ThermomechanicalModel
// and test ghost data exchange with fake data, and check if we received correct values
void checkMPICommunicatorAndExchange(const std::shared_ptr<model::ThermomechanicalModel>& model_p) {
  if (!model_p || !model_p->d_cm_p) {
    std::cout << "ThermomechanicalModel or its communicator is not initialized." << std::endl;
    return;
  }
  auto& cm = *model_p->d_cm_p;
  int rank = cm.d_comm.rank();
  int nproc = cm.d_comm.size();

  // Print communicator info
  if (rank == 0) {
    std::cout << "=== MPICommunicator info ===" << std::endl;
    std::cout << "Number of owned nodes: " << cm.d_owned_size << std::endl;
    std::cout << "Total owned_and_ghost_ids: " << cm.d_owned_and_ghost_ids.size() << std::endl;
    std::cout << "Sample owned_and_ghost_ids: ";
    for (size_t i = 0; i < std::min<size_t>(5, cm.d_owned_and_ghost_ids.size()); ++i) {
      std::cout << cm.d_owned_and_ghost_ids[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "===========================" << std::endl;
  }

  // Create fake data: assign to each node a value = rank*1000 + node_id
  std::vector<double> fake_data;
  fake_data.resize(model_p->d_mesh.n_nodes(), 0.0);
  for (size_t loc_i = 0; loc_i < cm.d_owned_and_ghost_ids.size(); ++loc_i) {
    auto node_id = cm.d_owned_and_ghost_ids[loc_i];
    const auto* node = model_p->d_mesh.node_ptr(node_id);
    int owner_rank = node->processor_id();
    fake_data[node_id] = owner_rank * 1000.0 + node_id;
  }

  // Save a copy of the data before exchange for checking
  std::vector<double> fake_data_before = fake_data;

  // Exchange fake data using communicator's syncScalarData
  cm.syncScalarData(fake_data);

  // Now, for each ghost node, check if the value matches the expected value from the owner rank
  // The expected value is: owner_rank * 1000 + ghost_id
  for (int proc = 0; proc < nproc; ++proc) {
    cm.d_comm.barrier();
    if (rank == proc) {
      std::cout << "[Rank " << rank << "] Ghost node values after exchange and check:" << std::endl;
      size_t count = 0;
      for (size_t i = 0; i < cm.d_ghost_ids.size(); ++i) {
        for (auto ghost_id : cm.d_ghost_ids[i]) {
          // Get owner rank of ghost_id using node_ptr
          const auto* ghost_node = model_p->d_mesh.node_ptr(ghost_id);
          int owner_rank = ghost_node->processor_id();
          double expected = owner_rank * 1000.0 + ghost_id;
          double received = (ghost_id < fake_data.size() ? fake_data[ghost_id] : -9999.0);
          bool correct = (std::abs(received - expected) < 1e-8);
          std::cout << "  ghost_id = " << ghost_id
                    << ", owner_rank = " << owner_rank
                    << ", received = " << received
                    << ", expected = " << expected
                    << "  --> " << (correct ? "OK" : "ERROR") << std::endl;
          if (++count >= 5) break;
        }
        if (count >= 5) break;
      }
      if (count == 0) {
        std::cout << "  (No ghost nodes on this rank)" << std::endl;
      }
    }
    cm.d_comm.barrier();
  }
}




void checkMPIReduce(const std::shared_ptr<model::ThermomechanicalModel>& model_p) {

  auto mpi_comm = model_p->d_cm_p->d_comm.get();
  auto my_rank = model_p->d_cm_p->d_comm.rank();
  auto nproc = model_p->d_cm_p->d_comm.size();
  printf("MPI rank %d entered checkMPIReduce()\n", my_rank);

  std::vector<double> local_obs_T(model_p->d_obs_nodes.size(), 0.0);
  std::vector<double> obs_T(model_p->d_obs_nodes.size(), 0.0);
  std::vector<double> obs_T_true(model_p->d_obs_nodes.size(), 0.0);

  const auto& obs_nodes = model_p->d_obs_nodes;
  const auto& obs_points_owned = model_p->d_obs_points_owned;

  {
    auto n = obs_nodes.size();
    // we verify that obs_points_owned is correct
    // every obs point should be owned by exactly one rank
    // and the sum of obs_points_owned should be equal to the number of obs points
    std::vector<int> obs_points_owned_count(n, 0);
    for (size_t i = 0; i < n; ++i) {
      obs_points_owned_count[i] = obs_points_owned[i] ? 1 : 0;
    }
    
    std::vector<int> obs_points_owned_count_sum(n, 0);
    MPI_Allreduce(obs_points_owned_count.data(), obs_points_owned_count_sum.data(), n, MPI_INT, MPI_SUM, mpi_comm);
    for (size_t i = 0; i < n; ++i) {
      if (my_rank == 0) {
        printf("obs_points_owned_count_sum[%zu] = %d, rank = %d\n", i, obs_points_owned_count_sum[i], my_rank);
      }
      if (obs_points_owned_count_sum[i] != 1) {
        throw std::runtime_error("obs_points_owned_count_sum[i] != 1");
      }
    }
  }

  int check_method = 2;

  for (size_t i = 0; i < obs_nodes.size(); ++i) {
    const auto& i_node = obs_nodes[i];
    if (check_method == 1) {
      local_obs_T[i] = std::pow(10, my_rank + 1);
      for (int j=0; j<nproc; j++) {
        obs_T_true[i] += std::pow(10, j + 1);
      }
    } else if (check_method == 2) {
      if (obs_points_owned[i]) {
        local_obs_T[i] = i_node*i_node;
      }
      obs_T_true[i] = i_node*i_node;
    }
  }

  // Sum all vectors and store result in processor 0
  auto comm = model_p->d_cm_p->d_comm.get();
  
  MPI_Reduce(
      local_obs_T.data(),          // send buffer
      obs_T.data(),         // receive buffer at root
      obs_nodes.size(),                  // number of elements
      MPI_DOUBLE,                // data type
      MPI_SUM,                   // operation
      0,                         // root rank
      comm             // communicator
  );

  if (my_rank == 0) {
    int n_errors = 0;
    for (size_t i = 0; i < obs_nodes.size(); ++i) {
      if (std::abs(obs_T[i] - obs_T_true[i]) > 1e-8) {
        n_errors++;
      }
      printf("Obs point %zu, T = %f, true T = %f  --> %s\n", i, obs_T[i], obs_T_true[i], (std::abs(obs_T[i] - obs_T_true[i]) < 1e-8 ? "OK" : "ERROR"));
    }
    if (n_errors > 0) {
      printf("n_errors = %d in rank = %d\n", n_errors, my_rank);
      throw std::runtime_error("checkMPIReduce() failed");
    } else {
      printf("checkMPIReduce() passed\n");
    }
  }
} 

void debugLoading(const std::shared_ptr<model::ThermomechanicalModel>& model_p) {
  // collect fixity of nodes in root rank
  // clone d_model_p->d_displacement_fixed
  const auto& loading_p = model_p->d_loading_p;
  std::vector<int> displacement_fixed_global(model_p->d_displacement_fixed.size(), 0);
  std::vector<int> displacement_fixed_local(model_p->d_displacement_fixed.size(), 0);
  for (size_t i = 0; i < model_p->d_displacement_fixed.size(); ++i) {
    bool is_fixed = false;
    for (size_t j = 0; j < 3; ++j) {
      if (!util::isDofFree(i, j, model_p->d_displacement_fixed)) {
        is_fixed = true;
        break;
      }
    }
    displacement_fixed_local[i] = is_fixed ? 1 : 0;
  }

  MPI_Reduce(displacement_fixed_local.data(), displacement_fixed_global.data(), 
      model_p->d_displacement_fixed.size(), MPI_INT, 
      MPI_SUM, 0, model_p->d_cm_p->d_comm.get());

  if (model_p->d_cm_p->d_rank == 0) {
    std::cout << "\n\n\n++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "debugging loading" << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
    if (loading_p->d_force_bcs.size() > 0) {
      std::cout << "force bc type = " << loading_p->d_force_bcs[0].d_type << std::endl;
      loading_p->d_force_bcs[0].print();
      printf("number of force nodes = %zu\n", loading_p->d_force_nodes.size());
      std::cout << std::endl << std::endl << std::endl;
    }

    if (loading_p->d_disp_bcs.size() > 0) {
      std::cout << "disp bc type = " << loading_p->d_disp_bcs[0].d_type << std::endl;
      loading_p->d_disp_bcs[0].print();
      printf("number of disp nodes = %zu\n", loading_p->d_disp_nodes[0].size());
      std::cout << std::endl << std::endl << std::endl;
      

      std::cout << "Fixed nodes: " << std::endl;
      size_t n_fixed = 0;
      for (size_t i = 0; i < displacement_fixed_global.size(); ++i) {
        // use util::isDofFree to identify fixed dof and print node id and dof
        bool is_fixed = displacement_fixed_global[i] > 0;
        if (is_fixed) {
          n_fixed++;

          // check if the node is in the disp bc region
          const libMesh::Point& xi = *model_p->d_mesh.node_ptr(i);
          auto is_in = loading_p->d_disp_bcs[0].isInRegion(xi);
          std::string is_in_str = is_in ? "true" : "false";
          // print that node is fixed and its coordinates and is_in flag
          printf("Disp BC: BC id = %d, Node = %zu, Coordinates = (%f, %f, %f), Owner = %d, Rank = %d\n", 0, i, xi(0), xi(1), xi(2), model_p->d_mesh.node_ptr(i)->processor_id(), model_p->d_cm_p->d_rank);
        }
      }
      printf("n_disp_nodes = %zu\n", n_fixed);
      std::cout << std::endl;
    }
  }
}

void debugPeriForce(const std::shared_ptr<model::ThermomechanicalModel>& model_p) {

  std::cout << "\n\n\n++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "debugging thermal peri force" << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
  
  // debug force
  std::vector<libMesh::Point> local_thermal_peri_force(model_p->d_mesh.n_nodes(), libMesh::Point(0.0, 0.0, 0.0));
  std::vector<libMesh::Point> global_thermal_peri_force(model_p->d_mesh.n_nodes(), libMesh::Point(0.0, 0.0, 0.0));

  // Loop over owned nodes
  for (size_t loc_i = 0; loc_i < model_p->d_cm_p->d_owned_size; loc_i++) {
    const auto& i = model_p->d_cm_p->d_owned_and_ghost_ids[loc_i];  // Global node ID
    const libMesh::Point& xi = *model_p->d_mesh.node_ptr(i);

    libMesh::Point force_i(0.0, 0.0, 0.0);

    // Loop over neighbors using neighbor list
    for (size_t loc_j = 0; loc_j < model_p->d_neighbor_list[i].size(); loc_j++) {
      libMesh::dof_id_type j = model_p->d_neighbor_list[i][loc_j];
      const libMesh::Point& xj = *model_p->d_mesh.node_ptr(j);
      
      // Compute reference and current bond vectors
      libMesh::Point dx = xj - xi;
      libMesh::Point du = model_p->d_displacement[j] - model_p->d_displacement[i];
      double r = dx.norm();
      
      // Get bond state
      bool fs = model_p->d_fracture[i][loc_j] == 0;  // 0 = unbroken bond

      // Compute bond strain
      double s = model_p->d_material_p->getS(dx, du);
      
      // Get force magnitude
      auto f_i = model_p->d_material_p->getThermalBondForce(r, s, fs, model_p->d_mx[i], model_p->d_theta[i], model_p->d_temperature[i]);
      auto f_j = model_p->d_material_p->getThermalBondForce(r, s, fs, model_p->d_mx[j], model_p->d_theta[j], model_p->d_temperature[j]);
      
      // Add force contribution
      force_i += (f_i + f_j) * model_p->d_neighbor_volume[i][loc_j] * model_p->d_material_p->getBondForceDirection(dx, du);
    }

    local_thermal_peri_force[i] = force_i;
  }

  model_p->debugVector(local_thermal_peri_force, "thermal peri force");
}

void setObservationPoints(std::shared_ptr<model::ThermomechanicalModel> &model_p, const geom::Cuboid &domain) {

  auto Lx = domain.d_Lx, Ly = domain.d_Ly, Lz = domain.d_Lz;
  auto center = domain.center();
  std::vector<libMesh::Point> obs_points;

  obs_points.push_back(center);
  obs_points.push_back(center + libMesh::Point(-0.5*Lx, -0.5*Ly, -0.5*Lz)); // corner 1
  obs_points.push_back(center + libMesh::Point(-0.5*Lx, 0.5*Ly, -0.5*Lz)); // corner 2
  obs_points.push_back(center + libMesh::Point(0.5*Lx, -0.5*Ly, -0.5*Lz)); // corner 3
  obs_points.push_back(center + libMesh::Point(0.5*Lx, 0.5*Ly, -0.5*Lz)); // corner 4

  if (model_p->d_mesh.mesh_dimension() == 3) {
    obs_points.push_back(center + libMesh::Point(-0.5*Lx, -0.5*Ly, 0.5*Lz)); // corner 5
    obs_points.push_back(center + libMesh::Point(-0.5*Lx, 0.5*Ly, 0.5*Lz)); // corner 6
    obs_points.push_back(center + libMesh::Point(0.5*Lx, -0.5*Ly, 0.5*Lz)); // corner 7
    obs_points.push_back(center + libMesh::Point(0.5*Lx, 0.5*Ly, 0.5*Lz)); // corner 8
  }
  
  // set observation points
  model_p->setObservationPoints(obs_points);
}