
#include "geometry.h"
#include "thermomechanical_model.h"
#include "libmesh_includes.h"

#include "bc.h"
#include "heat_source.h"
#include "io.h"
#include "loading.h"
#include "material.h"
#include "mpi.h"

void applyTemperatureBoundaryConditions(libMesh::EquationSystems &equation_systems, const inp::MaterialDeck &material_deck, const double &T0, const double & T1);

void setDisplacementAndForceConditions(std::shared_ptr<model::ThermomechanicalModel> &model_p, const geom::Cuboid &domain, const double &tFinal, const double &load_rate);

void setObservationPoints(std::shared_ptr<model::ThermomechanicalModel> &model_p, const geom::Cuboid &domain);

int main(int argc, char** argv) {
  // Initialize libMesh
  libMesh::LibMeshInit init(argc, argv);
  util::io::setCommunicator(init.comm().rank());

  // Problem parameters
  int dim = 2;
  double Lx = 0.04, Ly = 0.01, Lz = 0.0;
  if (dim == 3) {
    Lz = 0.01;
  }
  libMesh::Point center(0.5*Lx, 0.5*Ly, 0.5*Lz);
  geom::Cuboid domain(Lx, Ly, Lz, center);
  const unsigned int nx = 160, ny = 40, nz = 40; // Number of elements in each direction
  
  // Setup input deck
  inp::MaterialDeck material_deck;
  material_deck.setDefaults();
  material_deck.d_dim = dim;
  material_deck.d_horizon = 3*(domain.d_Lx/nx);
  // material_deck.d_alpha = 0.0;

  double dt = 2e-4;
  double tFinal = 1.0;
  int nsteps = tFinal/dt;
  int write_interval = nsteps/100;
  
  const double T0 = 293.0; // Initial temperature in K
  const double T1 = 293.0; // Dirichlet BC at x = L (if no robinBC)
  const double load_rate = 0.0; // Load rate in N/s

  // Create a Gaussian heat source at the center of the domain
  std::string sfn_type = "gaussian";
  std::vector<double> sfn_params = {domain.d_Ly/10, 1.0, domain.d_x(0), domain.d_x(1), domain.d_x(2)};
  // std::string tfn_type = "linear_step_const_value";
  // std::vector<double> tfn_params = {1000.0/tFinal, 0.0, 0.2*tFinal, tFinal};
  std::string tfn_type = "sin";
  std::vector<double> tfn_params = {4000.0/tFinal, 2/tFinal};
  auto heat_source_p = std::make_unique<loading::HeatSource>(sfn_type, sfn_params, tfn_type, tfn_params);

  // Create mesh
  libMesh::Mesh mesh(init.comm());
  if (dim == 2) {
    libMesh::MeshTools::Generation::build_square(mesh, nx, ny, 0., domain.d_Lx, 0., domain.d_Ly, libMesh::QUAD4);
  } else {
    libMesh::MeshTools::Generation::build_cube(mesh, nx, ny, nz, 0., domain.d_Lx, 0., domain.d_Ly, 0., domain.d_Lz, libMesh::HEX8);
  }

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

    if (step % write_interval == 0)
      model_p->write(step/write_interval + 1, true);
  }

  return 0;
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
    {
      auto& bc = loading_p->d_disp_bcs[0];
      bc.d_type = "Displacement_BC";
      auto left_face_center = domain.d_x - libMesh::Point(0.5*domain.d_Lx, 0.0, 0.0);
      bc.d_region_p = std::make_shared<geom::Cuboid>(horizon, domain.d_Ly, domain.d_Lz, left_face_center + libMesh::Point(0.5*horizon, 0.0, 0.0));
      bc.d_direction = {0, 1, 2};
      bc.d_is_zero = true;
    }

    if (loading_p->d_disp_bcs.size() > 1) {
      auto& bc = loading_p->d_disp_bcs[1];
      bc.d_type = "Displacement_BC";
      auto right_face_center = domain.d_x + libMesh::Point(0.5*domain.d_Lx, 0.0, 0.0);
      bc.d_region_p = std::make_shared<geom::Cuboid>(horizon, domain.d_Ly, domain.d_Lz, right_face_center - libMesh::Point(0.5*horizon, 0.0, 0.0));
      bc.d_direction = {0, 1, 2};
      bc.d_is_zero = true;
    }
  }
  
  // force BC: volume of thickness horizon from x = L - horizon to x = L
  if (std::abs(load_rate) > 1e-10) {
    loading_p->d_force_bcs.resize(1);
    auto& bc = loading_p->d_force_bcs[0];
    bc.d_type = "Force_BC";
    auto right_face_center = domain.d_x + libMesh::Point(0.5*domain.d_Lx, 0.0, 0.0);
    bc.d_region_p = std::make_shared<geom::Cuboid>(horizon, domain.d_Ly, domain.d_Lz, right_face_center - libMesh::Point(0.5*horizon, 0.0, 0.0));
    bc.d_direction = {1};
    // bc.d_time_fn_type = "linear_step_const_value";
    // bc.d_time_fn_params = {load_rate, 0.0, 0.25*tFinal, tFinal};
    bc.d_time_fn_type = "sin";
    bc.d_time_fn_params = {load_rate, 4.0/tFinal}; // 4 cycles in [0, tFinal]
  }
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