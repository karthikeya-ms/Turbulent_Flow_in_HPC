#include "StdAfx.hpp"

#include "TurbulentSimulation.hpp"


TurbulentSimulation::TurbulentSimulation(Parameters& parameters, TurbulentFlowField& flowField):
  Simulation::Simulation(parameters, flowField),
  turbFlowField_(flowField),
  // and stencils / iterators using turbFlowField (WS2)
  wallhStencil_(parameters),
  wallhIterator_(turbFlowField_, parameters, wallhStencil_, 1, 0), // Thanks to Field Iterator Implementation
  turbViscStencil_(parameters),
  turbViscIterator_(turbFlowField_, parameters, turbViscStencil_, 1, 0), // keeping in mind the FGH implementation
  // not .., 0, 1) as we use turbVisc at i,j,k==2 for i,j,k==1 implicitly in
  // implementation because did not understand how to use globalBoundaryIterator
  maxTurbViscStencil_(parameters),
  maxTurbViscIterator_(turbFlowField_, parameters, maxTurbViscStencil_, 1, 0), // Thanks to Field Iterator Implementation
  turbFGHStencil_(parameters),
  turbFGHIterator_(turbFlowField_, parameters, turbFGHStencil_),

  // Stencils and iterators for spalart allmaras model
  ChViscosityInitStencil_(parameters),
  ChViscosityInitIterator_(turbFlowField_, parameters, ChViscosityInitStencil_, 1, 0), // Inside domain channel
  ChViscosityStencil_(parameters),
  ChViscosityIterator_(turbFlowField_, parameters, ChViscosityStencil_, 1, 0),
  QStencil_(parameters),
  QIterator_(turbFlowField_, parameters, QStencil_, 1, 0),
  NablaStencil_(parameters),
  NablaIterator_(turbFlowField_, parameters, NablaStencil_, 1, 0),
  parallel_manager_(parameters, flowField)
{
}

void TurbulentSimulation::initializeFlowField() {

  Simulation::initializeFlowField(); // Same init as Simulation
  wallhIterator_.iterate(); // Calculation needed only once, hence here.
  // ChViscosityInitIterator_.iterate();

  spdlog::info("turbModel: {}", parameters_.simulation.turbModel);
  if (parameters_.simulation.turbModel == "turbMixing"){
    spdlog::info("MixingOpt: {}", parameters_.turbMix.delta);
  }
  spdlog::info("scenario: {}", parameters_.simulation.scenario);
}

void TurbulentSimulation::solveTimestep() {

  // others maybe added (WS2)
  
  // Can use fghIterator_ and others as these memebres were "protected" in Simulation class. 
  // Such protected members can be directly accessed in the child class.

  // Determine and set max. timestep which is allowed in this simulation
  parallel_manager_.communicateVelocity();
  setTimeStep();
  // Compute FGH
  turbFGHIterator_.iterate();
  // Set global boundary values
  wallFGHIterator_.iterate();
  // TODO WS1: compute the right hand side (RHS)
  rhsIterator_.iterate();
  // Solve for pressure
  solver_->solve();
  parallel_manager_.communicatePressure();
  // TODO WS2: communicate pressure values
  // Compute velocity
  velocityIterator_.iterate();
  obstacleIterator_.iterate();
  parallel_manager_.communicateVelocity();
  // TODO WS2: communicate velocity values
  // Iterate for velocities on the boundary
  wallVelocityIterator_.iterate();
  
  if (parameters_.simulation.turbModel == "turbSA") {
    parallel_manager_.communicateChViscosity();
    //Compute source terms in Spalart-Allmaras model
    QIterator_.iterate();
    //Compute Nabla terms in Spalart-Allmaras model
    NablaIterator_.iterate();
    //Compute characteristic viscosity
    ChViscosityIterator_.iterate();
  }

  // Compute local viscosities
  turbViscIterator_.iterate();

  if (parameters_.simulation.turbModel != "turbSA") {
    parallel_manager_.communicateChViscosity();
  }

  if (parameters_.simulation.turbModel == "turbSA") {  
    // Update current ChVisc values
    turbFlowField_.updateChVis();
  }

}

void TurbulentSimulation::plotVTK(int timeStep, RealType simulationTime) {
  Stencils::VisualizeStencil         visualizeStencil(parameters_);
  FieldIterator<TurbulentFlowField>  visualizeIterator(turbFlowField_, parameters_, visualizeStencil, 1, 0);
  Stencils::TurbulentVTKStencil     vtkStencil(parameters_);
  FieldIterator<TurbulentFlowField> vtkIterator(turbFlowField_, parameters_, vtkStencil, 1, 0);
  
  visualizeIterator.iterate();
  vtkIterator.iterate();
  vtkStencil.write(turbFlowField_, timeStep, simulationTime);
}

void TurbulentSimulation::setTimeStep() {

  RealType localMin, globalMin, maxLocalVisc;
  ASSERTION(parameters_.geometry.dim == 2 || parameters_.geometry.dim == 3);
  RealType factor = 1.0 / (parameters_.meshsize->getDxMin() * parameters_.meshsize->getDxMin())
                    + 1.0 / (parameters_.meshsize->getDyMin() * parameters_.meshsize->getDyMin());

  // Determine maximum velocity
  maxUStencil_.reset();
  maxUFieldIterator_.iterate();
  maxUBoundaryIterator_.iterate();

  // Determine maximum turbVisc
  maxTurbViscStencil_.reset();
  maxTurbViscIterator_.iterate();

  //RealType u_mini = maxUStencil_.getMaxValues()[0] < 1e-12 ? 1e-12 : maxUStencil_.getMaxValues()[0]; commented after feedback from WS1
  //RealType v_mini = maxUStencil_.getMaxValues()[1] < 1e-12 ? 1e-12 : maxUStencil_.getMaxValues()[1];
  
  if (parameters_.geometry.dim == 3) {
    //RealType w_mini = maxUStencil_.getMaxValues()[2] < 1e-12 ? 1e-12 : maxUStencil_.getMaxValues()[2]; commented after feedback from WS1
    factor += 1.0 / (parameters_.meshsize->getDzMin() * parameters_.meshsize->getDzMin());
    //parameters_.timestep.dt = 1.0 / (w_mini); commented after feedback from WS1
      parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[2] + EPSILON);
  } else {
    //parameters_.timestep.dt = 1.0 / (u_mini);
      parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[0] + EPSILON);
  }

  // localMin = std::min(parameters_.timestep.dt, std::min(std::min(parameters_.flow.Re/(2 * factor), 1.0 /
  // maxUStencil_.getMaxValues()[0]), 1.0 / maxUStencil_.getMaxValues()[1]));

  maxLocalVisc = 1.0 / parameters_.flow.Re + maxTurbViscStencil_.getMaxValues();

  localMin = std::min(
     1.0 / (2 * maxLocalVisc * factor),
    std::min(
      parameters_.timestep.dt, std::min(1 / (maxUStencil_.getMaxValues()[0] + EPSILON), 1 / (maxUStencil_.getMaxValues()[1] +EPSILON))
    )
  );

  // Here, we select the type of operation before compiling. This allows to use the correct
  // data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
  // machines.

  globalMin = MY_FLOAT_MAX;
  MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

  parameters_.timestep.dt = globalMin;
  parameters_.timestep.dt *= parameters_.timestep.tau;
}
