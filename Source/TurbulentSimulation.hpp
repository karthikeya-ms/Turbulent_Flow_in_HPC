#pragma once

#include "Definitions.hpp"
#include "TurbulentFlowField.hpp"
#include "GlobalBoundaryFactory.hpp"
#include "Iterators.hpp"

#include "Simulation.hpp"
#include "Stencils/ComputeLocalViscosityStencil.hpp"
#include "Stencils/WallhStencil.hpp"
#include "Stencils/MaxTurbViscStencil.hpp"
#include "Stencils/TurbulentFGHStencil.hpp"
#include "Stencils/TurbulentVTKStencil.hpp" // turbulent visualisation

#include "Stencils/ChViscosityStencil.hpp" // Spalart-Allmaras turbulence model stencils
#include "Stencils/QStencil.hpp"
#include "Stencils/NablaStencil.hpp"

class TurbulentSimulation : public Simulation {
protected:
  TurbulentFlowField& turbFlowField_;

  Stencils::WallhStencil                  wallhStencil_;
  FieldIterator<TurbulentFlowField>       wallhIterator_;
  Stencils::ComputeLocalViscosityStencil  turbViscStencil_;
  FieldIterator<TurbulentFlowField>       turbViscIterator_;
  Stencils::MaxTurbViscStencil            maxTurbViscStencil_;
  FieldIterator<TurbulentFlowField>       maxTurbViscIterator_;
  Stencils::TurbulentFGHStencil           turbFGHStencil_;
  FieldIterator<TurbulentFlowField>       turbFGHIterator_;

  Stencils::ChViscosityStencil            ChViscosityStencil_;
  FieldIterator<TurbulentFlowField>       ChViscosityIterator_;
  Stencils::QStencil                      QStencil_;
  FieldIterator<TurbulentFlowField>       QIterator_;
  Stencils::NablaStencil                  NablaStencil_;
  FieldIterator<TurbulentFlowField>       NablaIterator_;


  // Set up the boundary conditions
  // GlobalBoundaryIterator<TurbulentFlowField>  turbWallFGHIterator_;
  // declare other stencils / iterators that work on turbFlowField here (WS2)

  void setTimeStep();
public:
  TurbulentSimulation(Parameters& parameters, TurbulentFlowField& flowField);
  ~TurbulentSimulation() = default;
  /** Initialises the flow field according to the scenario */
  void initializeFlowField();

  void solveTimestep(); // should use override? virtual? (WS2?)

  /** Plots the flow field */
  void plotVTK(int timeStep, RealType simulationTime);

  ParallelManagers::PetscTurbulentParallelManager parallel_manager_;
};
