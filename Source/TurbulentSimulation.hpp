#pragma once

#include "Definitions.hpp"
#include "TurbulentFlowField.hpp"
#include "GlobalBoundaryFactory.hpp"
#include "Iterators.hpp"

#include "Simulation.hpp"
#include "Stencils/ComputeLocalViscosityStencil.hpp"
#include "Stencils/WallhStencil.hpp"
#include "Stencils/TurbulentVTKStencil.hpp" // turbulent visualisation

class TurbulentSimulation : public Simulation {
protected:
  TurbulentFlowField& turbFlowField_;

  Stencils::ComputeLocalViscosityStencil  turbViscStencil_;
  FieldIterator<TurbulentFlowField>       turbViscIterator_;
  Stencils::WallhStencil                  wallhStencil_;
  FieldIterator<TurbulentFlowField>       wallhIterator_;
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
};
