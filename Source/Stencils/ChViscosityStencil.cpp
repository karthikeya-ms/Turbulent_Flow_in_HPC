#include "StdAfx.hpp"

#include "ChViscosityStencil.hpp"

Stencils::ChViscosityStencil::ChViscosityStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {}

void Stencils::ChViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j) {
    const RealType dt       = parameters_.timestep.dt;
    const int      obstacle = flowField.getFlags().getValue(i, j);

  // += stands for ChVis(n+1) = ChVis(n) + additional terms
  if ((obstacle & OBSTACLE_SELF) == 0) {    // If this is a fluid cell
    RealType ChVisc = flowField.getOldChVis().getScalar(i, j);
    // if (ChVisc < EPSILON) { 
    //   spdlog::info("At i{}, j{}", i, j);
    //   throw std::runtime_error("ChVisc at ComputeLocalViscStencil"); 
    // }
    flowField.getNewChVis().getScalar(i, j) = fmax( 0.0,
      ChVisc + dt * (flowField.getQ().getScalar(i,j) + flowField.getNabla().getScalar(i,j))
    );
  }
  // else{ throw std::runtime_error("Reaching out of domain at ComputeLocalViscStencil"); } 
  // Never executed for channel scenario.
}

void Stencils::ChViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
    const RealType dt       = parameters_.timestep.dt;
    const int      obstacle = flowField.getFlags().getValue(i, j, k);

  // += stands for ChVis(n+1) = ChVis(n) + additional terms
  if ((obstacle & OBSTACLE_SELF) == 0) {// If this is a fluid cell
    RealType ChVisc = flowField.getOldChVis().getScalar(i, j, k);
    flowField.getNewChVis().getScalar(i, j, k) = fmax( EPSILON,
      ChVisc + dt * (flowField.getQ().getScalar(i,j,k) + flowField.getNabla().getScalar(i,j,k))
    );
  }
}
