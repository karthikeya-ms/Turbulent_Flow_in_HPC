#include "StdAfx.hpp"

#include "ChViscosityStencil.hpp"

Stencils::ChViscosityStencil::ChViscosityStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {}

void Stencils::ChViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j) {
    const RealType dt       = parameters_.timestep.dt;
    const int      obstacle = flowField.getFlags().getValue(i, j);
    // ScalarField&   ChVis    = flowField.getChVis();

  // += stands for ChVis(n+1) = ChVis(n) + additional terms
  if ((obstacle & OBSTACLE_SELF) == 0) {    // If this is a fluid cell
    // if ((obstacle & OBSTACLE_RIGHT) == 0 && (obstacle & OBSTACLE_TOP) == 0) { 
      // Check whether the neighbor is also fluid
    
    // spdlog::info(" before: chVis{}", flowField.getChVis().getScalar(i, j));
    flowField.getChVis().getScalar(i, j) += dt * (flowField.getQ().getScalar(i,j) + flowField.getNabla().getScalar(i,j));
    // spdlog::info(" after: chVis{}", flowField.getChVis().getScalar(i, j));

    // throw std::runtime_error("Stop at ChViscStencil");


  }
  // else{
  //   throw std::runtime_error("Reaching out of domain at ComputeLocalViscStencil");
  // }
}

void Stencils::ChViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
    const RealType dt       = parameters_.timestep.dt;
    const int      obstacle = flowField.getFlags().getValue(i, j, k);
    // ScalarField&   ChVis    = flowField.getChVis();

  // += stands for ChVis(n+1) = ChVis(n) + additional terms
  if ((obstacle & OBSTACLE_SELF) == 0) {// If this is a fluid cell
    if ((obstacle & OBSTACLE_RIGHT) == 0 && (obstacle & OBSTACLE_TOP) == 0) { 
      // Check whether the neighbor is also fluid
      flowField.getChVis().getScalar(i, j, k) += dt * (flowField.getQ().getScalar(i,j,k) + flowField.getNabla().getScalar(i,j,k));
    } else { // Otherwise, set to zero.
      flowField.getChVis().getScalar(i, j, k) = 0;
    }
  }
}
