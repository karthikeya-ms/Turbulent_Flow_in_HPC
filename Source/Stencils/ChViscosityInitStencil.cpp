#include "StdAfx.hpp"

#include "ChViscosityInitStencil.hpp"

Stencils::ChViscosityInitStencil::ChViscosityInitStencil(const Parameters& parameters):
FieldStencil<TurbulentFlowField>(parameters){}

void Stencils::ChViscosityInitStencil::apply(TurbulentFlowField& flowField, int i, int j)
{
    flowField.getChVis().getScalar(i, j) = 1/(2 * parameters_.flow.Re);
    
    // spdlog::info("CheckChVisc: chVis{}", flowField.getChVis().getScalar(i, j));
    // throw std::runtime_error("Stop at ChViscStencil");
}

void Stencils::ChViscosityInitStencil::apply(TurbulentFlowField& flowField, int i, int j, int k)
{
    flowField.getChVis().getScalar(i, j, k) = 1/(2 * parameters_.flow.Re);
}
