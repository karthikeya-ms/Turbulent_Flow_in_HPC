#include "StdAfx.hpp"

#include "ChViscosityInitStencil.hpp"

Stencils::ChViscosityInitStencil::ChViscosityInitStencil(const Parameters& parameters):
FieldStencil<TurbulentFlowField>(parameters){}

void Stencils::ChViscosityInitStencil::apply(TurbulentFlowField& flowField, int i, int j)
{
    flowField.getOldChVis().getScalar(i, j) = 3.0 / (parameters_.flow.Re); //3.0 / parameters_.flow.Re;//
    flowField.getTurbVisc().getScalar(i, j) = 0.22 / (parameters_.flow.Re); 
    
    if (j==2) {
        flowField.getOldChVis().getScalar(i, j-1) =  -1*flowField.getOldChVis().getScalar(i, j);
        flowField.getTurbVisc().getScalar(i, j-1) =  -1*flowField.getTurbVisc().getScalar(i, j);
    }
    if (j == parameters_.geometry.sizeY + 1) {
        flowField.getOldChVis().getScalar(i, j+1) = -1*flowField.getOldChVis().getScalar(i, j);
        flowField.getTurbVisc().getScalar(i, j-1) = -1*flowField.getTurbVisc().getScalar(i, j);
    }
    // if (i==2 && j==2){
        // spdlog::info("CheckChVisc: chVis{}", flowField.getOldChVis().getScalar(i, j));
        // throw std::runtime_error("Stop at ChViscStencil");
    // }
    
}

void Stencils::ChViscosityInitStencil::apply(TurbulentFlowField& flowField, int i, int j, int k)
{
    flowField.getOldChVis().getScalar(i, j, k) = 3.0 / parameters_.flow.Re;
}
