#include "StdAfx.hpp"

#include "VisualizeStencil.hpp"
#include "TurbulentStencilFunctions.hpp"

Stencils::VisualizeStencil::VisualizeStencil(const Parameters& parameters):
FieldStencil<TurbulentFlowField>(parameters){}

void Stencils::VisualizeStencil::apply(TurbulentFlowField& flowField, int i, int j)
{
    RealType Visc_ = 1 / parameters_.flow.Re;
    RealType ChVis = flowField.getNewChVis().getScalar(i, j);

    loadLocalVelocity2D(flowField, localVelocity_, i, parameters_.geometry.sizeY);
    loadLocalMeshsize2D(parameters_, localMeshsize_, i, parameters_.geometry.sizeY);

    RealType dudy_wall = dudy(localVelocity_, localMeshsize_);
    RealType tau_wall = std::fabs(dudy_wall) * (Visc_ + ChVis);

    loadLocalVelocity2D(flowField, localVelocity_, i, j);
    loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);
    RealType dudy_net = dudy(localVelocity_, localMeshsize_);
    RealType tau_net = std::fabs(dudy_net) * (Visc_ + ChVis);

    RealType u_tau = std::sqrt(tau_wall);   // wall shear velocity
    
    if (dudy_wall == 0) {
        flowField.getYPlus().getScalar(i,j) = 0.0;
        flowField.getUPlus().getScalar(i,j) = 0.0; 
        flowField.getTau().getScalar(i,j)   = 0.0;
    } 
    else{
        flowField.getYPlus().getScalar(i,j) = flowField.getWallh().getScalar(i,j) * u_tau / Visc_;  // y_plus
        flowField.getUPlus().getScalar(i,j) = (
            0.5 * (flowField.getVelocity().getVector(i - 1, j)[0] + flowField.getVelocity().getVector(i, j)[0])
            ) / (u_tau);   // u_plus
        flowField.getTau().getScalar(i,j) = tau_net / tau_wall; // tau_norm
    }
    
    // if (isnanf(flowField.getTau().getScalar(i,j))){
    //     spdlog::info("At i{}, j{} tau_wall{}" ,i ,j, tau_wall);
    // }
}

void Stencils::VisualizeStencil::apply(TurbulentFlowField& flowField, int i, int j, int k)
{
    RealType Visc_ = 1 / parameters_.flow.Re;
    RealType ChVis = flowField.getNewChVis().getScalar(i, j, k);

    loadLocalVelocity3D(flowField, localVelocity_, i, parameters_.geometry.sizeY, k);
    loadLocalMeshsize3D(parameters_, localMeshsize_, i, parameters_.geometry.sizeY, k);
    RealType dudy_wall = dudy(localVelocity_, localMeshsize_);
    RealType tau_wall = std::fabs(dudy_wall) * (Visc_ + ChVis);

    loadLocalVelocity3D(flowField, localVelocity_, i, j, k);
    loadLocalMeshsize3D(parameters_, localMeshsize_, i, j, k);
    RealType dudy_net = dudy(localVelocity_, localMeshsize_);
    RealType tau_net = std::fabs(dudy_net) * (Visc_ + ChVis);

    RealType u_tau = std::sqrt(tau_wall);   // wall shear velocity
    
    if (dudy_wall == 0) {
        flowField.getYPlus().getScalar(i,j) = 0.0;
        flowField.getUPlus().getScalar(i,j) = 0.0; 
        flowField.getTau().getScalar(i,j)   = 0.0;
    } 
    else {
        int k_half = parameters_.geometry.sizeZ/2 - 1;
        flowField.getYPlus().getScalar(i, j, k) = flowField.getWallh().getScalar(i, j, k_half) * u_tau / Visc_;  // y_plus
        flowField.getUPlus().getScalar(i, j, k) = (
            0.5 * (flowField.getVelocity().getVector(i - 1, j, k)[0] + flowField.getVelocity().getVector(i, j, k)[0])
            ) / (u_tau);   // u_plus
        flowField.getTau().getScalar(i, j, k) = tau_net / tau_wall; // tau_norm
    }
}
