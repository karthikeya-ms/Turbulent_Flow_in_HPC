#include "StdAfx.hpp"

#include "ComputeLocalViscosityStencil.hpp"

#include "Definitions.hpp"
#include "TurbulentStencilFunctions.hpp"

Stencils::ComputeLocalViscosityStencil::ComputeLocalViscosityStencil(const Parameters& parameters):
FieldStencil<TurbulentFlowField>(parameters){}

void Stencils::ComputeLocalViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j) {
    RealType delta = 0.0;
    if (parameters_.turbMix.delta == 0) { //Case for no boundary layer
        delta = 0.0;
    } 
    else if (parameters_.turbMix.delta == 1) { //Case for laminar flat plate
        auto x = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i, j);
        delta  = 4.91 * x / sqrt(parameters_.walls.vectorLeft[0] * x * parameters_.flow.Re);
        //reynolds number here stand for reciprocal of viscosity
        // cause velocity and characteristic length are normalised.
    } 
    else if (parameters_.turbMix.delta == 2) {//Case for turbulent flat plate
        auto x = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i, j);
        delta  = 0.382 * x / pow(parameters_.walls.vectorLeft[0] * x * parameters_.flow.Re, 0.2);
    }

    //Calculate mixing length by Prandtl's model
    RealType lm = std::min(parameters_.turbMix.k * flowField.getWallh().getScalar(i, j), 0.09 * delta);
  
    loadLocalVelocity2D(flowField, localVelocity_, i, j);
    loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);

    RealType S11 = dudx(localVelocity_, localMeshsize_);
    RealType S22 = dvdy(localVelocity_, localMeshsize_);
    RealType S12 = 0.5 ;//* (dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_));

    RealType SijSij = S11 * S11 + S22 * S22 + 2 * S12 * S12;

    flowField.getTurbVisc().getScalar(i, j) = lm * lm * sqrt(2 * SijSij);
}

void Stencils::ComputeLocalViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
    RealType delta = 0.0;
    if (parameters_.turbMix.delta == 0) {
        delta = 0.0;
    } 
    else if (parameters_.turbMix.delta == 1) {
        auto x = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i, j);
        delta  = 4.91 * x / sqrt(parameters_.walls.vectorLeft[0] * x * parameters_.flow.Re);
    } 
    else if (parameters_.turbMix.delta == 2){
        auto x = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i, j);
        delta  = 0.382 * x / pow(parameters_.walls.vectorLeft[0] * x * parameters_.flow.Re, 0.2);
    }
    RealType lm = std::min(parameters_.turbMix.k * flowField.getWallh().getScalar(i, j, k), 0.09 * delta);

    loadLocalVelocity3D(flowField, localVelocity_, i, j, k);
    loadLocalMeshsize3D(parameters_, localMeshsize_, i, j, k);

    RealType S11 = dudx(localVelocity_, localMeshsize_);
    RealType S22 = dvdy(localVelocity_, localMeshsize_);
    RealType S33 = dwdz(localVelocity_, localMeshsize_);
    RealType S12 = 0.5 ;//* (dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_));
    RealType S13 = 0.5 ;//* (dudz(localVelocity_, localMeshsize_) + dwdx(localVelocity_, localMeshsize_));
    RealType S23 = 0.5 ;//* (dvdz(localVelocity_, localMeshsize_) + dwdy(localVelocity_, localMeshsize_));

    RealType SijSij = S11 * S11 + S22 * S22 + S33 * S33 + 2 * S12 * S12 + 2 * S13 * S13 + 2 * S23 * S23;

    flowField.getTurbVisc().getScalar(i, j, k) = lm * lm * sqrt(2 * SijSij);
}