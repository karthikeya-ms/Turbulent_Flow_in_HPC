#include "StdAfx.hpp"

#include "QStencil.hpp"
#include "TurbulentStencilFunctions.hpp"

Stencils::QStencil::QStencil(const Parameters& parameters):
FieldStencil<TurbulentFlowField>(parameters){}

void Stencils::QStencil::apply(TurbulentFlowField& flowField, int i, int j)
{
    const RealType d = flowField.getWallh().getScalar(i, j);
    const RealType ChVis = flowField.getChVis().getScalar(i, j);

    loadLocalVelocity2D(flowField, localVelocity_, i, j);
    loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);

    const RealType S11 = dudx(localVelocity_, localMeshsize_);
    const RealType S22 = dvdy(localVelocity_, localMeshsize_);
    const RealType S12 = 0.5 * (dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_));
    const RealType SijSij = S11 * S11 + S22 * S22 + 2 * S12 * S12;
    const RealType S = sqrt(2 * SijSij);

    const RealType X = ChVis * parameters_.flow.Re;
    const RealType fv1 = pow(X, 3) / (pow(X, 3) + pow(cv1, 3));
    const RealType fv2 = 1 - X / (1 + X * fv1);
    const RealType S_ = S + fv2 * ChVis / (pow(k, 2) * pow(d, 2));

    const RealType ft2 = ct3 * exp(ct4 * pow(X, 2) * (-1.0));

    const RealType r = ChVis / (S_ * pow(k, 2) * pow(d, 2));
    const RealType g = r + cw2 * (pow(r, 6) - r);
    const RealType fw = g * pow(((1 + pow(cw3, 6)) / (pow(g, 6) + pow(cw3, 6))), 0.1667);
    
    flowField.getQ().getScalar(i,j) = cb1 * S_ * ChVis * (1 - ft2) 
                                    + (cb1 * ft2 / pow(k, 2) - cw1 * fw) * pow(ChVis / d, 2);
}

void Stencils::QStencil::apply(TurbulentFlowField& flowField, int i, int j, int k)
{
    const RealType d = flowField.getWallh().getScalar(i, j, k);
    const RealType ChVis = flowField.getChVis().getScalar(i, j, k);

    loadLocalVelocity3D(flowField, localVelocity_, i, j, k);
    loadLocalMeshsize3D(parameters_, localMeshsize_, i, j, k);

    const RealType S11 = dudx(localVelocity_, localMeshsize_);
    const RealType S22 = dvdy(localVelocity_, localMeshsize_);
    const RealType S33 = dwdz(localVelocity_, localMeshsize_);
    const RealType S12 = 0.5 * (dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_));
    const RealType S13 = 0.5 * (dudz(localVelocity_, localMeshsize_) + dwdx(localVelocity_, localMeshsize_));
    const RealType S23 = 0.5 * (dvdz(localVelocity_, localMeshsize_) + dwdy(localVelocity_, localMeshsize_));
    const RealType SijSij = S11 * S11 + S22 * S22 + S33 * S33 + 2 * S12 * S12 + 2 * S13 * S13 + 2 * S23 * S23;
    const RealType S = sqrt(2 * SijSij);

    const RealType X = ChVis * parameters_.flow.Re;
    const RealType fv1 = pow(X, 3) / (pow(X, 3) + pow(cv1, 3));
    const RealType fv2 = 1 - X / (1 + X * fv1);
    const RealType S_ = S + fv2 * ChVis / (pow(k, 2) * pow(d, 2));

    const RealType ft2 = ct3 * exp(ct4 * pow(X, 2) * (-1.0));

    const RealType r = ChVis / (S_ * pow(k, 2) * pow(d, 2));
    const RealType g = r + cw2 * (pow(r, 6) - r);
    const RealType fw = g * pow(((1 + pow(cw3, 6)) / (pow(g, 6) + pow(cw3, 6))), 0.1667);
    
    flowField.getQ().getScalar(i,j,k) = cb1 * S_ * ChVis * (1 - ft2) 
                                    + (cb1 * ft2 / pow(k, 2) - cw1 * fw) * pow(ChVis / d, 2);
}
