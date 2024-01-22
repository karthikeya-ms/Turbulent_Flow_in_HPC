#include "StdAfx.hpp"

#include "QStencil.hpp"
#include "TurbulentStencilFunctions.hpp"

Stencils::QStencil::QStencil(const Parameters& parameters):
FieldStencil<TurbulentFlowField>(parameters){}

void Stencils::QStencil::apply(TurbulentFlowField& flowField, int i, int j)
{
    
    const int obstacle  = flowField.getFlags().getValue(i, j);
    if ((obstacle & OBSTACLE_SELF) == 0) { // If the cell is fluid
        const RealType d = flowField.getWallh().getScalar(i, j) + 1e-6;
        const RealType ChVis = flowField.getOldChVis().getScalar(i, j);

        loadLocalVelocity2D(flowField, localVelocity_, i, j);
        loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);

        // const RealType S11 = dudx(localVelocity_, localMeshsize_);
        // const RealType S22 = dvdy(localVelocity_, localMeshsize_);
        const RealType S12 = 0.5 * (dudy(localVelocity_, localMeshsize_) - dvdx(localVelocity_, localMeshsize_));
        // const RealType SijSij = S11 * S11 + S22 * S22 + 2 * S12 * S12;
        // RealType dudy = (flowField.getVelocity().getVector(i, j + 1)[0] - flowField.getVelocity().getVector(i, j - 1)[0])
        //               / (0.5*parameters_.meshsize->getDy(i, j-1) + parameters_.meshsize->getDy(i,j) + 0.5*parameters_.meshsize->getDy(i,j+1));

        // RealType dvdx = (flowField.getVelocity().getVector(i + 1, j)[1] - flowField.getVelocity().getVector(i - 1, j)[1])
        //               / (0.5*parameters_.meshsize->getDx(i-1, j) + parameters_.meshsize->getDx(i,j) + 0.5*parameters_.meshsize->getDx(i+1,j));

        // const RealType S12 = 0.5 * (dudy - dvdx);

        const RealType S = fabs(S12) * sqrt(2.0); // vorticity

        const RealType X = ChVis * parameters_.flow.Re;
        const RealType fv1 = pow(X, 3.0) / (pow(X, 3.0) + pow(cv1, 3.0));
        const RealType fv2 = 1.0 - X / (1.0 + X * fv1);
        const RealType S_ = fmax(S + fv2 * ChVis / pow(_k*d, 2.0), 0.3 * S); // Openfoam

        const RealType ft2 = ct3 * exp(ct4 * pow(X, 2.0) * (-1.0));

        const RealType r = fmin(10.0, ChVis / (fmax(S_ , 1e-6) * pow(_k*d, 2.0))); // Openfoam
        const RealType g = r + cw2 * (pow(r, 6.0) - r);
        const RealType fw = g * pow(((1.0 + pow(cw3, 6.0)) / (pow(g, 6.0) + pow(cw3, 6.0))), 1.0/6.0);
    
        flowField.getQ().getScalar(i,j) = cb1 * S_ * ChVis * (1.0 - ft2) 
                                        + (cb1 * ft2 / pow(_k, 2.0) - cw1 * fw) * pow(ChVis / d, 2.0);

        
        // spdlog::info("CheckQ: d{} X{} fv1{} fv2{} S_{} ft2{} r{} g{} fw{}", d, X, fv1, fv2, S_, ft2, r, g, fw);
        // throw std::runtime_error("Stop at QStencil");
    }
}

void Stencils::QStencil::apply(TurbulentFlowField& flowField, int i, int j, int k)
{
    
    const int obstacle  = flowField.getFlags().getValue(i, j, k);
    if ((obstacle & OBSTACLE_SELF) == 0) { // If the cell is fluid
        const RealType d = flowField.getWallh().getScalar(i, j, k);
        const RealType ChVis = flowField.getOldChVis().getScalar(i, j, k);

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
        const RealType S_ = S + fv2 * ChVis / (pow(_k, 2) * pow(d, 2));

        const RealType ft2 = ct3 * exp(ct4 * pow(X, 2) * (-1.0));

        const RealType r = ChVis / (S_ * pow(_k, 2) * pow(d, 2) + EPSILON);
        const RealType g = r + cw2 * (pow(r, 6) - r);
        const RealType fw = g * pow(((1 + pow(cw3, 6)) / (pow(g, 6) + pow(cw3, 6))), 0.1667);
    
        flowField.getQ().getScalar(i,j,k) = cb1 * S_ * ChVis * (1 - ft2) 
                                        + (cb1 * ft2 / pow(_k, 2) - cw1 * fw) * pow(ChVis / d, 2);

        // spdlog::info("CheckQ: d{} X{} fv1{} fv2{} S_{} ft2{} r{} g{} fw{}", d, X, fv1, fv2, S_, ft2, r, g, fw);
        // throw std::runtime_error("Stop at QStencil");
    }
}
