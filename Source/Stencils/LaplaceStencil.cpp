#include "StdAfx.hpp"

#include "LaplaceStencil.hpp"

Stencils::LaplaceStencil::LaplaceStencil(const Parameters& parameters):
FieldStencil<TurbulentFlowField>(parameters){}

void Stencils::LaplaceStencil::apply(TurbulentFlowField& flowField, int i, int j)
{
    ScalarField ChVis = flowField.getChVis(); 
    const RealType* velocity = flowField.getVelocity().getVector(i, j);

    const RealType dx_0  = parameters_.meshsize->getDx(i, j);
    const RealType dx_M1 = parameters_.meshsize->getDx(i-1, j);
    const RealType dx_P1 = parameters_.meshsize->getDx(i+1, j);
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType dy_0  = parameters_.meshsize->getDy(i, j);
    const RealType dy_M1 = parameters_.meshsize->getDy(i, j-1);
    const RealType dy_P1 = parameters_.meshsize->getDy(i, j+1);
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    RealType dVidx0 = (ChVis.getScalar(i, j) - ChVis.getScalar(i-1, j)) / dx0;
    RealType dVidy0 = (ChVis.getScalar(i, j) - ChVis.getScalar(i, j-1)) / dy0;

    RealType dVidx1 = (ChVis.getScalar(i+1, j) - ChVis.getScalar(i, j)) / dx1;
    RealType dVidy1 = (ChVis.getScalar(i, j+1) - ChVis.getScalar(i, j)) / dy1;

    RealType d2Vidx2 = (dVidx1 - dVidx0) / dx_0;
    RealType d2Vidy2 = (dVidy1 - dVidy0) / dy_0;
    //ChVis location for implementation can be wrong, recheck here.

    RealType cbe = (cb2 + 1) / cb3; 
    RealType order1st =  cbe * (pow(dVidx0, 2) + pow(dVidy0, 2)) - (velocity[0] * dVidx0 + velocity[1] * dVidy0);
    RealType order2nd = (d2Vidx2 + d2Vidy2) / cb3;

    flowField.getLaplace().getScalar(i, j) = order1st + order2nd;
}

void Stencils::LaplaceStencil::apply(TurbulentFlowField& flowField, int i, int j, int k)
{
    ScalarField ChVis = flowField.getChVis(); 
    const RealType* velocity = flowField.getVelocity().getVector(i, j, k);

    const RealType dx_0  = parameters_.meshsize->getDx(i, j, k);
    const RealType dx_M1 = parameters_.meshsize->getDx(i-1, j, k);
    const RealType dx_P1 = parameters_.meshsize->getDx(i+1, j, k);
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType dy_0  = parameters_.meshsize->getDy(i, j, k);
    const RealType dy_M1 = parameters_.meshsize->getDy(i, j-1, k);
    const RealType dy_P1 = parameters_.meshsize->getDy(i, j+1, k);
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType dz_0  = parameters_.meshsize->getDz(i, j, k);
    const RealType dz_M1 = parameters_.meshsize->getDz(i, j, k-1);
    const RealType dz_P1 = parameters_.meshsize->getDz(i, j, k+1);
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);

    RealType dVidx0 = (ChVis.getScalar(i, j, k) - ChVis.getScalar(i-1, j, k)) / dx0;
    RealType dVidy0 = (ChVis.getScalar(i, j, k) - ChVis.getScalar(i, j-1, k)) / dy0;
    RealType dVidz0 = (ChVis.getScalar(i, j, k) - ChVis.getScalar(i, j, k-1)) / dz0;

    RealType dVidx1 = (ChVis.getScalar(i+1, j, k) - ChVis.getScalar(i, j, k)) / dx1;
    RealType dVidy1 = (ChVis.getScalar(i, j+1, k) - ChVis.getScalar(i, j, k)) / dy1;
    RealType dVidz1 = (ChVis.getScalar(i, j, k+1) - ChVis.getScalar(i, j, k)) / dz1;

    RealType d2Vidx2 = (dVidx1 - dVidx0) / dx_0;
    RealType d2Vidy2 = (dVidy1 - dVidy0) / dy_0;
    RealType d2Vidz2 = (dVidz1 - dVidz0) / dz_0;
    //ChVis location for implementation can be wrong, recheck here.

    RealType cbe = (cb2 + 1) / cb3; 
    RealType order1st =  cbe * (pow(dVidx0, 2) + pow(dVidy0, 2) + pow(dVidz0, 2))
                        - (velocity[0] * dVidx0 + velocity[1] * dVidy0 + velocity[2] * dVidz0);
    RealType order2nd = (d2Vidx2 + d2Vidy2 + d2Vidz2) / cb3;

    flowField.getLaplace().getScalar(i, j, k) = order1st + order2nd;
}