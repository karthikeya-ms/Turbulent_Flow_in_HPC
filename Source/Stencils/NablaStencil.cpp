#include "StdAfx.hpp"

#include "NablaStencil.hpp"
#include "TurbulentStencilFunctions.hpp"

Stencils::NablaStencil::NablaStencil(const Parameters& parameters):
FieldStencil<TurbulentFlowField>(parameters){}

void Stencils::NablaStencil::apply(TurbulentFlowField& flowField, int i, int j)
{
    ScalarField& ChVis = flowField.getChVis(); 
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
    //SUMMARY: AddressSanitizer: heap-use-after-free
    RealType dVidy0 = (ChVis.getScalar(i, j) - ChVis.getScalar(i, j-1)) / dy0;

    RealType dVidx1 = (ChVis.getScalar(i+1, j) - ChVis.getScalar(i, j)) / dx1;
    RealType dVidy1 = (ChVis.getScalar(i, j+1) - ChVis.getScalar(i, j)) / dy1;

    RealType Vi_R = 0.5 * (ChVis.getScalar(i, j) * dx_P1 + ChVis.getScalar(i+1, j) * dx_0) / dx1;
    RealType Vi_L = 0.5 * (ChVis.getScalar(i, j) * dx_M1 + ChVis.getScalar(i-1, j) * dx_0) / dx0;
    RealType Vi_U = 0.5 * (ChVis.getScalar(i, j) * dy_P1 + ChVis.getScalar(i, j+1) * dy_0) / dy1;
    RealType Vi_D = 0.5 * (ChVis.getScalar(i, j) * dy_M1 + ChVis.getScalar(i, j-1) * dy_0) / dy0;

    RealType d2Vidx2 = (dVidx1 * Vi_R - dVidx0 * Vi_L) / dx_0;
    RealType d2Vidy2 = (dVidy1 * Vi_U - dVidy0 * Vi_D) / dy_0;
    //ChVis gradients can be optimised with increased accuray requirement.

    RealType cbe = cb2 / cb3; 
    RealType order1st =  cbe * (pow(dVidx0, 2) + pow(dVidy0, 2)) - (velocity[0] * dVidx0 + velocity[1] * dVidy0);
    RealType order2nd = (d2Vidx2 + d2Vidy2) / cb3;

    flowField.getNabla().getScalar(i, j) = order1st + order2nd;
}

void Stencils::NablaStencil::apply(TurbulentFlowField& flowField, int i, int j, int k)
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

    RealType Vi_R = 0.5 * (ChVis.getScalar(i, j, k) * dx_P1 + ChVis.getScalar(i+1, j, k) * dx_0) / dx1 + 1 / parameters_.flow.Re;
    RealType Vi_L = 0.5 * (ChVis.getScalar(i, j, k) * dx_M1 + ChVis.getScalar(i-1, j, k) * dx_0) / dx0 + 1 / parameters_.flow.Re;
    RealType Vi_U = 0.5 * (ChVis.getScalar(i, j, k) * dy_P1 + ChVis.getScalar(i, j+1, k) * dy_0) / dy1 + 1 / parameters_.flow.Re;
    RealType Vi_D = 0.5 * (ChVis.getScalar(i, j, k) * dy_M1 + ChVis.getScalar(i, j-1, k) * dy_0) / dy0 + 1 / parameters_.flow.Re;
    RealType Vi_F = 0.5 * (ChVis.getScalar(i, j, k) * dz_P1 + ChVis.getScalar(i, j, k+1) * dz_0) / dz1 + 1 / parameters_.flow.Re;
    RealType Vi_B = 0.5 * (ChVis.getScalar(i, j, k) * dz_M1 + ChVis.getScalar(i, j, k-1) * dz_0) / dz0 + 1 / parameters_.flow.Re;

    RealType d2Vidx2 = (dVidx1 * Vi_R - dVidx0 * Vi_L) / dx_0;
    RealType d2Vidy2 = (dVidy1 * Vi_U - dVidy0 * Vi_D) / dy_0;
    RealType d2Vidz2 = (dVidz1 * Vi_F - dVidz0 * Vi_B) / dz_0;
    //ChVis gradients can be optimised with increased accuray requirement.

    RealType cbe = cb2 / cb3; 
    RealType order1st =  cbe * (pow(dVidx0, 2) + pow(dVidy0, 2) + pow(dVidz0, 2))
                        - (velocity[0] * dVidx0 + velocity[1] * dVidy0 + velocity[2] * dVidz0);
    RealType order2nd = (d2Vidx2 + d2Vidy2 + d2Vidz2) / cb3;

    flowField.getNabla().getScalar(i, j, k) = order1st + order2nd;
}