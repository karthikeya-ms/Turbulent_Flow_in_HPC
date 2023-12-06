#include "StdAfx.hpp"

#include "ComputeLocalViscosityStencil.hpp"

Stencils::ComputeLocalViscosityStencil::ComputeLocalViscosityStencil(const Parameters& parameters):
FieldStencil<TurbulentFlowField>(parameters){}

void Stencils::ComputeLocalViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j)
{
    RealType& value = flowField.getTurbVisc().getScalar(i, j);
    // const RealType*  wallh = flowField.getWallh().getScalar(i, j);

    // Refer Velocity, FGH, RHS and others stencil
    // auto dx = parameters_.meshsize->getDx(i,j);
    // auto dy = parameters_.meshsize->getDy(i,j);
    // auto dt = parameters_.timestep.dt;

    value = 1;
}

void Stencils::ComputeLocalViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j, int k)
{
    RealType& value = flowField.getTurbVisc().getScalar(i, j, k);
    // const RealType*  wallh = flowField.getWallh().getScalar(i, j);

    // Refer Velocity, FGH, RHS and others stencil
    // auto dx = parameters_.meshsize->getDx(i, j, k);
    // auto dy = parameters_.meshsize->getDy(i, j, k);
    // auto dz = parameters_.meshsize->getDz(i, j, k);
    // auto dt = parameters_.timestep.dt;

    value = 1;
}