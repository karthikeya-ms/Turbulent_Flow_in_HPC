#include "StdAfx.hpp"

#include "WallhStencil.hpp"

Stencils::WallhStencil::WallhStencil(const Parameters& parameters):
FieldStencil<TurbulentFlowField>(parameters){}

void Stencils::WallhStencil::apply(TurbulentFlowField& flowField, int i, int j)
{
    RealType& value     = flowField.getWallh().getScalar(i, j);
    const int obstacle  = flowField.getFlags().getValue(i, j);

    // parameters_.meshsize->getPosX(i,j,k) exists
    // MeshSize.hpp::UniformMesh::getPosX - Line 66
    // Stencils::TurbulentVTKStencil::writePoints - Line 64
    // Stencils::VelocityStencil/ObstacleStencil::apply - using OBSTACLE
    // BoundaryType.hpp - enum for obstacle

    if ((obstacle & OBSTACLE_SELF) == 0) {    // If this is a fluid cell  
        const RealType posX     = parameters_.meshsize->getPosX(i, j) + 0.5*parameters_.meshsize->getDx(i, j);    
        const RealType posY     = parameters_.meshsize->getPosY(i, j) + 0.5*parameters_.meshsize->getDy(i, j); 
        const RealType wallX    = parameters_.geometry.lengthX;
        const RealType wallY    = parameters_.geometry.lengthY;

        RealType distX_nearest = (posX) < (wallX - posX) ? (posX) : (wallX - posX) ;
        RealType distY_nearest = (posY) < (wallY - posY) ? (posY) : (wallY - posX) ;

        if (parameters_.simulation.scenario == "pressure-channel"){
            const RealType xLimit = parameters_.bfStep.xRatio * parameters_.geometry.lengthX;
            const RealType yLimit = parameters_.bfStep.yRatio * parameters_.geometry.lengthY;

            if ((xLimit < posX) && (yLimit > posY)){
                distX_nearest =  distX_nearest < (xLimit - posX) ? distX_nearest : (xLimit - posX) ;
            }
            if ((xLimit > posX) && (yLimit < posY)){
                distY_nearest =  distY_nearest < (yLimit - posY) ? distY_nearest : (yLimit - posY) ;
            }
        }

        // Only checking normal wall distance even in BFStep case, 
        // even if the pythagorean distance to the corner can be the smallest one.
        value = distX_nearest < distY_nearest ? distX_nearest : distY_nearest ;
    }
    else { value = 0.0; }
    ASSERTION(value >= 0.0);
}

void Stencils::WallhStencil::apply(TurbulentFlowField& flowField, int i, int j, int k)
{
    RealType& value     = flowField.getWallh().getScalar(i, j, k);
    const int obstacle  = flowField.getFlags().getValue(i, j, k);
    
    if ((obstacle & OBSTACLE_SELF) == 0) {    // If this is a fluid cell
        const RealType& value2D = flowField.getWallh().getScalar(i, j);

        const RealType posZ     = parameters_.meshsize->getPosZ(i, j, k) + 0.5*parameters_.meshsize->getDz(i, j, k);
        const RealType wallZ    = parameters_.geometry.lengthZ;

        RealType distZ_nearest  = (posZ) < (wallZ - posZ) ? (posZ) : (wallZ - posZ) ;

        value = value2D < distZ_nearest ? value2D : distZ_nearest;
    }
    else { value = 0.0; }
    ASSERTION(value >= 0.0);
}