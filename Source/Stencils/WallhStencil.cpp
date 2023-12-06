#include "StdAfx.hpp"

#include "WallhStencil.hpp"

Stencils::WallhStencil::WallhStencil(const Parameters& parameters):
FieldStencil<TurbulentFlowField>(parameters){}

void Stencils::WallhStencil::apply(TurbulentFlowField& flowField, int i, int j)
{
    RealType& value     = flowField.getWallh().getScalar(i, j);
    const int obstacle  = flowField.getFlags().getValue(i, j);

    const RealType posX     = parameters_.meshsize->getPosX(i, j) 
                            + 0.5*parameters_.meshsize->getDx(i, j);    
    const RealType posY     = parameters_.meshsize->getPosY(i, j) 
                            + 0.5*parameters_.meshsize->getDy(i, j); 

    // parameters_.meshsize->getPosX(i,j,k) exists
    // MeshSize.hpp::UniformMesh::getPosX - Line 66
    // Stencils::TurbulentVTKStencil::writePoints - Line 64
    // Stencils::VelocityStencil/ObstacleStencil::apply - using OBSTACLE
    // BoundaryType.hpp - enum for obstacle

    if ((obstacle & OBSTACLE_SELF) == 0) {    // If this is a fluid cell  
        const RealType posX_lastWall = NormalWallDistance(flowField, OBSTACLE_LEFT, i, j);      // -ve dir
        const RealType posX_nextWall = NormalWallDistance(flowField, OBSTACLE_RIGHT, i, j);     // +ve dir
        const RealType posY_lastWall = NormalWallDistance(flowField, OBSTACLE_BOTTOM, i, j);    // -ve dir
        const RealType posY_nextWall = NormalWallDistance(flowField, OBSTACLE_TOP, i, j);       // +ve dir
        
        // Take the difference and check, otherwise only Pos*_lastWall will be taken 
        const RealType distX_nearest =  (posX - posX_lastWall) < (posX_nextWall - posX) ? 
                                        (posX - posX_lastWall) : (posX_nextWall - posX) ;
        const RealType distY_nearest =  (posY - posY_lastWall) < (posY_nextWall - posY) ? 
                                        (posY - posY_lastWall) : (posY_nextWall - posY) ;

        // Only checking normal wall distance even in BFStep case, 
        // even if the pythagorean distance to the corner can be the smallest one.
        value = distX_nearest < distY_nearest ? distX_nearest : distY_nearest ;

        // if (parameters_.simulation.scenario == "pressure-channel"){
        //     const RealType xLimit = parameters.bfStep.xRatio * parameters.geometry.lengthX;
        //     const RealType yLimit = parameters.bfStep.yRatio * parameters.geometry.lengthY;

        //     if ((xLimit < posX) && (yLimit < posY)){
        //         RealType cornerDist = sqrt((xLimit-posX)*(xLimit-posX) + (yLimit-posY)*(yLimit-posY));
        //         value = value < cornerDist ? value : cornerDist;
        //     }
        // }
    else { value = 0.0; }
    ASSERTION(value >= 0.0);
}

void Stencils::WallhStencil::apply(TurbulentFlowField& flowField, int i, int j, int k)
{
    RealType& value     = flowField.getWallh().getScalar(i, j, k);
    RealType value2D    = apply(flowField, i, j);
    const int obstacle  = flowField.getFlags().getValue(i, j, k);
    
    const RealType posZ     = parameters_.meshsize->getPosZ(i, j, k) 
                            + 0.5*parameters_.meshsize->getDz(i, j, k);

    if ((obstacle & OBSTACLE_SELF) == 0) {    // If this is a fluid cell  
        const RealType posZ_lastWall = NormalWallDistance(flowField, OBSTACLE_FRONT, i, j, k); 
        const RealType posZ_nextWall = NormalWallDistance(flowField, OBSTACLE_BACK, i, j, k);

        const RealType posZ_nearest = posZ_lastWall < posZ_nextWall ? posZ_lastWall : posZ_nextWall;

        value = value2D < posZ_nearest ? value2D : posZ_nearest;
    }
    else { value = 0.0; }
    ASSERTION(value >= 0.0);
}

RealType Stencils::WallhStencil::NormalWallDistance(TurbulentFlowField& flowField, int OBSTACLE_DIR, int i, int j, int k){
    int iter_obstacle   = flowField.getFlags().getValue(i, j, k);
    int iter_i          = i;
    int iter_j          = j;
    int iter_k          = k;
    RealType result     = 0.0;
        
    while((iter_obstacle & OBSTACLE_DIR) != 0){
        if      ((OBSTACLE_DIR &  OBSTACLE_LEFT)==0)    { iter_i -= 1; }
        else if ((OBSTACLE_DIR &  OBSTACLE_RIGHT)==0)   { iter_i += 1; }
        else if ((OBSTACLE_DIR &  OBSTACLE_BOTTOM)==0)  { iter_j -= 1; }
        else if ((OBSTACLE_DIR &  OBSTACLE_TOP)==0)     { iter_j += 1; }
        else if ((OBSTACLE_DIR &  OBSTACLE_FRONT)==0)   { iter_k -= 1; }
        else if ((OBSTACLE_DIR &  OBSTACLE_BACK)==0)    { iter_k += 1; }
        else {}
        iter_obstacle = flowField.getFlags().getValue(iter_i, iter_j, iter_k);
    }

    if ((OBSTACLE_DIR &  OBSTACLE_LEFT)==0){
        result  = parameters_.meshsize->getPosX(iter_i, iter_j, iter_k);
    }
    else if ((OBSTACLE_DIR &  OBSTACLE_RIGHT)==0){
        iter_i += 1;
        result  = parameters_.meshsize->getPosX(iter_i, iter_j, iter_k);
    }
    else if ((OBSTACLE_DIR &  OBSTACLE_BOTTOM)==0){
        result  = parameters_.meshsize->getPosY(iter_i, iter_j, iter_k);
    }
    else if ((OBSTACLE_DIR &  OBSTACLE_TOP)==0){
        iter_j += 1;
        result  = parameters_.meshsize->getPosY(iter_i, iter_j, iter_k);
    }
    else if ((OBSTACLE_DIR &  OBSTACLE_FRONT)==0){
        result  = parameters_.meshsize->getPosZ(iter_i, iter_j, iter_k);
    }
    else if ((OBSTACLE_DIR &  OBSTACLE_BACK)==0){
        iter_k += 1;
        result  = parameters_.meshsize->getPosZ(iter_i, iter_j, iter_k);
    }
    else {}

    return result;
}