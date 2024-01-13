#include "StdAfx.hpp"

#include "WallhStencil.hpp"

Stencils::WallhStencil::WallhStencil(const Parameters& parameters):
FieldStencil<TurbulentFlowField>(parameters){}

void Stencils::WallhStencil::apply(TurbulentFlowField& flowField, int i, int j)
{
    RealType& value     = flowField.getWallh().getScalar(i, j);
    const int obstacle  = flowField.getFlags().getValue(i, j);

    if (((obstacle & OBSTACLE_SELF) == 0)) { // If this is a fluid cell
        value = get2Dwallh(i, j);
    }
    else { value = 0.0; }
    ASSERTION(value >= 0.0);
}

void Stencils::WallhStencil::apply(TurbulentFlowField& flowField, int i, int j, int k)
{
    RealType& value     = flowField.getWallh().getScalar(i, j, k);
    const int obstacle  = flowField.getFlags().getValue(i, j, k);
    
    if (((obstacle & OBSTACLE_SELF) == 0)) {  // If this is a fluid cell
        const RealType posZ     = parameters_.meshsize->getPosZ(i, j, k) + 0.5*parameters_.meshsize->getDz(i, j, k);
        const RealType wallZ    = parameters_.geometry.lengthZ;
        // get nearest wall dist in 2D
        const RealType value2D = get2Dwallh(i, j);
        // nearest wall dist in z-dir
        RealType distZ_nearest  = (posZ) < (wallZ - posZ) ? (posZ) : (wallZ - posZ) ;
        // compare all dir for final value
        value = value2D < distZ_nearest ? value2D : distZ_nearest;
        // if(i==2 && j==2 && k==2){
        //     spdlog::info("value2D: {}, distZ: {}, value: {}, bin: {}", value2D, distZ_nearest, value, flowField.getWallh().getScalar(i,j,0));
        // }
        // if(i==4 && j==4 && k==4){
        //     spdlog::info("value2D: {}, distZ: {}, value: {}, bin: {}", value2D, distZ_nearest, value, flowField.getWallh().getScalar(i,j,0));
        //     throw std::runtime_error("to stop");
        // }
    } // quick fix
    else { value = 0.0; }
    ASSERTION(value >= 0.0);
}

RealType Stencils::WallhStencil::get2Dwallh(int i, int j){
    RealType value = 0.0;
    const RealType posX     = parameters_.meshsize->getPosX(i, j) + 0.5*parameters_.meshsize->getDx(i, j); 
    const RealType wallX    = parameters_.geometry.lengthX;     
    const RealType posY     = parameters_.meshsize->getPosY(i, j) + 0.5*parameters_.meshsize->getDy(i, j); 
    const RealType wallY    = parameters_.geometry.lengthY;
    // neareast wall dist in y-dir
    RealType nearWallY = (posY) < (wallY - posY) ? (posY) : (wallY - posY) ; 
    // check scenarios
    if (parameters_.simulation.scenario == "cavity"){
        // neareast wall dist in x-dir
        RealType nearWallX = (posX) < (wallX - posX) ? (posX) : (wallX - posX) ; 
        value = nearWallX < nearWallY ? nearWallX : nearWallY ;
    }
    if (parameters_.simulation.scenario == "channel"){ 
        value = nearWallY ; // only y-dir are walls
        // check if step exists
        const RealType xLimit = parameters_.bfStep.xRatio * parameters_.geometry.lengthX;
        const RealType yLimit = parameters_.bfStep.yRatio * parameters_.geometry.lengthY;
        if (xLimit > 0 && yLimit > 0){
            // spdlog::info("Stick with simple Channel2D case");
            // throw std::runtime_error("BF step detected");
            // Only checking normal wall distance even in BFStep case, 
            // even if the pythagorean distance to the corner can be the smallest one.
            if (posX < xLimit){ // Q2
                value =  (posY - yLimit) < (wallY - posY) ? (posY - yLimit) : (wallY - posY) ; // new wall dist in y-dir
            }
            if (posY < yLimit){ // Q4
                value = value < (posX - xLimit) ? value : (posX - xLimit) ; // walls in -ve x-dir too
            }
        }
    }
    // if (value < 0) { value = 0; } // quick fix
    return value;
}