#include "StdAfx.hpp"

#include "ComputeLocalViscosityStencil.hpp"

#include "Definitions.hpp"
#include "TurbulentStencilFunctions.hpp"

Stencils::ComputeLocalViscosityStencil::ComputeLocalViscosityStencil(const Parameters& parameters):
FieldStencil<TurbulentFlowField>(parameters){}

void Stencils::ComputeLocalViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j) {
    // To account for turbVisc value at the booundaries we have some if statements
    // NOTE: this whole implementation is made with parameters_.simulation.scenario ="channel" in mind
   
    const int obstacle  = flowField.getFlags().getValue(i, j);
    if ((obstacle & OBSTACLE_SELF) == 0) { // If the cell is fluid
        
        RealType delta = 0.0; 
        if (parameters_.turbMix.delta == 0) { //Case for no boundary layer
            delta = 0.0;
        } 
        else if (parameters_.turbMix.delta == 1) { //Case for laminar flat plate
            auto x = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i, j);
            delta  = 4.91 * x / sqrt(x * parameters_.flow.Re);
        } 
        else if (parameters_.turbMix.delta == 2) { //Case for turbulent flat plate
            auto x = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i, j);
            delta  = 0.382 * x / pow(x * parameters_.flow.Re, 0.2);
        }

        //Calculate mixing length by Prandtl's model
        RealType lm = std::min(parameters_.turbMix.k * flowField.getWallh().getScalar(i, j), 0.09 * delta);
        
        loadLocalVelocity2D(flowField, localVelocity_, i, j);
        loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);

        RealType S11 = dudx(localVelocity_, localMeshsize_);
        RealType S22 = dvdy(localVelocity_, localMeshsize_);
        RealType S12 = 0.5 * (dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_));

        RealType SijSij = S11 * S11 + S22 * S22 + 2 * S12 * S12;

        RealType turbViscVal = lm * lm * sqrt(2 * SijSij);
        flowField.getTurbVisc().getScalar(i, j) = turbViscVal;

        // check if step exists
        const RealType xLimit = parameters_.bfStep.xRatio * parameters_.geometry.lengthX;
        const RealType yLimit = parameters_.bfStep.yRatio * parameters_.geometry.lengthY;
        if (xLimit > 0 && yLimit > 0){
            // for turbVisc at Inner walls for BFS case, we "reflect"
            if ((obstacle & OBSTACLE_LEFT) == 2) {
                flowField.getTurbVisc().getScalar(i-1, j) = -1*turbViscVal;
            }
            if ((obstacle & OBSTACLE_BOTTOM) == 8) {
                flowField.getTurbVisc().getScalar(i, j-1) = -1*turbViscVal;
            }
        }

        // Boundary assignments
        // We know we are inside a fluid cell

        // In x-dir boundaries, we use same value of turbVisc as the left cell (we "mirror")
        // to maintain the same value when CentralDiff is taken at the wall since the lm is the same
        if (i==2) {
            flowField.getTurbVisc().getScalar(i-1, j) = -1*turbViscVal;
        }
        if (i == parameters_.geometry.sizeX + 1) {
            flowField.getTurbVisc().getScalar(i+1, j) = -1*turbViscVal;
        }
        // In y-dir boundaries, we use -ve turbVisc of the bottom cell (we "reflect")
        // to maintain a value of zero when CentralDiff is taken at the wall
        if (j==2) {
            flowField.getTurbVisc().getScalar(i, j-1) = -1*turbViscVal;
        }
        if (j == parameters_.geometry.sizeY + 1) {
            flowField.getTurbVisc().getScalar(i, j+1) = -1*turbViscVal;
        }
        // if(i==parameters_.geometry.sizeX + 1 && j>=2){
        //     spdlog::info("at i:{}, j:{} lm ={}",i,j,lm);
        //     if (j ==parameters_.geometry.sizeY + 1){ 
        //         throw std::runtime_error("one iter completed");
        //     }
        // } 
    }
    else{ flowField.getTurbVisc().getScalar(i, j) = 0.0; }
}

void Stencils::ComputeLocalViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
    
    const int obstacle  = flowField.getFlags().getValue(i, j, k);
    if ((obstacle & OBSTACLE_SELF) == 0) { // If the cell is fluid

        RealType delta = 0.0;
        if (parameters_.turbMix.delta == 0) {
            delta = 0.0;
        } 
        else if (parameters_.turbMix.delta == 1) {
            auto x = parameters_.meshsize->getPosX(i, j, k) + 0.5 * parameters_.meshsize->getDx(i, j, k);
            delta  = 4.91 * x / sqrt(x * parameters_.flow.Re);
        } 
        else if (parameters_.turbMix.delta == 2){
            auto x = parameters_.meshsize->getPosX(i, j, k) + 0.5 * parameters_.meshsize->getDx(i, j, k);
            // auto y = parameters_.meshsize->getPosY(i, j) + 0.5 * parameters_.meshsize->getDy(i, j);
            delta  = 0.382 * x / pow(x * parameters_.flow.Re, 0.2);
            // spdlog::info("delta: {} for x,y: {},{} at {},{}",delta, x, y, i, j);
        }

        RealType lm = std::min(parameters_.turbMix.k * flowField.getWallh().getScalar(i, j, k), 0.09 * delta);
        if(lm<=0){
           spdlog::info("lm: {}",lm);
            throw std::runtime_error("to stop");
        }

        loadLocalVelocity3D(flowField, localVelocity_, i, j, k);
        loadLocalMeshsize3D(parameters_, localMeshsize_, i, j, k);

        RealType S11 = dudx(localVelocity_, localMeshsize_);
        RealType S22 = dvdy(localVelocity_, localMeshsize_);
        RealType S33 = dwdz(localVelocity_, localMeshsize_);
        RealType S12 = 0.5 * (dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_));
        RealType S13 = 0.5 * (dudz(localVelocity_, localMeshsize_) + dwdx(localVelocity_, localMeshsize_));
        RealType S23 = 0.5 * (dvdz(localVelocity_, localMeshsize_) + dwdy(localVelocity_, localMeshsize_));

        RealType SijSij = S11 * S11 + S22 * S22 + S33 * S33 + 2 * S12 * S12 + 2 * S13 * S13 + 2 * S23 * S23;

        RealType turbViscVal = lm * lm * sqrt(2 * SijSij);
        flowField.getTurbVisc().getScalar(i, j, k) = turbViscVal;

        // check if step exists
        const RealType xLimit = parameters_.bfStep.xRatio * parameters_.geometry.lengthX;
        const RealType yLimit = parameters_.bfStep.yRatio * parameters_.geometry.lengthY;
        if (xLimit > 0 && yLimit > 0){
            // for turbVisc at Inner walls for BFS case, we "reflect"
            if ((obstacle & OBSTACLE_LEFT) == 2) {
                flowField.getTurbVisc().getScalar(i-1, j, k) = -1*turbViscVal;
            }
            if ((obstacle & OBSTACLE_BOTTOM) == 8) {
                flowField.getTurbVisc().getScalar(i, j-1, k) = -1*turbViscVal;
            }
        }

        // Boundary assignments
        // In x-dir boundaries, we "mirror"
        if (i==2) {
            flowField.getTurbVisc().getScalar(i-1, j, k) = -1*turbViscVal;
        }
        if (i == parameters_.geometry.sizeX + 1) {
            flowField.getTurbVisc().getScalar(i+1, j, k) = -1*turbViscVal;
        }
        // In y-dir boundaries, we "reflect"
        if (j==2) {
            flowField.getTurbVisc().getScalar(i, j-1, k) = -1*turbViscVal;
        }
        if (j == parameters_.geometry.sizeY + 1) {
            flowField.getTurbVisc().getScalar(i, j+1, k) = -1*turbViscVal;
        }
        // In z-dir boundaries, we "reflect" due to walls, same reason as y-dir
        if (k==2) {
            flowField.getTurbVisc().getScalar(i, j, k-1) = -1*turbViscVal;
        }
        if (k == parameters_.geometry.sizeZ + 1) {
            flowField.getTurbVisc().getScalar(i, j, k+1) = -1*turbViscVal;
        }
    }
    else { flowField.getTurbVisc().getScalar(i, j, k) = 0.0; }
}