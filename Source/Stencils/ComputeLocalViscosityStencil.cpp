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
        RealType turbViscVal = 0.0;
        RealType ChViscVal   = 0.0;
        if (parameters_.simulation.turbModel == "turbSA"){
            // Spalart allmaras turbulence model
            ChViscVal       = flowField.getNewChVis().getScalar(i, j);
            RealType X      = ChViscVal * parameters_.flow.Re;
            RealType fv1    = pow(X, 3.0) / (pow(X, 3.0) + pow(parameters_.turbSA.cv1, 3.0));

            turbViscVal     = fv1 * ChViscVal;
            if (turbViscVal < 0.0) { 
              spdlog::info("At i{}, j{}: turbVisc{}, ChVisc{}", i, j, turbViscVal, ChViscVal);
              throw std::runtime_error("turbVisc at ComputeLocalViscStencil"); 
            }
            flowField.getTurbVisc().getScalar(i, j) = turbViscVal;
            
            // RealType dudy = (flowField.getVelocity().getVector(i, j + 1)[0] - flowField.getVelocity().getVector(i, j - 1)[0])
            //           / (0.5*parameters_.meshsize->getDy(i, j-1) + parameters_.meshsize->getDy(i,j) + 0.5*parameters_.meshsize->getDy(i,j+1));

            // RealType dvdx = (flowField.getVelocity().getVector(i + 1, j)[1] - flowField.getVelocity().getVector(i - 1, j)[1])
            //                 / (0.5*parameters_.meshsize->getDx(i-1, j) + parameters_.meshsize->getDx(i,j) + 0.5*parameters_.meshsize->getDx(i+1,j));

            // RealType W_12 = 0.5 * (dudy - dvdx);

            // RealType chi = flowField.getOldChVis().getScalar(i, j) / (1 / parameters_.flow.Re);

            // RealType c_v2 = 0.7;
            // RealType c_v3 = 0.9;

            // RealType f_t2 = 1.2 * std::exp(-0.5 * chi * chi);

            // RealType temp1 = std::pow(7.1, 3.0); // c_v1 = 7.1
            // RealType temp2 = std::pow(chi, 3.0);

            // RealType f_v1 = temp2 / (temp2 + temp1);

            // RealType f_v2 = 1.0 - (chi / (1 + (chi * f_v1)));

            // RealType temp3 = flowField.getOldChVis().getScalar(i, j)/(parameters_.turbSA.k*parameters_.turbSA.k*(flowField.getWallh().getScalar(i, j)+1e-6)*(flowField.getWallh().getScalar(i, j)+1e-6));

            // RealType S = 2.0 * std::sqrt((W_12 * W_12));

            // RealType S_bar = temp3 * f_v2;

            // RealType S_hat = std::max(0.3 * S, S + S_bar);

            // RealType r = std::min(10.0, temp3 / (fmax(S_hat, 1e-6)));

            // RealType g = r + 0.3 * (std::pow(r, 6.0) - r);

            // RealType cw3_pow6 = 64.0;

            // RealType f_w = g * std::pow(((1 + cw3_pow6) / (std::pow(g, 6.0) + cw3_pow6)), 1.0 / 6.0);

            // // ********************************************************************************
            // // Term 1 (ADVECTION)
            // // *********************************************************************************
            // // Delete either of the below two variables
            // RealType dx0 = parameters_.meshsize->getDx(i, j);

            // RealType term1
            //     = 0.5 * (1 / dx0) * (flowField.getVelocity().getVector(i, j)[0])
            //     * (flowField.getOldChVis().getScalar(i, j) + flowField.getOldChVis().getScalar(i + 1, j));

            // term1 -=  0.5*(1 / dx0) * (flowField.getVelocity().getVector(i-1, j)[0])
            //     * (flowField.getOldChVis().getScalar(i-1, j) + flowField.getOldChVis().getScalar(i, j));

            // term1
            //     += 0.5 * parameters_.solver.gamma * (1 / dx0) * std::fabs(flowField.getVelocity().getVector(i, j)[0])
            //     * (flowField.getOldChVis().getScalar(i, j) - flowField.getOldChVis().getScalar(i + 1, j));

            // term1 -= 0.5 * parameters_.solver.gamma * (1 / dx0) * std::fabs(flowField.getVelocity().getVector(i - 1, j)[0])
            //         * ((
            //             flowField.getOldChVis().getScalar(i - 1, j)
            //             - flowField.getOldChVis().getScalar(i, j)
            //         ));

            // RealType dy0 = parameters_.meshsize->getDy(i, j);

            // term1
            //     += 0.5 * (1 / dy0) * (flowField.getVelocity().getVector(i, j)[1])
            //     * (flowField.getOldChVis().getScalar(i, j) + flowField.getOldChVis().getScalar(i, j + 1));

            // term1 = term1 - 0.5*(1 / dy0) * (flowField.getVelocity().getVector(i, j-1)[1])
            //     * (flowField.getOldChVis().getScalar(i, j-1) + flowField.getOldChVis().getScalar(i, j));

            // term1
            //     += 0.5 * parameters_.solver.gamma * (1 / dy0) * std::fabs(flowField.getVelocity().getVector(i, j)[1])
            //     * (flowField.getOldChVis().getScalar(i, j) - flowField.getOldChVis().getScalar(i, j + 1));

            // term1 -= 0.5 * parameters_.solver.gamma * (1 / dy0) * std::fabs(flowField.getVelocity().getVector(i, j - 1)[1])
            //         * ((
            //             flowField.getOldChVis().getScalar(i, j - 1)
            //             - flowField.getOldChVis().getScalar(i, j)
            //         ));

            // //**********************************************************************
            // // Term 2 (PRODUCTION/SOURCE Term)
            // //**********************************************************************
            // RealType term2 = 0.1355 * (1 - f_t2) * S_hat * flowField.getOldChVis().getScalar(i, j);
            // //**********************************************************************
            // // Term 3 (Wall DESTRUCTION Source Term)
            // //*************************************************************************
            // RealType C_w1 = 3.239067817;

            // RealType temp4
            //     = (flowField.getOldChVis().getScalar(i, j) / (flowField.getWallh().getScalar(i, j) + 1e-6));

            // RealType term3 = (C_w1 * f_w - (0.1355 * f_t2 / (parameters_.turbSA.k * parameters_.turbSA.k)))
            //                 * temp4 * temp4;
            // //************************************************************
            // // Term 4 (DIFFUSION Term)
            // //************************************************************

            // RealType dx          = parameters_.meshsize->getDx(i, j);
            // RealType dx_iplus1j  = parameters_.meshsize->getDx(i + 1, j);
            // RealType dx_iminus1j = parameters_.meshsize->getDx(i - 1, j);

            // RealType dy          = parameters_.meshsize->getDy(i, j);
            // RealType dy_ijplus1  = parameters_.meshsize->getDy(i, j + 1);
            // RealType dy_ijminus1 = parameters_.meshsize->getDy(i, j - 1);

            // RealType nu_ij       = flowField.getOldChVis().getScalar(i, j);
            // RealType nu_iplus1j  = flowField.getOldChVis().getScalar(i + 1, j);
            // RealType nu_iminus1j = flowField.getOldChVis().getScalar(i - 1, j);
            // RealType nu_ijplus1  = flowField.getOldChVis().getScalar(i, j + 1);
            // RealType nu_ijminus1 = flowField.getOldChVis().getScalar(i, j - 1);
            // RealType nu          = 1 / parameters_.flow.Re;

            // RealType viscosity_laplacian_x = (nu + 0.5 * (nu_ij + nu_iplus1j))
            //                                 * ((nu_iplus1j - nu_ij) / (0.5 * (dx_iplus1j + dx)));
            // viscosity_laplacian_x -= (nu + 0.5 * (nu_iminus1j + nu_ij)) * ((nu_ij - nu_iminus1j) / (0.5 * (dx_iminus1j + dx)));

            // viscosity_laplacian_x = viscosity_laplacian_x / dx;

            // RealType viscosity_laplacian_y = (nu + 0.5 * (nu_ij + nu_ijplus1))
            //                                 * ((nu_ijplus1 - nu_ij) / (0.5 * (dy_ijplus1 + dy)));
            // viscosity_laplacian_y -= (nu + 0.5 * (nu_ijminus1 + nu_ij)) * ((nu_ij - nu_ijminus1) / (0.5 * (dy_ijminus1 + dy)));

            // viscosity_laplacian_y = viscosity_laplacian_y / dy;

            // RealType viscosity_laplacian = viscosity_laplacian_x + viscosity_laplacian_y;

            // RealType dnudx = (nu_iplus1j - nu_iminus1j) / ((0.5 * dx_iplus1j) + dx + (0.5 * dx_iminus1j));
            // RealType dnudy = (nu_ijplus1 - nu_ijminus1) / ((0.5 * dy_ijplus1) + dy + (0.5 * dy_ijminus1));

            // RealType advec_term = 0.622 * ((dnudx * dnudx) + (dnudy * dnudy));

            // RealType term4 = (1.5) * (viscosity_laplacian + advec_term);

            // flowField.getNewChVis().getScalar(
            //     i, j
            // ) = flowField.getOldChVis().getScalar(i, j)
            //     + parameters_.timestep.dt * (term2 - term3 + term4 - term1);

            // flowField.getNewChVis().getScalar(i, j) = std::max(
            //     0.0, flowField.getNewChVis().getScalar(i, j)
            // );

            //***********************************************
            // Boundary Conditions for nu_transport:
            //***********************************************

            // Boundary condition for ONLY channel flow
            // if ((parameters_.bfStep.xRatio * parameters_.geometry.sizeX <= 0) && (parameters_.bfStep.yRatio * parameters_.geometry.sizeY <= 0)) {
            //     // at the lower wall of the channel (no slip)
            //     if (((j + parameters_.parallel.firstCorner[1]) == 2) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
            //     flowField.getNewChVis().getScalar(
            //         i, j - 1
            //     ) = -flowField.getNewChVis().getScalar(i, j);
            //     }
            //     // at the top wall of the channel (no slip)
            //     if (((j + parameters_.parallel.firstCorner[1]) == parameters_.geometry.sizeY + 1) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
            //     flowField.getNewChVis().getScalar(
            //         i, j + 1
            //     ) = -flowField.getNewChVis().getScalar(i, j);
            //     }
            //     // at the inlet of the channel
            //     if (((i + parameters_.parallel.firstCorner[0]) == 2) && ((j + parameters_.parallel.firstCorner[1]) > 1)) {
            //     flowField.getNewChVis().getScalar(i - 1, j) = 3.0 / parameters_.flow.Re;
            //     }
            //     // //at the outlet of the channel
            //     if (((i + parameters_.parallel.firstCorner[0]) == parameters_.geometry.sizeX + 1) && ((j + parameters_.parallel.firstCorner[1]) > 1)) {
            //     flowField.getNewChVis().getScalar(
            //         i + 1, j
            //     ) = flowField.getNewChVis().getScalar(i, j);
            //     // std::cout<<"here"<<std::endl;
            //     }

            // }

            // chi  = flowField.getNewChVis().getScalar(i, j) * parameters_.flow.Re;
            // f_v1 = std::pow(chi, 3.0) / (std::pow(chi, 3.0) + std::pow(7.1, 3.0));
            // // calculating eddy viscosity from nu_tilda
            // flowField.getTurbVisc().getScalar(i, j) = f_v1 * flowField.getNewChVis().getScalar(i, j);
            
            // if (i == 2 && j == 2){
                // spdlog::info("term3{}", term3);
                // spdlog::info("CheckVisc: OldCh{} NewCh{} TurbVisc{} S_{} r{}", flowField.getOldChVis().getScalar(i, j), flowField.getNewChVis().getScalar(i, j), flowField.getTurbVisc().getScalar(i, j), S, r);
            // }
            // throw std::runtime_error("Stop at QStencil");
        }
        
        if (parameters_.simulation.turbModel =="turbMixing") {
            // Mixing length turbulence model
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

            turbViscVal = lm * lm * sqrt(2 * SijSij);
            flowField.getTurbVisc().getScalar(i, j) = turbViscVal;
        }

        if (parameters_.simulation.scenario == "channel"){
            // Boundary assignments
            // We know we are inside a fluid cell

            // check if bfStep exists
            const RealType xLimit = parameters_.bfStep.xRatio;
            const RealType yLimit = parameters_.bfStep.yRatio;
            if (xLimit > 0 && yLimit > 0){
                // for turbVisc at Inner walls for BFS case, we "reflect"
                // if (parameters_.simulation.turbModel == "turbMixing"){
                    if ((obstacle & OBSTACLE_LEFT) == OBSTACLE_LEFT) {
                        flowField.getTurbVisc().getScalar(i-1, j) = -1*turbViscVal;
                    }
                    if ((obstacle & OBSTACLE_BOTTOM) == OBSTACLE_BOTTOM) {
                        flowField.getTurbVisc().getScalar(i, j-1) = -1*turbViscVal;
                    }
                // }
                if (parameters_.simulation.turbModel == "turbSA"){
                    if ((obstacle & OBSTACLE_LEFT) == OBSTACLE_LEFT) {
                        flowField.getNewChVis().getScalar(i-1, j) = -1*ChViscVal;
                    }
                    if ((obstacle & OBSTACLE_BOTTOM) == OBSTACLE_BOTTOM) {
                        flowField.getNewChVis().getScalar(i, j-1) = -1*ChViscVal;
                    }
                }
            }

            // In x-dir boundaries, we use same value of turbVisc as the left cell (we "mirror")
            // to maintain the same value when CentralDiff is taken at the wall since the lm is the same
            // if (parameters_.simulation.turbModel == "turbMixing"){
                if (i==2) {
                    flowField.getTurbVisc().getScalar(i-1, j) = 0.22 / (parameters_.flow.Re);
                }
                if (i == parameters_.geometry.sizeX + 1) {
                    flowField.getTurbVisc().getScalar(i+1, j) = turbViscVal;
                }
                // In y-dir boundaries, we use -ve turbVisc of the bottom cell (we "reflect")
                // to maintain a value of zero when CentralDiff is taken at the wall
                if (j==2) {
                    flowField.getTurbVisc().getScalar(i, j-1) =  -1*turbViscVal;
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
            // } 

            if (parameters_.simulation.turbModel == "turbSA"){
                if (i==2) {
                    flowField.getNewChVis().getScalar(i-1, j) = 3.0 / (parameters_.flow.Re);
                }
                if (i == parameters_.geometry.sizeX + 1) {
                    flowField.getNewChVis().getScalar(i+1, j) = ChViscVal;
                }
                if (j==2) {
                    flowField.getNewChVis().getScalar(i, j-1) = -1*ChViscVal;
                }
                if (j == parameters_.geometry.sizeY + 1) {
                    flowField.getNewChVis().getScalar(i, j+1) = -1*ChViscVal;
                }
            }
        }
    }    
    // else{ flowField.getTurbVisc().getScalar(i, j) = 0.0; }
}

void Stencils::ComputeLocalViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
    
    
    const int obstacle  = flowField.getFlags().getValue(i, j, k);
    if ((obstacle & OBSTACLE_SELF) == 0) { // If the cell is fluid
        RealType turbViscVal = 0.0;
        RealType ChViscVal   = 0.0;
        if (parameters_.simulation.turbModel == "turbSA"){
            //Spalart allmaras turbulence model
            RealType X      = flowField.getNewChVis().getScalar(i, j, k) * parameters_.flow.Re;
            RealType fv1    = pow(X, 3) / (pow(X, 3) + pow(parameters_.turbSA.cv1, 3));

            turbViscVal     = fv1 * flowField.getNewChVis().getScalar(i, j, k);
            flowField.getTurbVisc().getScalar(i, j, k) = turbViscVal;
        }

        if (parameters_.simulation.turbModel =="turbMixing") {
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

            loadLocalVelocity3D(flowField, localVelocity_, i, j, k);
            loadLocalMeshsize3D(parameters_, localMeshsize_, i, j, k);

            RealType S11 = dudx(localVelocity_, localMeshsize_);
            RealType S22 = dvdy(localVelocity_, localMeshsize_);
            RealType S33 = dwdz(localVelocity_, localMeshsize_);
            RealType S12 = 0.5 * (dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_));
            RealType S13 = 0.5 * (dudz(localVelocity_, localMeshsize_) + dwdx(localVelocity_, localMeshsize_));
            RealType S23 = 0.5 * (dvdz(localVelocity_, localMeshsize_) + dwdy(localVelocity_, localMeshsize_));

            RealType SijSij = S11 * S11 + S22 * S22 + S33 * S33 + 2 * S12 * S12 + 2 * S13 * S13 + 2 * S23 * S23;

            turbViscVal = lm * lm * sqrt(2 * SijSij);
            flowField.getTurbVisc().getScalar(i, j, k) = turbViscVal;
        }

        if (parameters_.simulation.scenario == "channel"){
            // Boundary assignments
            // We know we are inside a fluid cell

            // check if bfStep exists
            const RealType xLimit = parameters_.bfStep.xRatio;
            const RealType yLimit = parameters_.bfStep.yRatio;
            if (xLimit > 0 && yLimit > 0){
                // for turbVisc at Inner walls for BFS case, we "reflect"
                // if (parameters_.simulation.turbModel == "turbMixing"){
                    if ((obstacle & OBSTACLE_LEFT) == OBSTACLE_LEFT) {
                        flowField.getTurbVisc().getScalar(i-1, j, k) = -1*turbViscVal;
                    }
                    if ((obstacle & OBSTACLE_BOTTOM) == OBSTACLE_BOTTOM) {
                        flowField.getTurbVisc().getScalar(i, j-1, k) = -1*turbViscVal;
                    }
                // }
                if (parameters_.simulation.turbModel == "turbSA"){
                    if ((obstacle & OBSTACLE_LEFT) == OBSTACLE_LEFT) {
                        flowField.getNewChVis().getScalar(i-1, j, k) = -1*ChViscVal;
                    }
                    if ((obstacle & OBSTACLE_BOTTOM) == OBSTACLE_BOTTOM) {
                        flowField.getNewChVis().getScalar(i, j-1, k) = -1*ChViscVal;
                    }
                }
            }

            // In x-dir boundaries, we use same value of turbVisc as the left cell (we "mirror")
            // to maintain the same value when CentralDiff is taken at the wall since the lm is the same
            // if (parameters_.simulation.turbModel == "turbMixing"){
                if (i==2) {
                    flowField.getTurbVisc().getScalar(i-1, j, k) = 0.22 / (parameters_.flow.Re);
                }
                if (i == parameters_.geometry.sizeX + 1) {
                    flowField.getTurbVisc().getScalar(i+1, j, k) = turbViscVal;
                }
                // In y-dir and z-dir boundaries, we use -ve turbVisc of the bottom cell (we "reflect")
                // to maintain a value of zero when CentralDiff is taken at the wall
                if (j==2) {
                    flowField.getTurbVisc().getScalar(i, j-1, k) =  -1*turbViscVal;
                }
                if (j == parameters_.geometry.sizeY + 1) {
                    flowField.getTurbVisc().getScalar(i, j+1, k) = -1*turbViscVal;
                }
                if (k==2) {
                    flowField.getTurbVisc().getScalar(i, j, k-1) =  -1*turbViscVal;
                }
                if (k == parameters_.geometry.sizeZ + 1) {
                    flowField.getTurbVisc().getScalar(i, j, k+1) = -1*turbViscVal;
                }
            // } 

            if (parameters_.simulation.turbModel == "turbSA"){
                if (i==2) {
                    flowField.getNewChVis().getScalar(i-1, j, k) = 3.0 / (parameters_.flow.Re);
                }
                if (i == parameters_.geometry.sizeX + 1) {
                    flowField.getNewChVis().getScalar(i+1, j, k) = ChViscVal;
                }
                if (j==2) {
                    flowField.getNewChVis().getScalar(i, j-1, k) = -1*ChViscVal;
                }
                if (j == parameters_.geometry.sizeY + 1) {
                    flowField.getNewChVis().getScalar(i, j+1, k) = -1*ChViscVal;
                }
                if (k==2) {
                    flowField.getNewChVis().getScalar(i, j, k-1) = -1*ChViscVal;
                }
                if (k == parameters_.geometry.sizeZ + 1) {
                    flowField.getNewChVis().getScalar(i, j, k+1) = -1*ChViscVal;
                }
            }
        }
    }
}