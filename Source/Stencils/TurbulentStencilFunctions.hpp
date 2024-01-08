#pragma once

#include "Definitions.hpp"
#include "Parameters.hpp"
#include "StencilFunctions.hpp" 
// only for access to common functions. Beware of same functionalities

namespace Stencils {

/*Functions below are used for turbulent FGH calculation*/

//Load local viscosity 2D
  inline void loadLocalViscosity2D(TurbulentFlowField& flowField, RealType* const localViscosity, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        localViscosity[39 + 9 * row + 3 * column] = flowField.getTurbVisc().getScalar(i + column, j + row);
        //Viscosity in each cell is just scalar, we only have one component, so just assigned it to x.
      }
    }
  }

//Load local viscosity 3D
  inline void loadLocalViscosity3D(TurbulentFlowField& flowField, RealType* const localViscosity, int i, int j, int k) {
    for (int layer = -1; layer <= 1; layer++) {
      for (int row = -1; row <= 1; row++) {
        for (int column = -1; column <= 1; column++) {
          localViscosity[39 + 27 * layer + 9 * row + 3 * column] = flowField.getTurbVisc().getScalar(
            i + column, j + row, k + layer); 
          //Viscosity in each cell is just scalar, we only have one component, so just  assigned it to x.
        }
      }
    }
  }

//For shear rate tensor calculation (upwind differencing scheme)
  inline RealType dudy(const RealType* const lv, const RealType* const lm) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);

    return (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, -1, 0, 0)]) / dy0;
  }

  inline RealType dudz(const RealType* const lv, const RealType* const lm) {
    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_M1 = lm[mapd(0, 0, -1, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);

    return (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, 0, -1, 0)]) / dz0;
  }

  inline RealType dvdx(const RealType* const lv, const RealType* const lm) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);

    return (lv[mapd(0, 0, 0, 1)] - lv[mapd(-1, 0, 0, 1)]) / dx0;
  }

  inline RealType dvdz(const RealType* const lv, const RealType* const lm) {
    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_M1 = lm[mapd(0, 0, -1, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);

    return (lv[mapd(0, 0, 0, 1)] - lv[mapd(0, 0, -1, 1)]) / dz0;
  }

  inline RealType dwdx(const RealType* const lv, const RealType* const lm) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);

    return (lv[mapd(0, 0, 0, 2)] - lv[mapd(-1, 0, 0, 2)]) / dx0;
  }

  inline RealType dwdy(const RealType* const lv, const RealType* const lm) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);

    return (lv[mapd(0, 0, 0, 2)] - lv[mapd(0, -1, 0, 2)]) / dy0;
  }

  // Second derivative of turbulent u-component
  inline RealType d2udx2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    // Evaluate the second derivative at the location of the u-component of the velocity field;
    // we therefore use the two neighbouring u-components and assume arbitrary mesh sizes in both
    // directions -> the formula arises from a straight-forward taylor expansion
    // -> for equal meshsizes, we obtain the usual [1 -2 1]-like stencil.
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);//Average mesh size for streched mesh 

    RealType dudxR = (lv[mapd(1, 0, 0, 0)] - lv[mapd(0, 0, 0, 0)]) / dx_P1;
    RealType dudxL = (lv[mapd(0, 0, 0, 0)] - lv[mapd(-1, 0, 0, 0)]) / dx_0;

    RealType ViscR = lvis[mapd(1, 0, 0, 0)] + 1 / parameters.flow.Re;
    RealType ViscL = lvis[mapd(0, 0, 0, 0)] + 1 / parameters.flow.Re;
  //Reynolds number here stand for reciprocal of viscosity cause velocity and characteristic length are been normalised.

    return (ViscR * dudxR - ViscL * dudxL) / dx1;
  }

  inline RealType d2udy2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);//Average mesh size for streched mesh 

    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);//Average mesh size for streched mesh 
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);//Average mesh size for streched mesh 
    
    RealType dudyT = (lv[mapd(0, 1, 0, 0)] - lv[mapd(0, 0, 0, 0)]) / dy1;
    RealType dudyB = (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, -1, 0, 0)]) / dy0;

    RealType ViscT1 = lvis[mapd(0, 1, 0, 0)];
    RealType ViscT2 = lvis[mapd(1, 1, 0, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(1, 0, 0, 0)];
    RealType ViscB1 = lvis[mapd(0, -1, 0, 0)];
    RealType ViscB2 = lvis[mapd(1, -1, 0, 0)];

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dx_P1 + ViscT2 * dx_0) / dx1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dy_P1
                     ) / dy1
                     + 1 / parameters.flow.Re; 

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dx_P1 + ViscB2 * dx_0) / dx1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dy_M1
                     ) / dy0
                     + 1 / parameters.flow.Re;

    return (ViscT * dudyT - ViscB * dudyB) / dy_0;
  }

  inline RealType d2udz2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType dz_M1 = lm[mapd(0, 0, -1, 2)];
    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);
    
    RealType dudzT = (lv[mapd(0, 0, 1, 0)] - lv[mapd(0, 0, 0, 0)]) / dz1;
    RealType dudzB = (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, 0, -1, 0)]) / dz0;

    RealType ViscT1 = lvis[mapd(0, 0, 1, 0)];
    RealType ViscT2 = lvis[mapd(1, 0, 1, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(1, 0, 0, 0)];
    RealType ViscB1 = lvis[mapd(0, 0, -1, 0)];
    RealType ViscB2 = lvis[mapd(1, 0, -1, 0)];

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dx_P1 + ViscT2 * dx_0) / dx1) * dz_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dz_P1
                     ) / dz1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dx_P1 + ViscB2 * dx_0) / dx1) * dz_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dz_M1
                     ) / dz0
                     + 1 / parameters.flow.Re;

    return (ViscT * dudzT - ViscB * dudzB) / dz_0;
  }

inline RealType d2udydx(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType ViscR1 = lvis[mapd(1, 0, 0, 0)];
    const RealType ViscR2 = lvis[mapd(1, 1, 0, 0)];
    const RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    const RealType ViscM2 = lvis[mapd(0, 1, 0, 0)];
    const RealType ViscL1 = lvis[mapd(-1, 0, 0, 0)];
    const RealType ViscL2 = lvis[mapd(-1, 1, 0, 0)];

    RealType dudyR = (lv[mapd(0, 1, 0, 0)] - lv[mapd(0, 0, 0, 0)]) / dy1;
    RealType dudyL = (lv[mapd(-1, 1, 0, 0)] - lv[mapd(-1, 0, 0, 0)]) / dy1;

    RealType ViscR = (0.5 * (0.5 * (ViscR1 * dx_0 + ViscM1 * dx_P1) / dx1) * dy_P1
                      + 0.5 * (0.5 * (ViscR2 * dx_0 + ViscM2 * dx_P1) / dx1) * dy_0
                     ) / dy1
                     + 1 / parameters.flow.Re;

    RealType ViscL = (0.5 * (0.5 * (ViscL1 * dx_0 + ViscM1 * dx_M1) / dx0) * dy_P1
                      + 0.5 * (0.5 * (ViscL2 * dx_0 + ViscM2 * dx_M1) / dx0) * dy_0
                     ) / dy1
                     + 1 / parameters.flow.Re;

    return ((ViscR * dudyR) - (ViscL * dudyL)) / dx_0;
  }

  inline RealType d2udzdx(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);

    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    RealType ViscT1 = lvis[mapd(1, 0, 0, 0)];
    RealType ViscT2 = lvis[mapd(1, 0, 1, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(0, 0, 1, 0)];
    RealType ViscB1 = lvis[mapd(-1, 0, 0, 0)];
    RealType ViscB2 = lvis[mapd(-1, 0, 1, 0)];

    RealType dudzT = (lv[mapd(0, 0, 1, 0)] - lv[mapd(0, 0, 0, 0)]) / dz1;
    RealType dudzB = (lv[mapd(-1, 0, 1, 0)] - lv[mapd(-1, 0, 0, 0)]) / dz1;

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dx_0 + ViscM1 * dx_P1) / dx1) * dz_P1
                      + 0.5 * (0.5 * (ViscT2 * dx_0 + ViscM2 * dx_P1) / dx1) * dz_0
                     ) / dz1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dx_0 + ViscM1 * dx_M1) / dx0) * dz_P1
                      + 0.5 * (0.5 * (ViscB2 * dx_0 + ViscM2 * dx_M1) / dx0) * dz_0
                     ) / dz1
                     + 1 / parameters.flow.Re;

    return (ViscT * dudzT - ViscB * dudzB) / dx_0;
  }

  //Second derivative for v component
inline RealType d2vdx2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType ViscR1 = lvis[mapd(1, 0, 0, 0)];
    const RealType ViscR2 = lvis[mapd(1, 1, 0, 0)];
    const RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    const RealType ViscM2 = lvis[mapd(0, 1, 0, 0)];
    const RealType ViscL1 = lvis[mapd(-1, 0, 0, 0)];
    const RealType ViscL2 = lvis[mapd(-1, 1, 0, 0)];

    RealType dvdxR = (lv[mapd(1, 0, 0, 1)] - lv[mapd(0, 0, 0, 1)]) / dx1;
    RealType dvdxL = (lv[mapd(0, 0, 0, 1)] - lv[mapd(-1, 0, 0, 1)]) / dx0;

    RealType ViscR = (0.5 * (0.5 * (ViscR1 * dx_0 + ViscM1 * dx_P1) / dx1) * dy_P1
                      + 0.5 * (0.5 * (ViscR2 * dx_0 + ViscM2 * dx_P1) / dx1) * dy_0
                     ) / dy1
                     + 1 / parameters.flow.Re;

    RealType ViscL = (0.5 * (0.5 * (ViscL1 * dx_0 + ViscM1 * dx_M1) / dx0) * dy_P1
                      + 0.5 * (0.5 * (ViscL2 * dx_0 + ViscM2 * dx_M1) / dx0) * dy_0
                     ) / dy1
                     + 1 / parameters.flow.Re;
                     
    return (ViscR * dvdxR - ViscL * dvdxL) / dx_0;
  }

inline RealType d2vdy2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    RealType dvdyT = (lv[mapd(0, 1, 0, 1)] - lv[mapd(0, 0, 0, 1)]) / dy_P1;
    RealType dvdyB = (lv[mapd(0, 0, 0, 1)] - lv[mapd(0, -1, 0, 1)]) / dy_0;

    RealType ViscT = lvis[mapd(0, 1, 0, 0)] + 1 / parameters.flow.Re;
    RealType ViscB = lvis[mapd(0, 0, 0, 0)] + 1 / parameters.flow.Re;

    return (ViscT * dvdyT - ViscB * dvdyB) / dy1;
  }

  inline RealType d2vdz2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(1, 0, 0, 2)];
    const RealType dz_M1 = lm[mapd(-1, 0, 0, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);

    const RealType ViscR1 = lvis[mapd(0, 0, 1, 0)];
    const RealType ViscR2 = lvis[mapd(0, 1, 1, 0)];
    const RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    const RealType ViscM2 = lvis[mapd(0, 1, 0, 0)];
    const RealType ViscL1 = lvis[mapd(0, 0, -1, 0)];
    const RealType ViscL2 = lvis[mapd(0, 1, -1, 0)];

    RealType dvdzR = (lv[mapd(1, 0, 0, 1)] - lv[mapd(0, 0, 0, 1)]) / dz1;
    RealType dvdzL = (lv[mapd(0, 0, 0, 1)] - lv[mapd(-1, 0, 0, 1)]) / dz0;

    RealType ViscR = (0.5 * (0.5 * (ViscR1 * dz_0 + ViscM1 * dz_P1) / dz1) * dy_P1
                      + 0.5 * (0.5 * (ViscR2 * dz_0 + ViscM2 * dz_P1) / dz1) * dy_0
                     ) / dy1
                     + 1 / parameters.flow.Re;

    RealType ViscL = (0.5 * (0.5 * (ViscL1 * dz_0 + ViscM1 * dz_M1) / dz0) * dy_P1
                      + 0.5 * (0.5 * (ViscL2 * dz_0 + ViscM2 * dz_M1) / dz0) * dy_0
                     ) / dy1
                     + 1 / parameters.flow.Re;

    return (ViscR * dvdzR - ViscL * dvdzL) / dz_0;
  }

  inline RealType d2vdxdy(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    RealType dvdxT = lv[mapd(1, 0, 0, 1)] - lv[mapd(0, 0, 0, 1)] / dx1;
    RealType dvdxB = lv[mapd(1, -1, 0, 1)] - lv[mapd(0, -1, 0, 1)] / dx1;

    RealType ViscT1 = lvis[mapd(0, 1, 0, 0)];
    RealType ViscT2 = lvis[mapd(1, 1, 0, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(1, 0, 0, 0)];
    RealType ViscB1 = lvis[mapd(0, -1, 0, 0)];
    RealType ViscB2 = lvis[mapd(1, -1, 0, 0)];

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dx_P1 + ViscT2 * dx_0) / dx1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dy_P1
                     ) / dy1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dx_P1 + ViscB2 * dx_0) / dx1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dy_M1
                     ) / dy0
                     + 1 / parameters.flow.Re;

    return (ViscT * dvdxT - ViscB * dvdxB) / dy_0;
  }

  inline RealType d2vdzdy(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(1, 0, 0, 2)];
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);

    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    RealType dvdzT = lv[mapd(0, 0, 1, 1)] - lv[mapd(0, 0, 0, 1)] / dz1;
    RealType dvdzB = lv[mapd(0, -1, 1, 1)] - lv[mapd(0, -1, 0, 1)] / dz1;

    RealType ViscT1 = lvis[mapd(0, 1, 0, 0)];
    RealType ViscT2 = lvis[mapd(0, 1, 1, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(0, 0, 1, 0)];
    RealType ViscB1 = lvis[mapd(0, -1, 0, 0)];
    RealType ViscB2 = lvis[mapd(0, -1, 1, 0)];

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dz_P1 + ViscT2 * dz_0) / dz1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dz_P1 + ViscM2 * dz_0) / dz1) * dy_P1
                     ) / dy1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dz_P1 + ViscB2 * dz_0) / dz1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dz_P1 + ViscM2 * dz_0) / dz1) * dy_M1
                     ) / dy0
                     + 1 / parameters.flow.Re;

    return (ViscT * dvdzT - ViscB * dvdzB) / dy_0;
  }

//Second derivative for turbulent w component 
inline RealType d2wdx2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);

    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType ViscT1 = lvis[mapd(1, 0, 0, 0)];
    const RealType ViscT2 = lvis[mapd(1, 0, 1, 0)];
    const RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    const RealType ViscM2 = lvis[mapd(0, 0, 1, 0)];
    const RealType ViscB1 = lvis[mapd(-1, 0, 0, 0)];
    const RealType ViscB2 = lvis[mapd(-1, 0, 1, 0)];

    RealType dwdxT = (lv[mapd(1, 0, 0, 2)] - lv[mapd(0, 0, 0, 2)]) / dx1;
    RealType dwdxB = (lv[mapd(0, 0, 0, 2)] - lv[mapd(-1, 0, 0, 2)]) / dx0;

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dx_0 + ViscM1 * dx_P1) / dx1) * dz_P1
                      + 0.5 * (0.5 * (ViscT2 * dx_0 + ViscM2 * dx_P1) / dx1) * dz_0
                     ) / dz1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dx_0 + ViscM1 * dx_M1) / dx0) * dz_P1
                      + 0.5 * (0.5 * (ViscB2 * dx_0 + ViscM2 * dx_M1) / dx0) * dz_0
                     ) / dz1
                     + 1 / parameters.flow.Re;

    return (ViscT * dwdxT - ViscB * dwdxB) / dx_0;
  }

  inline RealType d2wdy2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);

    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType ViscT1 = lvis[mapd(0, 1, 0, 0)];
    const RealType ViscT2 = lvis[mapd(0, 1, 1, 0)];
    const RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    const RealType ViscM2 = lvis[mapd(0, 0, 1, 0)];
    const RealType ViscB1 = lvis[mapd(0, -1, 0, 0)];
    const RealType ViscB2 = lvis[mapd(0, -1, 1, 0)];

    RealType dwdyT = (lv[mapd(0, 1, 0, 2)] - lv[mapd(0, 0, 0, 2)]) / dy1;
    RealType dwdyB = (lv[mapd(0, 0, 0, 2)] - lv[mapd(0, -1, 0, 2)]) / dy0;

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dy_0 + ViscM1 * dy_P1) / dy1) * dz_P1
                      + 0.5 * (0.5 * (ViscT2 * dy_0 + ViscM2 * dy_P1) / dy1) * dz_0
                     ) / dz1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dy_0 + ViscM1 * dy_M1) / dy0) * dz_P1
                      + 0.5 * (0.5 * (ViscB2 * dy_0 + ViscM2 * dy_M1) / dy0) * dz_0
                     ) / dz1
                     + 1 / parameters.flow.Re;

    return (ViscT * dwdyT - ViscB * dwdyB) / dy_0;
    }

  inline RealType d2wdz2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);

    RealType dwdzT = (lv[mapd(0, 0, 1, 2)] - lv[mapd(0, 0, 0, 2)]) / dz_P1;
    RealType dwdzB = (lv[mapd(0, 0, 0, 2)] - lv[mapd(0, 0, -1, 2)]) / dz_0;

    RealType ViscT = lvis[mapd(0, 0, 1, 0)] + 1 / parameters.flow.Re;
    RealType ViscB = lvis[mapd(0, 0, 0, 0)] + 1 / parameters.flow.Re;

    return (ViscT * dwdzT - ViscB * dwdzB) / dz1;
  }

inline RealType d2wdxdz(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz_M1 = lm[mapd(0, 0, -1, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);

    RealType dwdxT = lv[mapd(1, 0, 0, 2)] - lv[mapd(0, 0, 0, 2)] / dx1;
    RealType dwdxB = lv[mapd(1, 0, -1, 2)] - lv[mapd(0, 0, -1, 2)] / dx1;

    RealType ViscT1 = lvis[mapd(0, 0, 1, 0)];
    RealType ViscT2 = lvis[mapd(1, 0, 1, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(1, 0, 0, 0)];
    RealType ViscB1 = lvis[mapd(0, 0, -1, 0)];
    RealType ViscB2 = lvis[mapd(1, 0, -1, 0)];

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dx_P1 + ViscT2 * dx_0) / dx1) * dz_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dz_P1
                     ) / dz1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dx_P1 + ViscB2 * dx_0) / dx1) * dz_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dz_M1
                     ) / dz0
                     + 1 / parameters.flow.Re;

    return (ViscT * dwdxT - ViscB * dwdxB) / dz_0;
  }

inline RealType d2wdydz(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(1, 0, 0, 1)];
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz_M1 = lm[mapd(0, 0, -1, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);

    RealType dwdyT = lv[mapd(0, 1, 0, 2)] - lv[mapd(0, 0, 0, 2)] / dy1;
    RealType dwdyB = lv[mapd(0, 1, -1, 2)] - lv[mapd(0, 0, -1, 2)] / dy1;

    RealType ViscT1 = lvis[mapd(0, 0, 1, 0)];
    RealType ViscT2 = lvis[mapd(0, 1, 1, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(0, 1, 0, 0)];
    RealType ViscB1 = lvis[mapd(0, 0, -1, 0)];
    RealType ViscB2 = lvis[mapd(0, 1, -1, 0)];

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dy_P1 + ViscT2 * dy_0) / dy1) * dz_0
                      + 0.5 * (0.5 * (ViscM1 * dy_P1 + ViscM2 * dy_0) / dy1) * dz_P1
                     ) / dz1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dy_P1 + ViscB2 * dy_0) / dy1) * dz_0
                      + 0.5 * (0.5 * (ViscM1 * dy_P1 + ViscM2 * dy_0) / dy1) * dz_M1
                     ) / dz0
                     + 1 / parameters.flow.Re;

    return (ViscT * dwdyT - ViscB * dwdyB) / dz_0;
  }

//Turbulent FGH functions calculation
inline RealType computeTurbulentF2D(
    const RealType* const localVelocity, const RealType* const localViscosity, 
    const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    return localVelocity[mapd(0, 0, 0, 0)]
      + dt * (2*d2udx2(localVelocity, localViscosity, parameters, localMeshsize)
          + d2udy2(localVelocity, localViscosity, parameters, localMeshsize) 
          + d2vdxdy(localVelocity, localViscosity, parameters, localMeshsize)
          - du2dx(localVelocity, parameters, localMeshsize)
          - duvdy(localVelocity, parameters, localMeshsize) 
          + parameters.environment.gx);
  }

inline RealType computeTurbulentG2D(
    const RealType* const localVelocity, const RealType* const localViscosity, 
    const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    return localVelocity[mapd(0, 0, 0, 1)]
      + dt * (d2vdx2(localVelocity, localViscosity, parameters, localMeshsize)
          + d2udydx(localVelocity, localViscosity, parameters, localMeshsize)
          + 2*d2vdy2(localVelocity, localViscosity, parameters, localMeshsize)
          - duvdx(localVelocity, parameters, localMeshsize)
          - dv2dy(localVelocity, parameters, localMeshsize) 
          + parameters.environment.gy);
  }

inline RealType computeTurbulentF3D(
    const RealType* const localVelocity, const RealType* const localViscosity, 
    const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    return localVelocity[mapd(0, 0, 0, 0)]
      + dt * (2*d2udx2(localVelocity, localViscosity, parameters, localMeshsize)
          + d2udy2(localVelocity, localViscosity, parameters, localMeshsize) 
          + d2vdxdy(localVelocity, localViscosity, parameters, localMeshsize)
          + d2udz2(localVelocity, localViscosity, parameters, localMeshsize)
          + d2wdxdz(localVelocity, localViscosity, parameters, localMeshsize)
          - du2dx(localVelocity, parameters, localMeshsize) 
          - duvdy(localVelocity, parameters, localMeshsize)
          - duwdz(localVelocity, parameters, localMeshsize) + parameters.environment.gx);
  }

inline RealType computeTurbulentG3D(
    const RealType* const localVelocity, const RealType* const localViscosity, 
    const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    return localVelocity[mapd(0, 0, 0, 1)]
      + dt * (d2vdx2(localVelocity, localViscosity, parameters, localMeshsize)
          + d2udydx(localVelocity, localViscosity, parameters, localMeshsize)
          + 2*d2vdy2(localVelocity, localViscosity, parameters, localMeshsize)
          + d2vdz2(localVelocity, localViscosity, parameters, localMeshsize)
          + d2wdydz(localVelocity, localViscosity, parameters, localMeshsize)
          - dv2dy(localVelocity, parameters, localMeshsize) 
          - duvdx(localVelocity, parameters, localMeshsize)
          - dvwdz(localVelocity, parameters, localMeshsize) + parameters.environment.gy);
  }

inline RealType computeTurbulentH3D(
    const RealType* const localVelocity, const RealType* const localViscosity, 
    const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    return localVelocity[mapd(0, 0, 0, 2)]
      + dt * (d2wdx2(localVelocity, localViscosity, parameters, localMeshsize)
          + d2udzdx(localVelocity, localViscosity, parameters, localMeshsize)
          + d2wdy2(localVelocity, localViscosity, parameters, localMeshsize)
          + d2vdzdy(localVelocity, localViscosity, parameters, localMeshsize)
          + 2*d2wdz2(localVelocity, localViscosity, parameters, localMeshsize)
          - dw2dz(localVelocity, parameters, localMeshsize) 
          - duwdx(localVelocity, parameters, localMeshsize)
          - dvwdy(localVelocity, parameters, localMeshsize) + parameters.environment.gz);
  }

} // namespace Stencils
