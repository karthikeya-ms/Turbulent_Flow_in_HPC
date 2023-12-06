#include "StdAfx.hpp"

#include "TurbulentFlowField.hpp"

TurbulentFlowField::TurbulentFlowField(int Nx, int Ny):
  FlowField::FlowField(Nx, Ny),
  turb_visc_(ScalarField(Nx + 3, Ny + 3)),
  wall_h_(ScalarField(Nx + 3, Ny + 3)) {

  ASSERTION(Nx > 0);
  ASSERTION(Ny > 0);
}

TurbulentFlowField::TurbulentFlowField(int Nx, int Ny, int Nz):
  FlowField::FlowField(Nx, Ny, Nz),
  turb_visc_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  wall_h_(ScalarField(Nx + 3, Ny + 3, Nz + 3)) {
  
  ASSERTION(Nx > 0);
  ASSERTION(Ny > 0);
  ASSERTION(Nz > 0);
}

TurbulentFlowField::TurbulentFlowField(const Parameters& parameters):
  FlowField::FlowField(parameters),
  turb_visc_(
    parameters.geometry.dim == 2 ? 
    ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3) : 
    ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3, parameters.parallel.localSize[2] + 3)
  ),
  wall_h_(
    parameters.geometry.dim == 2 ? 
    ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3) : 
    ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3, parameters.parallel.localSize[2] + 3)
  )
  {}

ScalarField& TurbulentFlowField::getTurbVisc() { return turb_visc_; }

ScalarField& TurbulentFlowField::getWallh() { return wall_h_; }

void TurbulentFlowField::getPressureVelocityTurbVisc(RealType& pressure, RealType* const velocity, RealType& turbVisc, int i, int j) {
  // RealType* vHere = getVelocity().getVector(i, j);
  // RealType* vLeft = getVelocity().getVector(i - 1, j);
  // RealType* vDown = getVelocity().getVector(i, j - 1);

  // velocity[0] = (vHere[0] + vLeft[0]) / 2;
  // velocity[1] = (vHere[1] + vDown[1]) / 2;

  // pressure = getPressure().getScalar(i, j);

  getPressureAndVelocity(pressure, velocity, i, j);
  turbVisc = getTurbVisc().getScalar(i, j);
  // needs changes as per VTK visualisation (WS2)
}

void TurbulentFlowField::getPressureVelocityTurbVisc(RealType& pressure, RealType* const velocity, RealType& turbVisc, int i, int j, int k) {
  // RealType* vHere = getVelocity().getVector(i, j, k);
  // RealType* vLeft = getVelocity().getVector(i - 1, j, k);
  // RealType* vDown = getVelocity().getVector(i, j - 1, k);
  // RealType* vBack = getVelocity().getVector(i, j, k - 1);

  // velocity[0] = (vHere[0] + vLeft[0]) / 2;
  // velocity[1] = (vHere[1] + vDown[1]) / 2;
  // velocity[2] = (vHere[2] + vBack[2]) / 2;

  // pressure = getPressure().getScalar(i, j, k);

  getPressureAndVelocity(pressure, velocity, i, j, k);
  turbVisc = getTurbVisc().getScalar(i, j, k);
  // needs changes as per VTK visualisation (WS2)
}
