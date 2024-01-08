#include "StdAfx.hpp"

#include "TurbulentFlowField.hpp"

TurbulentFlowField::TurbulentFlowField(int Nx, int Ny):
  FlowField::FlowField(Nx, Ny),
  turb_visc_(ScalarField(Nx + 3, Ny + 3)),
  wall_h_(ScalarField(Nx + 3, Ny + 3)),
  Q_(ScalarField(Nx + 3, Ny + 3)),
  Nabla_(ScalarField(Nx + 3, Ny + 3)),
  ChVis_(ScalarField(Nx + 3, Ny + 3)) {

  ASSERTION(Nx > 0);
  ASSERTION(Ny > 0);
}

TurbulentFlowField::TurbulentFlowField(int Nx, int Ny, int Nz):
  FlowField::FlowField(Nx, Ny, Nz),
  turb_visc_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  wall_h_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  Q_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  Nabla_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  ChVis_(ScalarField(Nx + 3, Ny + 3, Nz + 3)) {
  
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
  ),
  Q_(
    parameters.geometry.dim == 2 ? 
    ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3) : 
    ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3, parameters.parallel.localSize[2] + 3)
  ),
  Nabla_(
    parameters.geometry.dim == 2 ? 
    ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3) : 
    ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3, parameters.parallel.localSize[2] + 3)
  ),
  ChVis_(
    parameters.geometry.dim == 2 ? 
    ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3) : 
    ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3, parameters.parallel.localSize[2] + 3)
  )
  {}

ScalarField& TurbulentFlowField::getTurbVisc() { return turb_visc_; }

ScalarField& TurbulentFlowField::getWallh() { return wall_h_; }

ScalarField& TurbulentFlowField::getQ() { return Q_; }

ScalarField& TurbulentFlowField::getNabla() { return Nabla_; }

ScalarField& TurbulentFlowField::getChVis() { return ChVis_; }

void TurbulentFlowField::getPressureVelocityTurbVisc(RealType& pressure, RealType* const velocity, RealType& turbVisc, int i, int j) {
  getPressureAndVelocity(pressure, velocity, i, j);
  turbVisc = getTurbVisc().getScalar(i, j);
}

void TurbulentFlowField::getPressureVelocityTurbVisc(RealType& pressure, RealType* const velocity, RealType& turbVisc, int i, int j, int k) {
  getPressureAndVelocity(pressure, velocity, i, j, k);
  turbVisc = getTurbVisc().getScalar(i, j, k);
}
