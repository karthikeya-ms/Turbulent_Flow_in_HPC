#pragma once

#include "FlowField.hpp"
#include "Parameters.hpp"

/** Turbulent Flow field
 *
 * Class intended to contain the state of the domain during turbulence.
 */
class TurbulentFlowField : public FlowField {
private:
  ScalarField turb_visc_; //! Scalar field representing the turbulent viscosity
  ScalarField wall_h_; //! Scalar field representing the nearest wall height
  // const ?

public:
  /** Constructor for the 2D turbulent flow field */
  TurbulentFlowField(int Nx, int Ny);

  /** Constructor for the 3D turbulent flow field */
  TurbulentFlowField(int Nx, int Ny, int Nz);

  /** Constructs turbulent flow field from parameters object */
  TurbulentFlowField(const Parameters& parameters);

  virtual ~TurbulentFlowField() = default;

  ScalarField& getTurbVisc(); 
  ScalarField& getWallh(); 

  void getPressureVelocityTurbVisc(RealType& pressure, RealType* const velocity, RealType& turbVisc, int i, int j);
  void getPressureVelocityTurbVisc(RealType& pressure, RealType* const velocity, RealType& turbVisc, int i, int j, int k);
};