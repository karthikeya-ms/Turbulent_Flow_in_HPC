#include "StdAfx.hpp"

#include "MinGlobalTimeStepStencil.hpp"

Stencils::MinGlobalTimeStepStencil::MinGlobalTimeStepStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {
  reset();
}

void Stencils::MinGlobalTimeStepStencil::apply(TurbulentFlowField& flowField, int i, int j) { dtMaxValue(flowField, i, j); }

void Stencils::MinGlobalTimeStepStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) { dtMaxValue(flowField, i, j, k); }

void Stencils::MinGlobalTimeStepStencil::dtMaxValue(TurbulentFlowField& flowField, int i, int j) {
  RealType localVisc  = 1.0 / parameters_.flow.Re + flowField.getTurbVisc().getScalar(i, j);
  RealType factor     = 1.0 / (parameters_.meshsize->getDx(i, j) * parameters_.meshsize->getDx(i, j)) +
                        1.0 / (parameters_.meshsize->getDy(i, j) * parameters_.meshsize->getDy(i, j));
  RealType local_dt   = 1.0 / (32 * localVisc * factor); 

  if (local_dt < minValue_) { minValue_ = local_dt; }
}

void Stencils::MinGlobalTimeStepStencil::dtMaxValue(TurbulentFlowField& flowField, int i, int j, int k) {
  RealType localVisc  = 1.0 / parameters_.flow.Re + flowField.getTurbVisc().getScalar(i, j, k);
  RealType factor     = 1.0 / (parameters_.meshsize->getDx(i, j, k) * parameters_.meshsize->getDx(i, j, k)) +
                        1.0 / (parameters_.meshsize->getDy(i, j, k) * parameters_.meshsize->getDy(i, j, k)) +
                        1.0 / (parameters_.meshsize->getDz(i, j, k) * parameters_.meshsize->getDz(i, j, k));
  RealType local_dt   = 1.0 / (32 * localVisc * factor);
  
  if (local_dt < minValue_) { minValue_ = local_dt; }
}

void Stencils::MinGlobalTimeStepStencil::reset() { minValue_ = 0; }

RealType Stencils::MinGlobalTimeStepStencil::getMinValue() const { 
  return minValue_ > EPSILON ? minValue_ : EPSILON; }
