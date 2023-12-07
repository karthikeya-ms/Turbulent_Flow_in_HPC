#include "StdAfx.hpp"

#include "MaxTurbViscStencil.hpp"

Stencils::MaxTurbViscStencil::MaxTurbViscStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {
  reset();
}

void Stencils::MaxTurbViscStencil::apply(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxTurbViscStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxTurbViscStencil::cellMaxValue(TurbulentFlowField& flowField, int i, int j) {
  RealType& turbVisc = flowField.getTurbVisc().getScalar(i, j);
  if (turbVisc > maxValues_) { maxValues_ = turbVisc; }
}

void Stencils::MaxTurbViscStencil::cellMaxValue(TurbulentFlowField& flowField, int i, int j, int k) {
  RealType& turbVisc = flowField.getTurbVisc().getScalar(i, j, k);
  if (turbVisc > maxValues_) { maxValues_ = turbVisc; }
}

void Stencils::MaxTurbViscStencil::reset() { maxValues_ = 0; }

const RealType& Stencils::MaxTurbViscStencil::getMaxValues() const { return maxValues_; }
