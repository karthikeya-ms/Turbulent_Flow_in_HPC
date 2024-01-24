#pragma once

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class NablaStencil: public FieldStencil<TurbulentFlowField> {
  private:
  //Load coefficients for Nabla term of spalart allmaras turbulence model
    RealType cb2 = parameters_.turbSA.cb2;
    RealType cb3 = parameters_.turbSA.cb3;

  public:
    NablaStencil(const Parameters& parameters);
    ~NablaStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
