#pragma once

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class LaplaceStencil: public FieldStencil<TurbulentFlowField> {
  private:
    RealType cb2 = parameters_.turbSA.cb2;
    RealType cb3 = parameters_.turbSA.cb3;

  public:
    LaplaceStencil(const Parameters& parameters);
    ~LaplaceStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils