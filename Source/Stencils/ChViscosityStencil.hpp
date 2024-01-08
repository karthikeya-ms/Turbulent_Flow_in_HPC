#pragma once

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /*Stencil to compute the characteristic viscosity once the source term Q and Nabla term has been calculated.*/
  class ChViscosityStencil: public FieldStencil<TurbulentFlowField> {
  public:
    ChViscosityStencil(const Parameters& parameters);
    ~ChViscosityStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils