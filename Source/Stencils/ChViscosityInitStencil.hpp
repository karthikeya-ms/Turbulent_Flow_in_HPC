#pragma once

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

#include "Assertion.hpp"

namespace Stencils {

  /** Stencil to compute the wall_h_ only once during initialisation.
   */
  class ChViscosityInitStencil: public FieldStencil<TurbulentFlowField>{
    public:
    ChViscosityInitStencil(const Parameters& parameters);
    ~ChViscosityInitStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
