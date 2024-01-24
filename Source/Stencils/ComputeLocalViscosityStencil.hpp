#pragma once

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Stencil to compute the velocity once the pressure has been found.
   */
  class ComputeLocalViscosityStencil: public FieldStencil<TurbulentFlowField> {
    private:
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];
    
    public:
    ComputeLocalViscosityStencil(const Parameters& parameters);
    ~ComputeLocalViscosityStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
