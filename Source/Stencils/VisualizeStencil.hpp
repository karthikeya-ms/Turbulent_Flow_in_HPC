#pragma once

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

#include "Assertion.hpp"

namespace Stencils {

  /** Stencil to compute the wall_h_ only once during initialisation.
   */
  class VisualizeStencil: public FieldStencil<TurbulentFlowField>{
    private: 
    // A local velocity variable that will be used to approximate derivatives. Size matches 3D
    // case, but can be used for 2D as well.
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];
    
    public:
    VisualizeStencil(const Parameters& parameters);
    ~VisualizeStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
