#pragma once

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

#include "Assertion.hpp"

namespace Stencils {

  /** Stencil to compute the wall_h_ only once during initialisation.
   */
  class WallhStencil: public FieldStencil<TurbulentFlowField>{
    public:
    WallhStencil(const Parameters& parameters);
    ~WallhStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;

    RealType NormalWallDistance(TurbulentFlowField& flowField, int OBSTACLE_DIR, int i, int j, int k=0);
  };

} // namespace Stencils
