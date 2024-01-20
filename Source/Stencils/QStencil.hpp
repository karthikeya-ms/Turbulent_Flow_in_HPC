#pragma once

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class QStencil: public FieldStencil<TurbulentFlowField> {
  private:
    // A local velocity variable that will be used to approximate derivatives. Size matches 3D
    // case, but can be used for 2D as well.
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];

    //Load coefficients for source term of spalart allmaras turbulence model
    RealType cb1 = parameters_.turbSA.cb1;
    RealType ct3 = parameters_.turbSA.ct3;
    RealType ct4 = parameters_.turbSA.ct4;
    RealType cv1 = parameters_.turbSA.cv1;
    RealType cw1 = parameters_.turbSA.cw1;
    RealType cw2 = parameters_.turbSA.cw2;
    RealType cw3 = parameters_.turbSA.cw3;
    RealType _k  = parameters_.turbSA.k;

  public:
    QStencil(const Parameters& parameters);
    ~QStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
