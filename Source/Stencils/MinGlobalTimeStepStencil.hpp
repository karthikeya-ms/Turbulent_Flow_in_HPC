#pragma once

#include "BoundaryStencil.hpp"
#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Computes the minimum value of all possible max timestep for all grid cells.
   */
  class MinGlobalTimeStepStencil: public FieldStencil<TurbulentFlowField> {
  private:
    RealType minValue_; //! Stores the maximum module of every component

    /** Sets the maximum possible value of dt for the cell. Replaces minValue_ if less than the current one.
     *
     * 2D version of the function
     * @param flowField Flow field
     * @param i Position in the X direction.
     * @param j Position in the Y direction.
     */
    void dtMaxValue(TurbulentFlowField& flowField, int i, int j);

    /** Sets the maximum possible value of dt for the cell. . Replaces minValue_ if less than the current one.
     *
     * 3D version of the function
     * @param flowField Flow field
     * @param i Position in the X direction.
     * @param j Position in the Y direction.
     * @param k Position in the Z direction.
     */
    void dtMaxValue(TurbulentFlowField& flowField, int i, int j, int k);

  public:
    MinGlobalTimeStepStencil(const Parameters& parameters);
    ~MinGlobalTimeStepStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;

    /** Resets the maximum values to zero before computing the timestep.
     */
    void reset();

    /** Returns the array with the maximum modules of the components of the velocity,
     *  divided by the respective local meshsize.
     */
    RealType getMinValue() const;
  };

} // namespace Stencils
