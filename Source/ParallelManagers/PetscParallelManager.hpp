#pragma once

#include "Definitions.hpp"
#include "Iterators.hpp"
#include "Parameters.hpp"
#include "Stencils/PressureBufferFillStencil.hpp"
#include "Stencils/PressureBufferReadStencil.hpp"
#include "Stencils/VelocityBufferFillStencil.hpp"
#include "Stencils/VelocityBufferReadStencil.hpp"

namespace ParallelManagers {

/** Class used to manage communication between processes.
 */
class PetscParallelManager {
private:
  Parameters &parameters_; //! Reference to the parameters
  FlowField &flowField_;
  int lowOffset_;
  int highOffset_;

public:
  PetscParallelManager(Parameters &parameters, FlowField &flowField,
                       int lowOffset, int highOffset);
  ~PetscParallelManager();
  void communicatePressure();
  void communicateVelocities();
};

} // namespace ParallelManagers
