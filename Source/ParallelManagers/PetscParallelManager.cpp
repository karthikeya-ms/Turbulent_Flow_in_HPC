#include "StdAfx.hpp"

#include "PetscParallelManager.hpp"

ParallelManagers::PetscParallelManager::PetscParallelManager(
    Parameters &parameters, FlowField &flowField, int lowOffset = 0,
    int highOffset = 0)
    : parameters_(parameters), flowField_(flowField), lowOffset_(lowOffset),
      highOffset_(highOffset) {}

ParallelManagers::PetscParallelManager::~PetscParallelManager() = default;

void ParallelManagers::PetscParallelManager::communicatePressure() {
  auto fillStencil = Stencils::PressureBufferFillStencil(parameters_);
  auto fillIterator = ParallelBoundaryIterator<FlowField>(
      flowField_, parameters_, fillStencil, lowOffset_, highOffset_);
  auto readStencil = Stencils::PressureBufferReadStencil(parameters_);
  auto readIterator = ParallelBoundaryIterator<FlowField>(
      flowField_, parameters_, readStencil, lowOffset_, highOffset_);
  fillIterator.iterate();
  MPI_Status* status;
  MPI_Sendrecv(fillStencil.topPressureFillBuffer, *fillStencil.localSize, MY_MPI_FLOAT, parameters_.parallel.topNb,0, readStencil.topPressureReadBuffer, *readStencil.localSize, MY_MPI_FLOAT, parameters_.parallel.rank, 0, MPI_COMM_WORLD, status);
  readIterator.iterate();
}

void ParallelManagers::PetscParallelManager::communicateVelocities() {
  auto fillStencil = Stencils::VelocityBufferFillStencil(parameters_);
  auto fillIterator = ParallelBoundaryIterator<FlowField>(
      flowField_, parameters_, fillStencil, lowOffset_, highOffset_);
  auto readStencil = Stencils::VelocityBufferReadStencil(parameters_);
  auto readIterator = ParallelBoundaryIterator<FlowField>(
      flowField_, parameters_, readStencil, lowOffset_, highOffset_);
  fillIterator.iterate();
  readIterator.iterate();
}