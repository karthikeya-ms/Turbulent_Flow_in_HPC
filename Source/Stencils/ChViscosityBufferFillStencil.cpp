#include "StdAfx.hpp"

#include "ChViscosityBufferFillStencil.hpp"

#include "Definitions.hpp"

Stencils::ChViscosityBufferFillStencil::ChViscosityBufferFillStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters),
  localSize(parameters.parallel.localSize) {

  if (parameters.geometry.dim == 3) {

    leftChViscosityFillBuffer  = std::make_unique<RealType[]>(localSize[1] * localSize[2]);
    rightChViscosityFillBuffer = std::make_unique<RealType[]>(localSize[1] * localSize[2]);

    bottomChViscosityFillBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[2]);
    topChViscosityFillBuffer    = std::make_unique<RealType[]>(localSize[0] * localSize[2]);

    frontChViscosityFillBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[1]);
    backChViscosityFillBuffer  = std::make_unique<RealType[]>(localSize[0] * localSize[1]);

  }

  else {
    leftChViscosityFillBuffer  = std::make_unique<RealType[]>(localSize[1]);
    rightChViscosityFillBuffer = std::make_unique<RealType[]>(localSize[1]);

    bottomChViscosityFillBuffer = std::make_unique<RealType[]>(localSize[0]);
    topChViscosityFillBuffer    = std::make_unique<RealType[]>(localSize[0]);
  }
}
// For 2D Cases
void Stencils::ChViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  if (j >= 2 && j <= (localSize[1] + 1)) {
    *(leftChViscosityFillBuffer.get() + (j - 2)) = (flowField.getNewChVis().getScalar(i + 2, j));
  }
}

void Stencils::ChViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  if (j >= 2 && j <= (localSize[1] + 1)) {
    // Need to verify indices
    *(rightChViscosityFillBuffer.get() + (j - 2)) = (flowField.getNewChVis().getScalar(i - 1, j));
  }
}

void Stencils::ChViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  if ((i >= 2 && i <= localSize[0] + 1)) {
    *(bottomChViscosityFillBuffer.get() + (i - 2)) = (flowField.getNewChVis().getScalar(i, j + 2));
  }
}

void Stencils::ChViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  if ((i >= 2 && i <= localSize[0] + 1)) {
    *(topChViscosityFillBuffer.get() + (i - 2)) = (flowField.getNewChVis().getScalar(i, j - 1));
  }
}

// For 3D Cases
void Stencils::ChViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
    *(leftChViscosityFillBuffer.get() + (j - 2) + (k - 2) * localSize[1]) = (flowField.getNewChVis().getScalar(i + 2, j, k)
    );
  }
}

void Stencils::ChViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
    *(rightChViscosityFillBuffer.get() + (j - 2) + (k - 2) * localSize[1]
    ) = (flowField.getNewChVis().getScalar(i - 1, j, k));
  }
}

void Stencils::ChViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
    *(bottomChViscosityFillBuffer.get() + (i - 2) * localSize[2] + (k - 2)
    ) = (flowField.getNewChVis().getScalar(i, j + 2, k));
  }
}

void Stencils::ChViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
    *(topChViscosityFillBuffer.get() + (i - 2) * localSize[2] + (k - 2)) = (flowField.getNewChVis().getScalar(i, j - 1, k)
    );
  }
}

void Stencils::ChViscosityBufferFillStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
    *(frontChViscosityFillBuffer.get() + (i - 2) * localSize[1] + (j - 2)
    ) = (flowField.getNewChVis().getScalar(i, j, k + 2));
  }
}

void Stencils::ChViscosityBufferFillStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
    *(backChViscosityFillBuffer.get() + (i - 2) * localSize[1] + (j - 2)) = (flowField.getNewChVis().getScalar(i, j, k - 1)
    );
  }
}