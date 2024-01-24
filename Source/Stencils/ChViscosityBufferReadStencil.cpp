#include "StdAfx.hpp"

#include "ChViscosityBufferReadStencil.hpp"

#include "Definitions.hpp"

Stencils::ChViscosityBufferReadStencil::ChViscosityBufferReadStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters),
  localSize(parameters.parallel.localSize) {

  if (parameters.geometry.dim == 3) {

    leftChViscosityReadBuffer  = std::make_unique<RealType[]>(localSize[1] * localSize[2]);
    rightChViscosityReadBuffer = std::make_unique<RealType[]>(localSize[1] * localSize[2]);

    bottomChViscosityReadBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[2]);
    topChViscosityReadBuffer    = std::make_unique<RealType[]>(localSize[0] * localSize[2]);

    frontChViscosityReadBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[1]);
    backChViscosityReadBuffer  = std::make_unique<RealType[]>(localSize[0] * localSize[1]);

  }

  else {
    leftChViscosityReadBuffer  = std::make_unique<RealType[]>(localSize[1]);
    rightChViscosityReadBuffer = std::make_unique<RealType[]>(localSize[1]);

    bottomChViscosityReadBuffer = std::make_unique<RealType[]>(localSize[0]);
    topChViscosityReadBuffer    = std::make_unique<RealType[]>(localSize[0]);
  }
} // End of constructor

// For 2D Cases
void Stencils::ChViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  if (parameters_.parallel.leftNb >= 0) {
    if (j >= 2 && j <= (localSize[1] + 1)) {
      flowField.getNewChVis().getScalar(i + 1, j) = *(leftChViscosityReadBuffer.get() + (j - 2));
    }
  }
}

void Stencils::ChViscosityBufferReadStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  if (parameters_.parallel.rightNb >= 0) {
    if (j >= 2 && j <= (localSize[1] + 1)) {
      flowField.getNewChVis().getScalar(i, j) = *(rightChViscosityReadBuffer.get() + (j - 2));
    }
  }
}

void Stencils::ChViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  if (parameters_.parallel.bottomNb >= 0) {
    if ((i >= 2 && i <= localSize[0] + 1)) {
      flowField.getNewChVis().getScalar(i, j + 1) = *(bottomChViscosityReadBuffer.get() + (i - 2));
    }
  }
}

void Stencils::ChViscosityBufferReadStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  if (parameters_.parallel.topNb >= 0) {
    if ((i >= 2 && i <= localSize[0] + 1)) {
      flowField.getNewChVis().getScalar(i, j) = *(topChViscosityReadBuffer.get() + (i - 2));
    }
  }
}

// For 3D Cases
void Stencils::ChViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.leftNb >= 0) {
    if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
      flowfield.getNewChVis().getScalar(
        i + 1, j, k
      ) = *(leftChViscosityReadBuffer.get() + (j - 2) + (k - 2) * localSize[1]);
    }
  }
}

void Stencils::ChViscosityBufferReadStencil::applyRightWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.rightNb >= 0) {
    if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
      flowfield.getNewChVis().getScalar(i, j, k) = *(rightChViscosityReadBuffer.get() + (j - 2) + (k - 2) * localSize[1]);
    }
  }
}

void Stencils::ChViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.bottomNb >= 0) {
    if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
      flowfield.getNewChVis().getScalar(
        i, j + 1, k
      ) = *(bottomChViscosityReadBuffer.get() + (k - 2) + (i - 2) * localSize[2]);
    }
  }
}

void Stencils::ChViscosityBufferReadStencil::applyTopWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.topNb >= 0) {
    if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
      flowfield.getNewChVis().getScalar(i, j, k) = *(topChViscosityReadBuffer.get() + (k - 2) + (i - 2) * localSize[2]);
    }
  }
}

void Stencils::ChViscosityBufferReadStencil::applyFrontWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.frontNb >= 0) {
    if ((i >= 2) & (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
      flowfield.getNewChVis().getScalar(
        i, j, k + 1
      ) = *(frontChViscosityReadBuffer.get() + (j - 2) + (i - 2) * localSize[1]);
    }
  }
}

void Stencils::ChViscosityBufferReadStencil::applyBackWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.backNb >= 0) {
    if ((i >= 2) & (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
      flowfield.getNewChVis().getScalar(i, j, k) = *(backChViscosityReadBuffer.get() + (j - 2) + (i - 2) * localSize[1]);
    }
  }
}

