#pragma once

#include "BoundaryStencil.hpp"
#include "DataStructures.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include <algorithm>
#include <memory>
#include <vector>

namespace Stencils {

  class VelocityBufferReadStencil: public BoundaryStencil<FlowField> {
  public:
    VelocityBufferReadStencil(const Parameters& parameters);
    ~VelocityBufferReadStencil() override = default;

    void applyLeftWall(FlowField& flowField, int i, int j) override;
    void applyRightWall(FlowField& flowField, int i, int j) override;
    void applyBottomWall(FlowField& flowField, int i, int j) override;
    void applyTopWall(FlowField& flowField, int i, int j) override;

    void applyLeftWall(FlowField& flowField, int i, int j, int k) override;
    void applyRightWall(FlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(FlowField& flowField, int i, int j, int k) override;
    void applyTopWall(FlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(FlowField& flowField, int i, int j, int k) override;
    void applyBackWall(FlowField& flowField, int i, int j, int k) override;

    std::unique_ptr<RealType[]> leftVelocityBuffer;
    std::unique_ptr<RealType[]> rightVelocityBuffer;
    std::unique_ptr<RealType[]> topVelocityBuffer;
    std::unique_ptr<RealType[]> bottomVelocityBuffer;
    std::unique_ptr<RealType[]> frontVelocityBuffer;
    std::unique_ptr<RealType[]> backVelocityBuffer;

    const int* localSize;

  };
} // namespace Stencils