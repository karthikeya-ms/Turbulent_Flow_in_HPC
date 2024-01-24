#pragma once

#include <algorithm>
#include <memory>

#include "BoundaryStencil.hpp"
#include "DataStructures.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {
  
  class ChViscosityBufferFillStencil: public BoundaryStencil<TurbulentFlowField> {
  public:
    
    ChViscosityBufferFillStencil(const Parameters& parameters);
    
    ~ChViscosityBufferFillStencil() override = default;

    

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j) override;
    
    void applyRightWall(TurbulentFlowField& flowField, int i, int j) override;
    

    void applyBottomWall(TurbulentFlowField& flowField, int i, int j) override;
    

    void applyTopWall(TurbulentFlowField& flowField, int i, int j) override;
    

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    

    void applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    
    void applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    
    void applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    
    void applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    
    void applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) override;

    
    std::unique_ptr<RealType[]> leftChViscosityFillBuffer;
    

    std::unique_ptr<RealType[]> rightChViscosityFillBuffer;
    

    std::unique_ptr<RealType[]> topChViscosityFillBuffer;
    

    std::unique_ptr<RealType[]> bottomChViscosityFillBuffer;
    

    std::unique_ptr<RealType[]> frontChViscosityFillBuffer;
    

    std::unique_ptr<RealType[]> backChViscosityFillBuffer;
    
    const int* localSize;
  };
} // namespace Stencils
