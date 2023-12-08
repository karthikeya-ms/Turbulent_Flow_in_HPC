#include "Definitions.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"
#include "Stencils/PressureBufferFillStencil.hpp"
#include "Stencils/PressureBufferReadStencil.hpp"
#include "Stencils/VelocityBufferFillStencil.hpp"
#include "Stencils/VelocityBufferReadStencil.hpp"
#include "Stencils/ViscosityBufferFillStencil.hpp"
#include "Stencils/ViscosityBufferReadStencil.hpp"
#include "Iterators.hpp"
#include "FlowField.hpp"
#include <limits>
#include <cmath>


namespace ParallelManagers {
    class PetscParallelManager{
    
    protected: 

        Parameters& parameters_;
        FlowField& flowfield_;
        Stencils::VelocityBufferFillStencil fillVelocityStencil;
        Stencils::VelocityBufferReadStencil readVelocityStencil;
        Stencils::PressureBufferFillStencil fillPressureStencil;
        Stencils::PressureBufferReadStencil readPressureStencil;

        ParallelBoundaryIterator<FlowField> velocityfillIterator;
        ParallelBoundaryIterator<FlowField> velocityreadIterator;
        ParallelBoundaryIterator<FlowField> pressurefillIterator;
        ParallelBoundaryIterator<FlowField> pressurereadIterator;

    public: 

        PetscParallelManager(Parameters& parameters, FlowField& flowfield);

        void communicateVelocity();
        void communicatePressure();

        virtual ~PetscParallelManager()= default;



    };

    class PetscTurbulentParallelManager: public PetscParallelManager {
    
    private:
        
        TurbulentFlowField& flowfield_;
        
        Stencils::ViscosityBufferFillStencil fillViscosityStencil;
        
        Stencils::ViscosityBufferReadStencil readViscosityStencil;
        
        ParallelBoundaryIterator<TurbulentFlowField> viscosityfillIterator;
        
        ParallelBoundaryIterator<TurbulentFlowField> viscosityreadIterator;

    public:
        
        PetscTurbulentParallelManager(Parameters& parameters, TurbulentFlowField& flowfield);
        
        void communicateViscosity();
        
        virtual ~PetscTurbulentParallelManager() = default;
  };

}
