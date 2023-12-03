/**
 * Code obtained from iterators.cpph for reference
 * ParallelBoundaryIterator(
    FlowFieldType&                            flowField,
    const Parameters&                         parameters,
    Stencils::BoundaryStencil<FlowFieldType>& stencil,
    int                                       lowOffset  = 0, //concept not clear
    int                                       highOffset = 0
  );
*/

#include "Definitions.hpp"
#include "Parameters.hpp"
#include "Stencils/PressureBufferFillStencil.hpp"
#include "Stencils/PressureBufferReadStencil.hpp"
#include "Stencils/VelocityBufferFillStencil.hpp"
#include "Stencils/VelocityBufferReadStencil.hpp"
#include "Iterators.hpp"
#include "FlowField.hpp"

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

}




