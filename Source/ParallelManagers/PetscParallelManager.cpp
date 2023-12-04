#include "PetscParallelManager.hpp"

#include <mpi.h>
#include <petsclog.h>

#include "FlowField.hpp"
#include "Parameters.hpp"

ParallelManagers::PetscParallelManager::PetscParallelManager(Parameters& parameters, FlowField& flowfield):
  parameters_(parameters),
  flowfield_(flowfield),
  fillVelocityStencil(parameters),
  readVelocityStencil(parameters),
  velocityfillIterator(flowfield, parameters, fillVelocityStencil, 0, 0),
  velocityreadIterator(flowfield, parameters, readVelocityStencil, 0, 0),
  fillPressureStencil(parameters),
  readPressureStencil(parameters),
  pressurefillIterator(flowfield, parameters, fillPressureStencil, 0, 0),
  pressurereadIterator(flowfield, parameters, readPressureStencil, 0, 0) {}

void ParallelManagers::PetscParallelManager::communicatePressure() {

  pressurefillIterator.iterate();

  if (parameters_.geometry.dim == 3) {
    const int*  localSize = fillPressureStencil.localSize;

    /* MPI_Sendrecv( 
      const void* sendbuf , 
      MPI_Count sendcount , 
      MPI_Datatype sendtype , 
      int dest , 
      int sendtag , 
      void* recvbuf , 
      MPI_Count recvcount , 
      MPI_Datatype recvtype , 
      int source , 
      int recvtag , 
      MPI_Comm comm , 
      MPI_Status* status);
    */

    MPI_Sendrecv( 
      //LEFT TO RIGHT
      fillPressureStencil.leftPressureFillBuffer.get(), 
      localSize[1] * localSize[2], 
      MY_MPI_FLOAT, 
      parameters_.parallel.leftNb, 
      201, 
      readPressureStencil.rightPressureReadBuffer.get(), 
      localSize[1] * localSize[2], 
      MY_MPI_FLOAT, 
      parameters_.parallel.rightNb, 
      201, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //RIGHT TO LEFT
      fillPressureStencil.rightPressureFillBuffer.get(), 
      localSize[1] * localSize[2], 
      MY_MPI_FLOAT, 
      parameters_.parallel.rightNb, 
      202, 
      readPressureStencil.leftPressureReadBuffer.get(), 
      localSize[1] * localSize[2], 
      MY_MPI_FLOAT, 
      parameters_.parallel.leftNb, 
      202, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //TOP TO BOTTOM
      fillPressureStencil.topPressureFillBuffer.get(), 
      localSize[0] * localSize[2], 
      MY_MPI_FLOAT, 
      parameters_.parallel.topNb, 
      203, 
      readPressureStencil.bottomPressureReadBuffer.get(), 
      localSize[0] * localSize[2], 
      MY_MPI_FLOAT, 
      parameters_.parallel.bottomNb, 
      203, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //BOTTOM TO TOP
      fillPressureStencil.bottomPressureFillBuffer.get(), 
      localSize[0] * localSize[2], 
      MY_MPI_FLOAT, 
      parameters_.parallel.bottomNb, 
      204, 
      readPressureStencil.topPressureReadBuffer.get(), 
      localSize[0] * localSize[2], 
      MY_MPI_FLOAT, 
      parameters_.parallel.topNb, 
      204, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //FRONT TO BACK
      fillPressureStencil.frontPressureFillBuffer.get(), 
      localSize[0] * localSize[1], 
      MY_MPI_FLOAT, 
      parameters_.parallel.frontNb, 
      205, 
      readPressureStencil.backPressureReadBuffer.get(), 
      localSize[0] * localSize[1], 
      MY_MPI_FLOAT, 
      parameters_.parallel.backNb, 
      205, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //BACK TO FRONT
      fillPressureStencil.backPressureFillBuffer.get(), 
      localSize[0] * localSize[1], 
      MY_MPI_FLOAT, 
      parameters_.parallel.backNb, 
      206, 
      readPressureStencil.frontPressureReadBuffer.get(), 
      localSize[0] * localSize[1], 
      MY_MPI_FLOAT, 
      parameters_.parallel.frontNb, 
      206, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);
    
  }

  if(parameters_.geometry.dim == 2){

    const int* localSize = fillPressureStencil.localSize;

    MPI_Sendrecv( 
      //LEFT TO RIGHT
      fillPressureStencil.leftPressureFillBuffer.get(), 
      localSize[1], 
      MY_MPI_FLOAT, 
      parameters_.parallel.leftNb, 
      201, 
      readPressureStencil.rightPressureReadBuffer.get(), 
      localSize[1], 
      MY_MPI_FLOAT, 
      parameters_.parallel.rightNb, 
      201, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //RIGHT TO LEFT
      fillPressureStencil.rightPressureFillBuffer.get(), 
      localSize[1], 
      MY_MPI_FLOAT, 
      parameters_.parallel.rightNb, 
      202, 
      readPressureStencil.leftPressureReadBuffer.get(), 
      localSize[1], 
      MY_MPI_FLOAT, 
      parameters_.parallel.leftNb, 
      202, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //TOP TO BOTTOM
      fillPressureStencil.topPressureFillBuffer.get(), 
      localSize[0], 
      MY_MPI_FLOAT, 
      parameters_.parallel.topNb, 
      203, 
      readPressureStencil.bottomPressureReadBuffer.get(), 
      localSize[0], 
      MY_MPI_FLOAT, 
      parameters_.parallel.bottomNb, 
      203, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //BOTTOM TO TOP
      fillPressureStencil.bottomPressureFillBuffer.get(), 
      localSize[0], 
      MY_MPI_FLOAT, 
      parameters_.parallel.bottomNb, 
      204, 
      readPressureStencil.topPressureReadBuffer.get(), 
      localSize[0], 
      MY_MPI_FLOAT, 
      parameters_.parallel.topNb, 
      204, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);
  }
  pressurereadIterator.iterate();
}  

void ParallelManagers::PetscParallelManager::communicateVelocity() {

  velocityfillIterator.iterate();

  if (parameters_.geometry.dim == 3) {
    const int*  localSize = fillVelocityStencil.localSize;

    MPI_Sendrecv( 
      //LEFT TO RIGHT
      fillVelocityStencil.leftVelocityFillBuffer.get(),
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.leftNb, 
      101, 
      readVelocityStencil.rightVelocityReadBuffer.get(), 
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.rightNb, 
      101, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //RIGHT TO LEFT
      fillVelocityStencil.rightVelocityFillBuffer.get(), 
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.rightNb, 
      102, 
      readVelocityStencil.leftVelocityReadBuffer.get(), 
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.leftNb, 
      102, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //TOP TO BOTTOM
      fillVelocityStencil.topVelocityFillBuffer.get(), 
      localSize[0] * localSize[2] + (localSize[0] + 1) * localSize[2] + localSize[0] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.topNb, 
      103, 
      readVelocityStencil.bottomVelocityReadBuffer.get(), 
      localSize[0] * localSize[2] + (localSize[0] + 1) * localSize[2] + localSize[0] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.bottomNb, 
      103, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //BOTTOM TO TOP
      fillVelocityStencil.bottomVelocityFillBuffer.get(), 
      localSize[0] * localSize[2] + (localSize[0] + 1) * localSize[2] + localSize[0] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.bottomNb, 
      104, 
      readVelocityStencil.topVelocityReadBuffer.get(), 
      localSize[0] * localSize[2] + (localSize[0] + 1) * localSize[2] + localSize[0] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.topNb, 
      104, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //FRONT TO BACK
      fillVelocityStencil.frontVelocityFillBuffer.get(), 
      localSize[0] * localSize[1] + (localSize[0] + 1) * localSize[1] + localSize[0] * (localSize[1] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.frontNb, 
      105, 
      readVelocityStencil.backVelocityReadBuffer.get(), 
      localSize[0] * localSize[1] + (localSize[0] + 1) * localSize[1] + localSize[0] * (localSize[1] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.backNb, 
      105, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //BACK TO FRONT
      fillVelocityStencil.backVelocityFillBuffer.get(), 
      localSize[0] * localSize[1] + (localSize[0] + 1) * localSize[1] + localSize[0] * (localSize[1] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.backNb, 
      106, 
      readVelocityStencil.frontVelocityReadBuffer.get(), 
      localSize[0] * localSize[1] + (localSize[0] + 1) * localSize[1] + localSize[0] * (localSize[1] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.frontNb, 
      106, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);
    
  }

  if(parameters_.geometry.dim == 2){

    const int* localSize = fillVelocityStencil.localSize;

    MPI_Sendrecv( 
      //LEFT TO RIGHT
      fillVelocityStencil.leftVelocityFillBuffer.get(), 
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.leftNb, 
      101, 
      readVelocityStencil.rightVelocityReadBuffer.get(), 
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.rightNb, 
      101, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //RIGHT TO LEFT
      fillVelocityStencil.rightVelocityFillBuffer.get(), 
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.rightNb, 
      102, 
      readVelocityStencil.leftVelocityReadBuffer.get(), 
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.leftNb, 
      102, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //TOP TO BOTTOM
      fillVelocityStencil.topVelocityFillBuffer.get(), 
      localSize[0] * localSize[2] + (localSize[0] + 1) * localSize[2] + localSize[0] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.topNb, 
      103, 
      readVelocityStencil.bottomVelocityReadBuffer.get(), 
      localSize[0] * localSize[2] + (localSize[0] + 1) * localSize[2] + localSize[0] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.bottomNb, 
      103, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);

    MPI_Sendrecv( 
      //BOTTOM TO TOP
      fillVelocityStencil.bottomVelocityFillBuffer.get(), 
      localSize[0] * localSize[2] + (localSize[0] + 1) * localSize[2] + localSize[0] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.bottomNb, 
      104, 
      readVelocityStencil.topVelocityReadBuffer.get(), 
      localSize[0] * localSize[2] + (localSize[0] + 1) * localSize[2] + localSize[0] * (localSize[2] + 1), 
      MY_MPI_FLOAT, 
      parameters_.parallel.topNb, 
      104, 
      PETSC_COMM_WORLD, 
      MPI_STATUS_IGNORE);
  }
  velocityreadIterator.iterate();


}
