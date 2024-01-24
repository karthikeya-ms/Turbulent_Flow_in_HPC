#include "PetscParallelManager.hpp"

#include <mpi.h>
#include <petsclog.h>

#include "FlowField.hpp"
#include "Parameters.hpp"
#include "Stencils/ViscosityBufferFillStencil.hpp"

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


 ParallelManagers::PetscTurbulentParallelManager::PetscTurbulentParallelManager
 (Parameters& parameters, TurbulentFlowField& flowfield):

  PetscParallelManager(parameters, flowfield),
  flowfield_(flowfield),
  fillViscosityStencil(parameters),
  readViscosityStencil(parameters),
  viscosityfillIterator(flowfield, parameters, fillViscosityStencil, 0, 0),
  viscosityreadIterator(flowfield, parameters, readViscosityStencil, 0, 0),
  fillChViscosityStencil(parameters),
  readChViscosityStencil(parameters),
  ChviscosityfillIterator(flowfield, parameters, fillChViscosityStencil, 0, 0),
  ChviscosityreadIterator(flowfield, parameters, readChViscosityStencil, 0, 0) {}

void ParallelManagers::PetscParallelManager::communicatePressure() {

  pressurefillIterator.iterate();

  if (parameters_.geometry.dim == 3) {
    MPI_Request request[12];
    MPI_Status  status[12];
    const int*  localSize = fillPressureStencil.localSize;

    MPI_Isend(
      // left to right
      fillPressureStencil.leftPressureFillBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      101,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readPressureStencil.rightPressureReadBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      101,
      PETSC_COMM_WORLD,
      &request[1]
    );

    MPI_Isend(
      // right to left
      fillPressureStencil.rightPressureFillBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      102,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readPressureStencil.leftPressureReadBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      102,
      PETSC_COMM_WORLD,
      &request[3]
    );
    MPI_Isend(
      // top to bottom
      fillPressureStencil.topPressureFillBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      103,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readPressureStencil.bottomPressureReadBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      103,
      PETSC_COMM_WORLD,
      &request[5]
    );
    MPI_Isend(
      // bottom to top
      fillPressureStencil.bottomPressureFillBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      104,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readPressureStencil.topPressureReadBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      104,
      PETSC_COMM_WORLD,
      &request[7]
    );
    MPI_Isend(
      // front to back
      fillPressureStencil.frontPressureFillBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      105,
      PETSC_COMM_WORLD,
      &request[8]
    );
    MPI_Irecv(
      readPressureStencil.backPressureReadBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      105,
      PETSC_COMM_WORLD,
      &request[9]
    );
    MPI_Isend(
      // back to front
      fillPressureStencil.backPressureFillBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      106,
      PETSC_COMM_WORLD,
      &request[10]
    );
    MPI_Irecv(
      readPressureStencil.frontPressureReadBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      106,
      PETSC_COMM_WORLD,
      &request[11]
    );

    for (size_t i = 0; i < 12; i++) {
      MPI_Wait(&request[i], &status[i]);
    }
  }

  if (parameters_.geometry.dim == 2) {
    MPI_Request request[8];
    MPI_Status  status[8];

    const int* localSize = fillPressureStencil.localSize;

    MPI_Isend(
      // left to right
      fillPressureStencil.leftPressureFillBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      101,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readPressureStencil.rightPressureReadBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      101,
      PETSC_COMM_WORLD,
      &request[1]
    );

    MPI_Isend(
      // right to left
      fillPressureStencil.rightPressureFillBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      102,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readPressureStencil.leftPressureReadBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      102,
      PETSC_COMM_WORLD,
      &request[3]
    );
    MPI_Isend(
      // top to bottom
      fillPressureStencil.topPressureFillBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      103,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readPressureStencil.bottomPressureReadBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      103,
      PETSC_COMM_WORLD,
      &request[5]
    );
    MPI_Isend(
      // bottom to top
      fillPressureStencil.bottomPressureFillBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      104,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readPressureStencil.topPressureReadBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      104,
      PETSC_COMM_WORLD,
      &request[7]
    );
    for (size_t i = 0; i < 8; i++) {
      MPI_Wait(&request[i], &status[i]);
    }
  }
  pressurereadIterator.iterate();
}

void ParallelManagers::PetscParallelManager::communicateVelocity() {
  velocityfillIterator.iterate();
  if (parameters_.geometry.dim == 3) {
    MPI_Request request[12];
    MPI_Status  status[12];

    const int* localSize = fillVelocityStencil.localSize;

    MPI_Isend(
      // left to right
      fillVelocityStencil.leftVelocityFillBuffer.get(),
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      201,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readVelocityStencil.rightVelocityReadBuffer.get(),
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      201,
      PETSC_COMM_WORLD,
      &request[1]
    );

    MPI_Isend(
      // right to left
      fillVelocityStencil.rightVelocityFillBuffer.get(),
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      202,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readVelocityStencil.leftVelocityReadBuffer.get(),
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      202,
      PETSC_COMM_WORLD,
      &request[3]
    );

    MPI_Isend(
      // top to bottom
      fillVelocityStencil.topVelocityFillBuffer.get(),
      (localSize[0] * localSize[2]) + localSize[0] * (localSize[2] + 1) + (localSize[0] + 1) * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      203,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readVelocityStencil.bottomVelocityReadBuffer.get(),
      (localSize[0] * localSize[2]) + localSize[0] * (localSize[2] + 1) + (localSize[0] + 1) * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      203,
      PETSC_COMM_WORLD,
      &request[5]
    );

    MPI_Isend(
      // bottom to top
      fillVelocityStencil.bottomVelocityFillBuffer.get(),
      (localSize[0] * localSize[2]) + localSize[0] * (localSize[2] + 1) + (localSize[0] + 1) * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      204,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readVelocityStencil.topVelocityReadBuffer.get(),
      (localSize[0] * localSize[2]) + localSize[0] * (localSize[2] + 1) + (localSize[0] + 1) * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      204,
      PETSC_COMM_WORLD,
      &request[7]
    );

    MPI_Isend(
      // front to back
      fillVelocityStencil.frontVelocityFillBuffer.get(),
      (localSize[1] * localSize[0]) + localSize[1] * (localSize[0] + 1) + (localSize[1] + 1) * localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      205,
      PETSC_COMM_WORLD,
      &request[8]
    );
    MPI_Irecv(
      readVelocityStencil.backVelocityReadBuffer.get(),
      (localSize[1] * localSize[0]) + localSize[1] * (localSize[0] + 1) + (localSize[1] + 1) * localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      205,
      PETSC_COMM_WORLD,
      &request[9]
    );
    MPI_Isend(
      // back to front
      fillVelocityStencil.backVelocityFillBuffer.get(),
      (localSize[1] * localSize[0]) + localSize[1] * (localSize[0] + 1) + (localSize[1] + 1) * localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      206,
      PETSC_COMM_WORLD,
      &request[10]
    );
    MPI_Irecv(
      readVelocityStencil.frontVelocityReadBuffer.get(),
      (localSize[1] * localSize[0]) + localSize[1] * (localSize[0] + 1) + (localSize[1] + 1) * localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      206,
      PETSC_COMM_WORLD,
      &request[11]
    );

    for (size_t i = 0; i < 12; i++) {
      MPI_Wait(&request[i], &status[i]);
    }
  }

  if (parameters_.geometry.dim == 2) {
    MPI_Request request[8];
    MPI_Status  status[8];

    const int* localSize = fillVelocityStencil.localSize;

    MPI_Isend(
      // left to right
      fillVelocityStencil.leftVelocityFillBuffer.get(),
      (localSize[1]) + (localSize[1] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      201,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readVelocityStencil.rightVelocityReadBuffer.get(),
      (localSize[1]) + (localSize[1] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      201,
      PETSC_COMM_WORLD,
      &request[1]
    );

    MPI_Isend(
      // right to left
      fillVelocityStencil.rightVelocityFillBuffer.get(),
      (localSize[1]) + (localSize[1] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      202,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readVelocityStencil.leftVelocityReadBuffer.get(),
      (localSize[1]) + (localSize[1] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      202,
      PETSC_COMM_WORLD,
      &request[3]
    );

    MPI_Isend(
      // top to bottom
      fillVelocityStencil.topVelocityFillBuffer.get(),
      (localSize[0]) + (localSize[0] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      203,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readVelocityStencil.bottomVelocityReadBuffer.get(),
      (localSize[0]) + (localSize[0] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      203,
      PETSC_COMM_WORLD,
      &request[5]
    );

    MPI_Isend(
      // bottom to top
      fillVelocityStencil.bottomVelocityFillBuffer.get(),
      (localSize[0]) + (localSize[0] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      204,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readVelocityStencil.topVelocityReadBuffer.get(),
      (localSize[0]) + (localSize[0] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      204,
      PETSC_COMM_WORLD,
      &request[7]
    );

    for (size_t i = 0; i < 8; i++) {
      MPI_Wait(&request[i], &status[i]);
    }
  }

  velocityreadIterator.iterate();
}

void ParallelManagers::PetscTurbulentParallelManager::communicateViscosity() {
  viscosityfillIterator.iterate();

  if (parameters_.geometry.dim == 3) {
    MPI_Request request[12];
    MPI_Status  status[12];
    const int*  localSize = fillViscosityStencil.localSize;

    MPI_Isend(
      // left to right
      fillViscosityStencil.leftViscosityFillBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      301,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readViscosityStencil.rightViscosityReadBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      301,
      PETSC_COMM_WORLD,
      &request[1]
    );

    MPI_Isend(
      // right to left
      fillViscosityStencil.rightViscosityFillBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      302,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readViscosityStencil.leftViscosityReadBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      302,
      PETSC_COMM_WORLD,
      &request[3]
    );
    MPI_Isend(
      // top to bottom
      fillViscosityStencil.topViscosityFillBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      303,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readViscosityStencil.bottomViscosityReadBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      303,
      PETSC_COMM_WORLD,
      &request[5]
    );
    MPI_Isend(
      // bottom to top
      fillViscosityStencil.bottomViscosityFillBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      304,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readViscosityStencil.topViscosityReadBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      304,
      PETSC_COMM_WORLD,
      &request[7]
    );
    MPI_Isend(
      // front to back
      fillViscosityStencil.frontViscosityFillBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      305,
      PETSC_COMM_WORLD,
      &request[8]
    );
    MPI_Irecv(
      readViscosityStencil.backViscosityReadBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      305,
      PETSC_COMM_WORLD,
      &request[9]
    );
    MPI_Isend(
      // back to front
      fillViscosityStencil.backViscosityFillBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      306,
      PETSC_COMM_WORLD,
      &request[10]
    );
    MPI_Irecv(
      readViscosityStencil.frontViscosityReadBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      306,
      PETSC_COMM_WORLD,
      &request[11]
    );

    for (size_t i = 0; i < 12; i++) {
      MPI_Wait(&request[i], &status[i]);
    }
  }

  if (parameters_.geometry.dim == 2) {
    MPI_Request request[8];
    MPI_Status  status[8];

    const int* localSize = fillViscosityStencil.localSize;

    MPI_Isend(
      // left to right
      fillViscosityStencil.leftViscosityFillBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      301,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readViscosityStencil.rightViscosityReadBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      301,
      PETSC_COMM_WORLD,
      &request[1]
    );

    MPI_Isend(
      // right to left
      fillViscosityStencil.rightViscosityFillBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      302,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readViscosityStencil.leftViscosityReadBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      302,
      PETSC_COMM_WORLD,
      &request[3]
    );
    MPI_Isend(
      // top to bottom
      fillViscosityStencil.topViscosityFillBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      303,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readViscosityStencil.bottomViscosityReadBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      303,
      PETSC_COMM_WORLD,
      &request[5]
    );
    MPI_Isend(
      // bottom to top
      fillViscosityStencil.bottomViscosityFillBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      304,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readViscosityStencil.topViscosityReadBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT, 
      parameters_.parallel.topNb,
      304,
      PETSC_COMM_WORLD,
      &request[7]
    );
    for (size_t i = 0; i < 8; i++) {
      MPI_Wait(&request[i], &status[i]);
    }
  }
  viscosityreadIterator.iterate();
}



void ParallelManagers::PetscTurbulentParallelManager::communicateChViscosity() {
  ChviscosityfillIterator.iterate();
//****************************************
// For 3D 
//***************************************
  if (parameters_.geometry.dim == 3) {
    // for Ch visc
    MPI_Request Ch_request[12];
    MPI_Status  Ch_status[12];

    const int*  localSize = fillChViscosityStencil.localSize; // localSize same for Ch visc

    //left to right Ch visc (IF CONDITION REQ FOR SA)
    MPI_Isend(
      fillChViscosityStencil.leftChViscosityFillBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      901,
      PETSC_COMM_WORLD,
      &Ch_request[0]
    );
    MPI_Irecv(
      readChViscosityStencil.rightChViscosityReadBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      901,
      PETSC_COMM_WORLD,
      &Ch_request[1]
    );

    // right to left Ch visc
    MPI_Isend(
      fillChViscosityStencil.rightChViscosityFillBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      902,
      PETSC_COMM_WORLD,
      &Ch_request[2]
    );
    MPI_Irecv(
      readChViscosityStencil.leftChViscosityReadBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      902,
      PETSC_COMM_WORLD,
      &Ch_request[3]
    );

    //top to bottom Ch visc
    MPI_Isend(
      fillChViscosityStencil.topChViscosityFillBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      903,
      PETSC_COMM_WORLD,
      &Ch_request[4]
    );
    MPI_Irecv(
      readChViscosityStencil.bottomChViscosityReadBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      903,
      PETSC_COMM_WORLD,
      &Ch_request[5]
    );

    // bottom to top Ch visc
    MPI_Isend(
      fillChViscosityStencil.bottomChViscosityFillBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      904,
      PETSC_COMM_WORLD,
      &Ch_request[6]
    );
    MPI_Irecv(
      readChViscosityStencil.topChViscosityReadBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      904,
      PETSC_COMM_WORLD,
      &Ch_request[7]
    );

    // front to back Ch visc
    MPI_Isend(
      fillChViscosityStencil.frontChViscosityFillBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      905,
      PETSC_COMM_WORLD,
      &Ch_request[8]
    );
    MPI_Irecv(
      readChViscosityStencil.backChViscosityReadBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      905,
      PETSC_COMM_WORLD,
      &Ch_request[9]
    );

    // back to front Ch visc
    MPI_Isend(
      fillChViscosityStencil.backChViscosityFillBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      906,
      PETSC_COMM_WORLD,
      &Ch_request[10]
    );
    MPI_Irecv(
      readChViscosityStencil.frontChViscosityReadBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      906,
      PETSC_COMM_WORLD,
      &Ch_request[11]
    );

    for (size_t i = 0; i < 12; i++) {
      MPI_Wait(&Ch_request[i], &Ch_status[i]);

    }
  }
  //********************************************************
  //for 2D
  //********************************************************
  if (parameters_.geometry.dim == 2) {
    
    // for Ch visc
    MPI_Request Ch_request[8];
    MPI_Status  Ch_status[8];

    const int* localSize = fillChViscosityStencil.localSize;

    // left to right Ch visc
    MPI_Isend(
      fillChViscosityStencil.leftChViscosityFillBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      901,
      PETSC_COMM_WORLD,
      &Ch_request[0]
    );
    MPI_Irecv(
      readChViscosityStencil.rightChViscosityReadBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      901,
      PETSC_COMM_WORLD,
      &Ch_request[1]
    );

    // right to left Ch visc
    MPI_Isend(
      fillChViscosityStencil.rightChViscosityFillBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      902,
      PETSC_COMM_WORLD,
      &Ch_request[2]
    );
    MPI_Irecv(
      readChViscosityStencil.leftChViscosityReadBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      902,
      PETSC_COMM_WORLD,
      &Ch_request[3]
    );

    //top to bottom Ch viscs
    MPI_Isend(
      fillChViscosityStencil.topChViscosityFillBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      903,
      PETSC_COMM_WORLD,
      &Ch_request[4]
    );
    MPI_Irecv(
      readChViscosityStencil.bottomChViscosityReadBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      903,
      PETSC_COMM_WORLD,
      &Ch_request[5]
    );

    // bottom to top Ch Visc
        MPI_Isend(
      fillChViscosityStencil.bottomChViscosityFillBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      904,
      PETSC_COMM_WORLD,
      &Ch_request[6]
    );
    MPI_Irecv(
      readChViscosityStencil.topChViscosityReadBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      904,
      PETSC_COMM_WORLD,
      &Ch_request[7]
    );
    for (size_t i = 0; i < 8; i++) {
      MPI_Wait(&Ch_request[i], &Ch_status[i]);

    }
  }
    ChviscosityreadIterator.iterate();



}

