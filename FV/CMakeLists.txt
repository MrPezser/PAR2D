cmake_minimum_required(VERSION 3.25)
project(MPI2DFV)

set(CMAKE_CXX_STANDARD 17)

#include_directories(SYSTEM MSMPI_INC)

#find_package(MPI REQUIRED)
add_executable(MPI2DFV main.cpp EulerFlux.cpp EulerFlux.h BoundaryConditions.cpp BoundaryConditions.h FileIO.cpp FileIO.h Indexing.h SpatialDiscretization.cpp SpatialDiscretization.h MeshModule.cpp MeshModule.h StateVariables.h LUtools.cpp LUtools.h Jacobian.cpp Jacobian.h Thermo.cpp Thermo.h DGP1Tools.cpp DGP1Tools.h)
#target_link_libraries(MPI2DFV PRIVATE MPI::MPI_CXX)