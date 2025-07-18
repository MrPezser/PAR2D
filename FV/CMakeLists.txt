cmake_minimum_required(VERSION 3.25)
project(MPI2DFV)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

find_package(MPI REQUIRED)
add_executable(MPI2DFV main.cpp EulerFlux.cpp EulerFlux.h BoundaryConditions.cpp BoundaryConditions.h FileIO.cpp FileIO.h Indexing.h SpatialDiscretization.cpp SpatialDiscretization.h MeshModule.cpp MeshModule.h StateVariables.h LUtools.cpp LUtools.h Jacobian.cpp Jacobian.h Thermo.cpp Thermo.h DGP1Tools.cpp DGP1Tools.h)
target_link_libraries(MPI2DFV PRIVATE MPI::MPI_CXX)
