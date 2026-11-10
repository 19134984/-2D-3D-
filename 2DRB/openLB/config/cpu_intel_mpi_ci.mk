# Build config for Intel oneAPI C++ (icpx) with MPI, CI coverage.
#
# mpiicpx is Intel MPI's icpx wrapper. Runs inside the stock public oneAPI
# hpckit image (intel/oneapi-hpckit).

CXX             := mpiicpx
CC              := icx

CXXFLAGS        := -O3 -Wall
CXXFLAGS        += -std=c++20

PARALLEL_MODE   := MPI

PLATFORMS       := CPU_SISD

USE_EMBEDDED_DEPENDENCIES := ON
