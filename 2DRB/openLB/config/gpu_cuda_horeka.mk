CXX             := mpic++
CC              := gcc

CXXFLAGS        := -O3 -Wall -march=native -mtune=native
CXXFLAGS        += -std=c++20

PARALLEL_MODE   := MPI

PLATFORMS       := CPU_SISD GPU_CUDA

CUDA_CXX        := nvcc
CUDA_CXXFLAGS   := -O3 -std=c++20

CUDA_ARCH       := 90

FLOATING_POINT_TYPE := float

USE_EMBEDDED_DEPENDENCIES := ON
