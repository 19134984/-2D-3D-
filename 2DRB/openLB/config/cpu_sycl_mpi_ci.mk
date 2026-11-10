# MPI + GPU_SYCL continuous-integration coverage on the OpenCL CPU runtime.
#
# Builds the SYCL backend (icpx) under PARALLEL_MODE=MPI so the masked
# SendTask/RecvTask gather/scatter path is exercised across ranks; the
# JIT-compiled SPIR-V kernels run on the OpenCL CPU device selected via
# ONEAPI_DEVICE_SELECTOR=opencl:cpu.
#
# Requires OpenMPI's mpicxx wrapping icpx (OMPI_CXX=icpx) — put /usr/bin ahead
# of the oneAPI Intel-MPI wrappers on PATH. The SYCL precompiled unit (icpx)
# needs OpenMPI's headers on its include path (added below) since the mixed
# build compiles it with SYCL_CXX rather than the mpicxx wrapper.
#
# Usage (inside the sycl-ci container):
#  - export PATH=/usr/bin:$PATH OMPI_CXX=icpx OMPI_CC=icx
#  - cp this file to config.mk; make
#  - run under mpirun with ONEAPI_DEVICE_SELECTOR=opencl:cpu

CXX             := mpicxx
CC              := mpicc

CXXFLAGS        := -O3 -std=c++20

PARALLEL_MODE   := MPI

PLATFORMS       := CPU_SISD GPU_SYCL

SYCL_CXX        := icpx
SYCL_CXXFLAGS   += -O3 -std=c++20 -fsycl -fsycl-device-code-split=per_kernel
SYCL_CXXFLAGS   += -isystem /usr/lib/x86_64-linux-gnu/openmpi/include
SYCL_CXXFLAGS   += -isystem /usr/lib/x86_64-linux-gnu/openmpi/include/openmpi
SYCL_LDFLAGS    += -fsycl -fsycl-device-code-split=per_kernel

FLOATING_POINT_TYPE := float

USE_EMBEDDED_DEPENDENCIES := ON
