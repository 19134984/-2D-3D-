# Build config for GPU_SYCL + MPI + ADIOS2 continuous-integration coverage.
#
# Mirrors config/gpu_sycl_aurora.mk but targets the OpenCL CPU device (no fixed
# -fsycl-targets, so JIT-compiled SPIR-V runs on the runtime-selected device).
#
# Used by .ci/build_sycl_windfarm.sh to keep the SC wind-farm application
# (apps/adrian/fsi/windFarm3d) continuously built and smoke-tested in the Intel
# oneAPI hpckit container.
#
# The distribution ADIOS2 in that container is built against OpenMPI, so the
# host MPI must also be OpenMPI for a matching MPI_Comm ABI. The build script
# selects OpenMPI's mpicxx wrapping icpx (OMPI_CXX=icpx) so that both the MPI
# headers and -fsycl reach the compiler. The SYCL device translation unit is
# compiled by icpx directly, so it needs OpenMPI's and ADIOS2's include paths
# in SYCL_CXXFLAGS as well (apps that fold their .cpp into the device TU via
# OLB_MIXED_MODE_INCLUDE_CPP compile their host-side MPI/ADIOS2 code there too).

CXX             := mpicxx
CC              := mpicc

CXXFLAGS        := -O3 -std=c++20
CXXFLAGS        += -DADIOS2_USE_MPI
CXXFLAGS        += -isystem /usr/include/adios2/build/mpi

PARALLEL_MODE   := MPI

PLATFORMS       := CPU_SISD GPU_SYCL

SYCL_CXX        := icpx
SYCL_CXXFLAGS   += -O3 -std=c++20 -fsycl -fsycl-device-code-split=per_kernel
SYCL_CXXFLAGS   += -DFEATURE_ADIOS2 -DADIOS2_USE_MPI -isystem /usr/include/adios2/build/mpi
SYCL_CXXFLAGS   += -isystem /usr/lib/x86_64-linux-gnu/openmpi/include
SYCL_CXXFLAGS   += -isystem /usr/lib/x86_64-linux-gnu/openmpi/include/openmpi
SYCL_LDFLAGS    += -fsycl -fsycl-device-code-split=per_kernel

# Distribution ADIOS2 (OpenMPI-built Debian package in the hpckit-with-adios
# image). The Debian packaging splits the C++ bindings: MPI-specific entry
# points (e.g. ADIOS(MPI_Comm)) live in the *_cxx11_mpi/_core_mpi variants while
# the common bindings (e.g. Engine::Flush) stay in *_cxx11/_core, so all four
# must be listed explicitly under --no-copy-dt-needed-entries.
LDFLAGS         += -ladios2_mpi_cxx11_mpi -ladios2_mpi_cxx11
LDFLAGS         += -ladios2_mpi_core_mpi -ladios2_mpi_core

FLOATING_POINT_TYPE := float

FEATURES        := ADIOS2

USE_EMBEDDED_DEPENDENCIES := ON
