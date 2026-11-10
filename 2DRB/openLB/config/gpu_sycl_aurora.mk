# OpenLB on the Aurora Supercomputer

CXX             := mpicxx
CC              := icx

CXXFLAGS := -O3 -std=c++20
CXXFLAGS += -DOLB_PARALLEL_MODE_MPI_INJECT_ENTROPY

PARALLEL_MODE   := MPI

PLATFORMS       := CPU_SISD GPU_SYCL

SYCL_CXX := icpx
SYCL_CXXFLAGS += -O3 -std=c++20 -fsycl -fsycl-device-code-split=per_kernel
SYCL_CXXFLAGS += -fsycl -fsycl-device-code-split=per_kernel -fsycl-targets=spir64_gen -Xsycl-target-backend "-device pvc"
SYCL_LDFLAGS  += -fsycl -fsycl-device-code-split=per_kernel -fsycl-targets=spir64_gen -Xsycl-target-backend "-device pvc"

SYCL_CXXFLAGS += -DOLB_PARALLEL_MODE_MPI_INJECT_ENTROPY

CXXFLAGS += -DADIOS2_USE_MPI -isystem /opt/aurora/26.26.0/spack/unified/1.1.1/install/linux-x86_64/adios2-2.11.0-vwdc5l7/include 
LDFLAGS += -Wl,-rpath,/opt/aurora/26.26.0/spack/unified/1.1.1/install/linux-x86_64/adios2-2.11.0-vwdc5l7/lib64 /opt/aurora/26.26.0/spack/unified/1.1.1/install/linux-x86_64/adios2-2.11.0-vwdc5l7/lib64/libadios2_cxx_mpi.so.2.11.0 /opt/aurora/26.26.0/spack/unified/1.1.1/install/linux-x86_64/adios2-2.11.0-vwdc5l7/lib64/libadios2_cxx.so.2.11.0 -Wl,-rpath-link,/opt/aurora/26.26.0/spack/unified/1.1.1/install/linux-x86_64/adios2-2.11.0-vwdc5l7/lib64 

FLOATING_POINT_TYPE := float

FEATURES := ADIOS2

USE_EMBEDDED_DEPENDENCIES := ON
